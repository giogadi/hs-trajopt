module Data.SQP
       ( Problem (..)
       , optimize
       ) where

import Numeric.LinearAlgebra.HMatrix
import Numeric.Minimization.QuadProgPP

import Debug.Trace

data Problem = Problem
               { _cost :: Vector Double -> Double
               , _approxCost :: Vector Double
                             -> (Matrix Double, Vector Double, Double)
               , _trueIneqs :: Vector Double -> Vector Double
               , _approxAffineIneqs :: Vector Double
                                    -> (Matrix Double, Vector Double)

               , _numVariables :: Int
               , _numIneqs :: Int
               }

stepAcceptanceThreshold :: Double
stepAcceptanceThreshold = 0.25

trustShrinkFactor :: Double
trustShrinkFactor = 0.1

trustExpandFactor :: Double
trustExpandFactor = 1.5

constraintPenaltyScalingFactor :: Double
constraintPenaltyScalingFactor = 10

minTrustSize :: Double
minTrustSize = 1e-4

minModelImprove :: Double
minModelImprove = 1e-4

minModelImproveRatio :: Double
minModelImproveRatio = negate $ read "Infinity"

constraintSatisfactionThreshold :: Double
constraintSatisfactionThreshold = 1e-4

initTrustSize :: Double
initTrustSize = 0.1

initPenalty :: Double
initPenalty = 10

evalMerit :: Problem -> Double -> Vector Double -> Double
evalMerit problem penaltyParam x =
  let cost = _cost problem x
      penalty =
        penaltyParam * sum (map (max 0.0) $ toList $ _trueIneqs problem x)
  in  cost + penalty

optimize :: Problem
         -> Vector Double -- xInitial
         -> (Vector Double, Double)
optimize problem xInit =
  let meritInit = evalMerit problem initPenalty xInit
  in  findSuitableConstraintPenalty
        problem xInit meritInit initTrustSize initPenalty

findSuitableConstraintPenalty :: Problem
                              -> Vector Double -- x
                              -> Double -- old merit
                              -> Double -- trust region size
                              -> Double -- constraint penalty
                              -- xNew, trueMerit
                              -> (Vector Double, Double)
findSuitableConstraintPenalty problem x merit trustSize penaltyParam =
  let (xNew, trueMerit, newTrustSize) =
        reconvexify problem x merit trustSize penaltyParam
      constraintsSatisfied = all (<= constraintSatisfactionThreshold) $
                               toList $ _trueIneqs problem xNew
  in  if constraintsSatisfied
      then (xNew, trueMerit)
      else let penalty' = penaltyParam * constraintPenaltyScalingFactor
               updatedMerit = evalMerit problem penalty' xNew
           in  findSuitableConstraintPenalty
                 problem xNew updatedMerit newTrustSize penalty'

reconvexify :: Problem
            -> Vector Double -- x
            -> Double -- old true merit
            -> Double -- trust region size
            -> Double -- constraint penalty
            -- xNew, newTrueMerit, newTrustSize
            -> (Vector Double, Double, Double)
reconvexify problem x oldTrueMerit trustSize penaltyParam =
  let convexCost = _approxCost problem x
      (convexIneqMat, convexIneqVec) = _approxAffineIneqs problem x
      trustResult = findSuitableTrustStep problem x convexCost
                      (convexIneqMat, convexIneqVec) oldTrueMerit trustSize
                      penaltyParam
  in  case trustResult of
        Reconvexify xNew newMerit newTrustSize ->
          reconvexify problem xNew newMerit newTrustSize penaltyParam
        Finished xNew newMerit newTrustSize ->
          (xNew, newMerit, newTrustSize)

data TrustStepResult = -- x merit trustSize
                       Reconvexify (Vector Double) Double Double
                     | Finished (Vector Double) Double Double
  deriving (Show)

findSuitableTrustStep :: Problem
                      -> Vector Double -- x
                      -- convexified cost
                      -> (Matrix Double, Vector Double, Double)
                      -- convexified inequality constraints
                      -> (Matrix Double, Vector Double)
                      -> Double -- old true merit
                      -> Double -- trust region size
                      -> Double -- constraint penalty parameter
                      -- xNew, newTrueMerit, newTrustSize
                      -> TrustStepResult
findSuitableTrustStep
  problem x convexCost convexIneq oldTrueMerit trustSize penaltyParam =
  let trustStep =
        trustRegionStep problem x convexCost convexIneq oldTrueMerit
          trustSize penaltyParam
  in  case trustStep of
        Reject -> let newTrustSize = trustSize * trustShrinkFactor
                  in  findSuitableTrustStep problem x convexCost convexIneq
                        oldTrueMerit newTrustSize penaltyParam
        Accept xNew newTrueMerit ->
          let newTrustSize = trustSize * trustExpandFactor
          in  Reconvexify xNew newTrueMerit newTrustSize
        Converged xNew newTrueMerit ->
          let newTrustSize =
                max trustSize $ (minTrustSize / trustShrinkFactor) * 1.5
          in  Finished xNew newTrueMerit newTrustSize

data IterationResult = Reject
                     -- xNew newMerit
                     | Accept (Vector Double) Double
                     | Converged (Vector Double) Double
  deriving (Show)

trustRegionStep :: Problem
                -> Vector Double -- x
                -- convexified cost
                -> (Matrix Double, Vector Double, Double)
                -- convexified inequality constraints
                -> (Matrix Double, Vector Double)
                -> Double -- old true merit
                -> Double -- trust region size
                -> Double -- constraint penalty parameter
                -> IterationResult
trustRegionStep
  problem x convexCost convexIneq oldTrueMerit trustSize penaltyParam =
    let (xNew, modelMerit) = solveQuadraticSubproblem problem x
                              convexCost convexIneq trustSize penaltyParam
        trueMerit = evalMerit problem penaltyParam xNew
        trueImprove = oldTrueMerit - trueMerit
        modelImprove = oldTrueMerit - modelMerit
    in  if modelImprove < -1e-6
        then error $ "Model improvement got worse: " ++
               show (modelImprove, modelMerit, oldTrueMerit, trueMerit)
        else
          if trustSize < minTrustSize ||
             modelImprove < minModelImprove ||
             modelImprove / oldTrueMerit < minModelImproveRatio
          then Converged x oldTrueMerit
          else
            if trueImprove > 0.0 ||
               trueImprove / modelImprove > stepAcceptanceThreshold
            then Accept xNew trueMerit
            else Reject

solveQuadraticSubproblem :: Problem
                         -> Vector Double
                         -- convex cost
                         -> (Matrix Double, Vector Double, Double)
                         -- convexified inequality constraints
                         -> (Matrix Double, Vector Double)
                         -> Double -- trust region size
                         -> Double -- constraint penalty parameter
                         -> (Vector Double, Double) -- new x, model merit
solveQuadraticSubproblem
  problem x (costMatrix, costVector, costConstant)
  (approxIneqMatrix, approxIneqVector) trustSize penaltyParam =
  -- We approximate each nonlinear inequality as |ax + b|^+. For each
  -- of these, we introduce a new optimization variable t (i.e., a
  -- slack variable) that comes with two inequalities:
  --
  -- 0 <= t
  --
  -- ax + b <= t
  --
  let numIneqs = _numIneqs problem
      numVariables = _numVariables problem

      costMatrixWithSlacks =
        diagBlock [ costMatrix
                  , konst 0.0 (numIneqs, numIneqs)]

      -- LET'S MAKE COST SPD BICTH
      (l, v) = eigSH costMatrixWithSlacks
      positiveEigVals = cmap (max 1e-12) l
      augmentedCostMatrix = v <> diag positiveEigVals <> tr v

      augmentedCostVector = vjoin [costVector, konst 1.0 numIneqs]

      augmentedIneqMatrix =
        fromBlocks [ [ scalar penaltyParam * approxIneqMatrix
                     , negate (ident numIneqs)]
                   , [0, negate (ident numIneqs)]]
      augmentedIneqVector =
        vjoin [scalar penaltyParam * approxIneqVector, konst 0.0 numIneqs]

      ineqMatrixWithTrustConstraints =
        augmentedIneqMatrix ===
        (ident numVariables ||| konst 0.0 (numVariables, numIneqs)) ===
        (negate (ident numVariables) ||| konst 0.0 (numVariables, numIneqs))
      ineqVectorWithTrustConstraints =
        vjoin [ augmentedIneqVector
              , negate (x + scalar trustSize)
              , x - scalar trustSize]

      -- inequalities are given as ax + b <= 0, but quadprog++ wants
      -- a'x + b' >= 0
      result = solveQuadProg
        (augmentedCostMatrix, augmentedCostVector)
        Nothing $
        Just ((-ineqMatrixWithTrustConstraints),
              (-ineqVectorWithTrustConstraints))

  in  case result of
        Left e -> error $ show e
        Right (xNewAugmented, newMerit) ->
          let xNew = subVector 0 numVariables xNewAugmented

              -- Compute the merit of the approximation without regularization
              cost = xNew `dot` ((0.5 * (costMatrix #> xNew)) + costVector)
              approxConstraints =
                (approxIneqMatrix #> xNew) + approxIneqVector
              penalty =
                penaltyParam * (sum $ map (max 0.0) $ toList approxConstraints)
              meritNoReg = cost + penalty

              -- DEBUG
              meritError = abs $ newMerit - meritNoReg
          in  if meritError / (abs meritNoReg) > 0.01
              then error $ "merits don't match!" ++ show (meritNoReg, newMerit)
              else (xNew, meritNoReg + costConstant)
