module Data.SQP
       ( Problem
       , optimize
       ) where

import Numeric.LinearAlgebra.HMatrix
import Numeric.Minimization.QuadProgPP

data Problem = Problem
               { _cost :: (Matrix Double, Vector Double)
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

evalCost :: Problem -> Vector Double -> Double
evalCost problem x =
  let (costMatrix, costVector) = _cost problem
      -- cost(x) = x'Ax + x'b = x'(Ax + b)
  in  x `dot` ((costMatrix #> x) + costVector)

optimize :: Problem
         -> Vector Double -- xInitial
         -> (Vector Double, Double)
optimize problem xInit =
  let costInit = evalCost problem xInit
  in  findSuitableConstraintPenalty
        problem xInit costInit initTrustSize initPenalty

findSuitableConstraintPenalty :: Problem
                              -> Vector Double -- x
                              -> Double -- old cost
                              -> Double -- trust region size
                              -> Double -- constraint penalty
                              -- xNew, trueCost
                              -> (Vector Double, Double)
findSuitableConstraintPenalty problem x cost trustSize penaltyParam =
  let (xNew, trueCost, newTrustSize) =
        reconvexify problem x cost trustSize penaltyParam
      constraintsSatisfied = all (<= constraintSatisfactionThreshold) $
                               toList $ _trueIneqs problem xNew
  in  if constraintsSatisfied
      then (xNew, trueCost)
      else let penalty' = penaltyParam * constraintPenaltyScalingFactor
           in  findSuitableConstraintPenalty
                 problem xNew cost newTrustSize penalty'

reconvexify :: Problem
            -> Vector Double -- x
            -> Double -- old true cost
            -> Double -- trust region size
            -> Double -- constraint penalty
            -- xNew, newTrueCost, newTrustSize
            -> (Vector Double, Double, Double)
reconvexify problem x oldTrueCost trustSize penaltyParam =
  let (convexIneqMat, convexIneqVec) = _approxAffineIneqs problem x
      trustResult = findSuitableTrustStep problem x
                      (convexIneqMat, convexIneqVec) oldTrueCost trustSize
                      penaltyParam
  in  case trustResult of
        Reconvexify xNew newCost newTrustSize ->
          reconvexify problem xNew newCost newTrustSize penaltyParam
        Finished xNew newCost newTrustSize ->
          (xNew, newCost, newTrustSize)

data TrustStepResult = -- x cost trustSize
                       Reconvexify (Vector Double) Double Double
                     | Finished (Vector Double) Double Double

findSuitableTrustStep :: Problem
                      -> Vector Double -- x
                      -- convexified inequality constraints
                      -> (Matrix Double, Vector Double)
                      -> Double -- old true cost
                      -> Double -- trust region size
                      -> Double -- constraint penalty parameter
                      -- xNew, newTrueCost, newTrustSize
                      -> TrustStepResult
findSuitableTrustStep
  problem x (ineqMat, ineqVec) oldTrueCost trustSize penaltyParam =
  let trustStep =
        trustRegionStep problem x (ineqMat, ineqVec) oldTrueCost
          trustSize penaltyParam
  in  case trustStep of
        Reject -> let newTrustSize = trustSize * trustShrinkFactor
                  in  findSuitableTrustStep problem x (ineqMat, ineqVec)
                        oldTrueCost newTrustSize penaltyParam
        Accept xNew newTrueCost ->
          let newTrustSize = trustSize * trustExpandFactor
          in  Reconvexify xNew newTrueCost newTrustSize
        Converged xNew newTrueCost ->
          let newTrustSize =
                max trustSize $ (minTrustSize / trustShrinkFactor) * 1.5
          in  Finished xNew newTrueCost newTrustSize

data IterationResult = Reject
                     -- xNew newCost
                     | Accept (Vector Double) Double
                     | Converged (Vector Double) Double

trustRegionStep :: Problem
                -> Vector Double -- x
                -- convexified inequality constraints
                -> (Matrix Double, Vector Double)
                -> Double -- old true cost
                -> Double -- trust region size
                -> Double -- constraint penalty parameter
                -> IterationResult
trustRegionStep
  problem x (ineqMat, ineqVec) oldTrueCost trustSize penaltyParam =
    let (xNew, modelCost) = solveQuadraticSubproblem problem x
                              (ineqMat, ineqVec) trustSize penaltyParam
        trueCost = evalCost problem x
        trueImprove = oldTrueCost - trueImprove
        modelImprove = modelCost - oldTrueCost
    in  if modelImprove < 0.0
        then error $ "Model improvement got worse: " ++ show modelImprove
        else
          if trustSize < minTrustSize ||
             modelImprove < minModelImprove ||
             modelImprove / oldTrueCost < minModelImproveRatio
          then Converged x oldTrueCost
          else
            if trueImprove < 0.0 ||
               trueImprove / modelImprove > stepAcceptanceThreshold
            then Accept xNew trueCost
            else Reject

solveQuadraticSubproblem :: Problem
                         -> Vector Double
                         -- convexified inequality constraints
                         -> (Matrix Double, Vector Double)
                         -> Double -- trust region size
                         -> Double -- constraint penalty parameter
                         -> (Vector Double, Double) -- new x, model cost
solveQuadraticSubproblem
  problem x (approxIneqMatrix, approxIneqVector) trustSize penaltyParam =
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

      (costMatrix, costVector) = _cost problem
      augmentedCostMatrix =
        costMatrix ||| konst 0.0 (rows costMatrix, numIneqs)
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
        Right (xNewAugmented, xCost) ->
          (subVector 0 numVariables xNewAugmented, xCost)
