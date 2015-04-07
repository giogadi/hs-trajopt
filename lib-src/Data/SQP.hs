module Data.SQP
       ( Problem

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
stepAcceptanceThreshold = undefined

trustShrinkFactor :: Double
trustShrinkFactor = undefined

trustExpandFactor :: Double
trustExpandFactor = undefined

initialTrust :: Double
initialTrust = undefined

xRelativeTolerance :: Double
xRelativeTolerance = undefined

costRelativeTolerance :: Double
costRelativeTolerance = undefined

constraintPenaltyScalingFactor :: Double
constraintPenaltyScalingFactor = undefined

findSuitableConstraintPenalty :: Problem
                              -> Vector Double -- x
                              -> Double -- trust region size
                              -> Double -- constraint penalty
                              -- xNew, trueCost
                              -> (Vector Double, Double)
findSuitableConstraintPenalty problem x trustSize penaltyParam =
  let (xNew, trueCost, newTrustSize) =
        reconvexify problem x trustSize penaltyParam
      constraintsSatisfied = all (<= 0.0) $ toList $ _trueIneqs problem xNew
  in  if constraintsSatisfied
      then (xNew, trueCost)
      else let newPenalty = penaltyParam * constraintPenaltyScalingFactor
           in  findSuitableConstraintPenalty problem xNew newTrustSize

reconvexify :: Problem
            -> Vector Double -- x
            -> Double -- old true cost
            -> Double -- old model cost
            -> Double -- trust region size
            -> Double -- constraint penalty
            -- xNew, newTrueCost, newModelCost, newTrustSize
            -> (Vector Double, Double, Double, Double)
reconvexify problem x oldTrueCost oldModelCost trustSize penaltyParam =
  let (convexIneqMat, convexIneqVec) = _approxAffineIneqs problem x
      (xNew, newTrueCost, newModelCost, newTrustSize) =
        findSuitableTrustStep problem x (convexIneqMat, convexIneqVec)
          oldTrueCost oldModelCost trustSize penaltyParam
      xDiff = xNew - x
      costDiff = newTrueCost - oldTrueCost
  in  if norm_2 xDiff < xRelativeTolerance * (1 + norm_2 x) || -- xtol
         norm_2 costDiff < costRelativeTolerance * (1 + norm_2 oldTrueCost) -- ftol
      then (xNew, newTrueCost, newModelCost, newTrustSize)
      else reconvexify problem
             xNew newTrueCost newModelCost newTrustSize penaltyParam


findSuitableTrustStep :: Problem
                      -> Vector Double -- x
                      -- convexified inequality constraints
                      -> (Matrix Double, Vector Double)
                      -> Double -- old true cost
                      -> Double -- old model cost
                      -> Double -- trust region size
                      -> Double -- constraint penalty parameter
                      -- xNew, newTrueCost, newModelCost, newTrustSize
                      -> (Vector Double, Double, Double, Double)
findSuitableTrustStep
  problem x (ineqMat, ineqVec) oldTrueCost oldModelCost trustSize penaltyParam =
  let trustStep =
        trustRegionStep problem x (ineqMat, ineqVec) oldTrueCost oldModelCost
          trustSize penaltyParam
  in  case trustStep of
        -- TODO handle case where s < xtol
        Reject -> let newTrustSize = trustSize * trustShrinkFactor
                  in  trustRegionOptimize
                        problem x oldTrueCost oldModelCost
                        newTrustSize penaltyParam
        Accept (xNew, newTrueCost, newModelCost) ->
          let newTrustSize = trustSize * trustExpandFactor
          in  (xNew, newTrueCost, newModelCost, newTrustSize)

data IterationResult = Reject
                     -- xNew, newTrueCost, newModelCost
                     | Accept (Vector Double, Double, Double)

trustRegionStep :: Problem
                -> Vector Double -- x
                -- convexified inequality constraints
                -> (Matrix Double, Vector Double)
                -> Double -- old true cost
                -> Double -- old model cost
                -> Double -- trust region size
                -> Double -- constraint penalty parameter
                -> IterationResult
trustRegionStep
  problem x (ineqMat, ineqVec) oldTrueCost oldModelCost trustSize penaltyParam =
  let (xNew, modelCost) = solveQuadraticSubproblem
                            problem x (ineqMat, ineqVec) trustSize penaltyParam
      (costMatrix, costVector) = _cost problem
      -- cost(x) = x'Ax + x'b = x'(Ax + b)
      trueCost = x `dot` ((costMatrix #> x) + costVector)
      trueImprove = oldTrueCost - trueImprove
      modelImprove = modelCost - oldModelCost
  in  if trueImprove / modelImprove > stepAcceptanceThreshold
      then Accept (xNew, trueCost, modelCost)
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
