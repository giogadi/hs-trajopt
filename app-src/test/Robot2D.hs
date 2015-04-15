import qualified Data.Foldable as F
import Test.HUnit
import Test.Framework
import Test.Framework.Providers.HUnit

import qualified Linear.V2 as V
import Numeric.LinearAlgebra.HMatrix

import Data.SQP

--import Data.Collision

import Debug.Trace

startState :: State2D
startState = vector [(-1), (-1)]

goalState :: State2D
goalState = vector [1, 1]

internalPathLength :: Int
internalPathLength = 1

costMatrix :: Int -> Matrix Double
costMatrix numInternalStates
  | numInternalStates < 1 =
      error "costMatrix can't be generated with fewer than 1 internal states!"
  | otherwise =
      let mainDiag = replicate (2 * numInternalStates) 2
          aboveDiag =
            map (\ix -> ((ix, ix + 2), (-1))) $
              take (2 * (numInternalStates - 1)) [0..]
          belowDiag =
            map (\ix -> ((ix + 2, ix), (-1))) $
              take (2 * (numInternalStates - 1)) [0..]
      in  diagl mainDiag +
            assoc (2 * numInternalStates, 2 * numInternalStates) 0
              (aboveDiag ++ belowDiag)

cost :: Vector Double -> Double
cost statesV =
  let states = startState : vectorToPath statesV ++ [goalState]
      squaredDists = zipWith (\x x' -> let dx = x - x' in dot dx dx)
        states $ tail states
  in  0.5 * sum squaredDists

approxCost :: Vector Double -> (Matrix Double, Vector Double, Double)
approxCost _ = ( costMatrix internalPathLength
               , vjoin [negate startState, konst 0.0 (2 * (internalPathLength - 1))] +
                 vjoin [konst 0.0 (2 * (internalPathLength - 1)), negate goalState]
               , 0.5 * (dot startState startState + dot goalState goalState) )

costTest :: Assertion
costTest =
  let pathV = fromList [(-10), 10]
      (costMat, costVec, costAffine) = approxCost pathV
      testApproxCost = 0.5 * (pathV `dot` (costMat #> pathV)) +
        dot pathV costVec + costAffine
  in  (testApproxCost, cost pathV) @?= (400, 400)

type State2D = Vector Double

pathToVector :: [State2D] -> Vector Double
pathToVector = vjoin

vectorToPath :: Vector Double -> [State2D]
vectorToPath = takesV $ replicate internalPathLength 2

obstacleRadius :: Double
obstacleRadius = 0.5

obstaclePosition :: Vector Double
obstaclePosition = vector [0.01, 0.01]

signedDistance :: Vector Double -> (Double, Vector Double)
signedDistance p =
  let dp = p - obstaclePosition
      dpNorm = norm_2 dp
      signedDist = dpNorm - obstacleRadius
  in  (signedDist, dp / scalar dpNorm)

ineqs :: Vector Double -> Vector Double
ineqs statesV =
  let states = vectorToPath statesV
  in  fromList $ map (negate . fst . signedDistance) states

approxIneqs :: Vector Double -> (Matrix Double, Vector Double)
approxIneqs statesV =
  let states = vectorToPath statesV
      (dists, normals) = unzip $ map signedDistance states
      assocList =
        concatMap (\(ix, n) -> [ ((ix, 2 * ix), negate $ n ! 0)
                               , ((ix, 2 * ix + 1), negate $ n ! 1)]) $
          zip [0..] normals
      ineqMat = assoc (internalPathLength, 2 * internalPathLength) 0 assocList
      affine =
        fromList $ zipWith (+) (map negate dists) $ zipWith dot states normals
  in  (ineqMat, affine)

-- robotSquare :: State2D -> Polygon
-- robotSquare state =
--   let [x, y] = toList state
--       points = [V.V2 1 1, V.V2 (-1) 1, V.V2 (-1) (-1), V.V2 1 (-1)]
--   in  map (+ V.V2 x y) points

-- obstacle :: Polygon
-- obstacle = [V.V2 5 5, V.V2 (-5) 5, V.V2 (-5) (-5), V.V2 5 (-5)]

-- ineqs :: Vector Double -> Vector Double
-- ineqs statesV =
--   let states = vectorToPath statesV
--   in  fromList $
--         map (negate . fst . (`signedDistance` obstacle) . robotSquare) states

-- approxIneqs :: Vector Double -> (Matrix Double, Vector Double)
-- approxIneqs statesV =
--   let states = vectorToPath statesV
--       (dists, normals) = unzip $
--         map ((`signedDistance` obstacle) . robotSquare) states
--       assocList = concatMap (\(ix, V.V2 nx ny) -> [ ((ix, 2 * ix), nx)
--                                                 , ((ix, 2 * ix + 1), ny) ]) $
--                     zip [0..] normals
--       ineqMat = assoc (internalPathLength, 2 * internalPathLength) 0 assocList
--       normalsVecs = map (fromList . F.toList) normals
--       affine =
--         negate $ fromList $ zipWith (+) dists $ zipWith dot normalsVecs states
--   in  (ineqMat, affine)

-- ineqsTest :: Assertion
-- ineqsTest =
--   let statesV = pathToVector $ replicate internalPathLength $ vector [(-7), 0]
--       (ineqsMat, ineqsVec) = approxIneqs statesV
--       ineqPenalty = (ineqsMat #> statesV) + ineqsVec
--   in  (ineqPenalty, ineqs statesV) @?=
--         (konst (-1.0) internalPathLength, konst (-1.0) internalPathLength)

planningTest :: Assertion
planningTest =
  let problem = Problem
                { _cost = cost
                , _approxCost = approxCost
                , _trueIneqs = ineqs
                , _approxAffineIneqs = approxIneqs
                , _numVariables = 2 * internalPathLength
                , _numIneqs = internalPathLength
                }
      initPath =
        map (\t -> startState + scalar (fromIntegral t / fromIntegral (internalPathLength + 1)) *
              (goalState - startState))
          [1..internalPathLength]
      (xResult, fResult) = optimize problem $ pathToVector initPath
      (xTrue, fTrue) = ((-50), (-50)) -- whatever
  in  (xResult, fResult) @?= (xTrue, fTrue)

main = defaultMain tests
  where tests = [ -- testCase "cost-test" costTest
                -- , testCase "ineq-test" ineqsTest
                 testCase "planning-test" planningTest ]
