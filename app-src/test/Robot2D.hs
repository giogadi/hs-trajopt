{-# LANGUAGE TypeFamilies #-}

import qualified Data.Foldable as F
import Test.HUnit
import Test.Framework
import Test.Framework.Providers.HUnit

import qualified Linear.V2 as V
import Numeric.AD
import Numeric.LinearAlgebra.HMatrix

import Data.SQP

import Debug.Trace

startState :: State2D
startState = vector [(-1), (-1)]

goalState :: State2D
goalState = vector [1, 1]

internalPathLength :: Int
internalPathLength = 10

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
      in  diagBlock [diagl mainDiag +
            assoc (2 * numInternalStates, 2 * numInternalStates) 0
              (aboveDiag ++ belowDiag), matrix 1 [0]]

cost :: Vector Double -> Double
cost pathV =
  let (internalStates, _) = vectorToPath pathV
      states = startState : internalStates ++ [goalState]
      squaredDists = zipWith (\x x' -> let dx = x - x' in dot dx dx)
        states $ tail states
  in  0.5 * sum squaredDists

approxCost :: Vector Double -> (Matrix Double, Vector Double, Double)
approxCost _ = ( costMatrix internalPathLength
               , vjoin [(vjoin [negate startState, konst 0.0 (2 * (internalPathLength - 1))] +
                 vjoin [konst 0.0 (2 * (internalPathLength - 1)), negate goalState]), vector [0]]
               , 0.5 * (dot startState startState + dot goalState goalState) )

costTest :: Assertion
costTest =
  let pathV = fromList [(-10), 10]
      (costMat, costVec, costAffine) = approxCost pathV
      testApproxCost = 0.5 * (pathV `dot` (costMat #> pathV)) +
        dot pathV costVec + costAffine
  in  (testApproxCost, cost pathV) @?= (400, 400)

type State2D = Vector Double

pathToVector :: [State2D] -> Double -> Vector Double
pathToVector states stepSize = vjoin $ states ++ [vector [stepSize]]

vectorToPath :: Vector Double -> ([State2D], Double)
vectorToPath pathV = ( takesV (replicate internalPathLength 2) pathV
                     , pathV ! (2 * internalPathLength) )

obstacleRadius :: Double
obstacleRadius = 0.5

obstaclePosition :: Vector Double
obstaclePosition = vector [0.001, 0.001]

signedDistance :: State2D -> (Double, Vector Double)
signedDistance p =
  let dp = p - obstaclePosition
      dpNorm = norm_2 dp
      signedDist = dpNorm - obstacleRadius
  in  (signedDist, dp / scalar dpNorm)

ineqs :: Vector Double -> Vector Double
ineqs pathV =
  let (states, stepSize) = vectorToPath pathV
      collisionIneqs = map (negate . fst . signedDistance) states
      boundIneqs = negate stepSize
  in  fromList $ collisionIneqs ++ [boundIneqs]

approxIneqs :: Vector Double -> (Matrix Double, Vector Double)
approxIneqs pathV =
  let (states, _) = vectorToPath pathV
      (dists, normals) = unzip $ map signedDistance states
      assocList =
        concatMap (\(ix, n) -> [ ((ix, 2 * ix), negate $ n ! 0)
                               , ((ix, 2 * ix + 1), negate $ n ! 1)]) $
          zip [0..] normals
      ineqMat = diagBlock
        [ assoc (internalPathLength, 2 * internalPathLength) 0 assocList
        , matrix 1 [(-1)] ]
      affine =
        fromList $ (zipWith (+) (map negate dists) $ zipWith dot states normals)
          ++ [0]
  in  (ineqMat, affine)

unconcat :: [Int] -> [a] -> [[a]]
unconcat [] [] = []
unconcat _ [] = []
unconcat (n : ns) xs = let (xh, xt) = splitAt n xs
                       in  xh : unconcat ns xt

eqsL :: (Mode a, Scalar a ~ Double) => [a] -> [a]
eqsL pathL =
  let stepSize = last pathL
      internalStatesFlat = take (2 * internalPathLength) pathL
      internalStates = unconcat (repeat 2) internalStatesFlat
      states = (map auto $ toList startState) : internalStates ++ [map auto $ toList goalState]
      dist2 s s' = sum $ map (^ (2 :: Int)) $ zipWith (-) s s'
      squaredDists = zipWith dist2 states $ tail states
      stepSizeSqrd = stepSize * stepSize
  in  map (\d2 -> d2 - stepSizeSqrd) squaredDists

eqs :: Vector Double -> Vector Double
eqs pathV =
  vector $ eqsL $ toList pathV

approxEqs :: Vector Double -> (Matrix Double, Vector Double)
approxEqs pathV =
  let f = fromList $ eqsL $ toList pathV
      f' = fromLists $ jacobian eqsL $ toList pathV
  in  ( f', f - (f' #> pathV) )

planningTest :: Assertion
planningTest =
  let problem = Problem
                { _cost = cost
                , _approxCost = approxCost
                , _trueIneqs = ineqs
                , _approxAffineIneqs = approxIneqs
                , _eqs = eqs
                , _approxEqs = approxEqs
                , _numVariables = 2 * internalPathLength + 1
                , _numIneqs = internalPathLength + 1
                , _numEqs = internalPathLength + 1
                }
      initPathStates =
        map (\t -> startState + scalar (fromIntegral t / fromIntegral (internalPathLength + 1)) *
              (goalState - startState))
          [1..internalPathLength]
      initPathStatesPerturbed = map (+ vector [0.0, 0.0001]) initPathStates
      shortestPathLength = norm_2 $ goalState - startState
      initPathVec = pathToVector initPathStatesPerturbed $
       shortestPathLength / fromIntegral (internalPathLength + 1)
      (xResult, fResult) = optimize problem initPathVec
      (xTrue, fTrue) = ((-50), (-50)) -- whatever
  in  (xResult, fResult) @?= (xTrue, fTrue)

main = defaultMain tests
  where tests = [ -- testCase "cost-test" costTest
                -- , testCase "ineq-test" ineqsTest
                 testCase "planning-test" planningTest ]
