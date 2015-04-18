import Test.HUnit
import Test.Framework
import Test.Framework.Providers.HUnit

import Numeric.LinearAlgebra.HMatrix
import Numeric.Minimization.QuadProgPP

import Numeric.AD

import Data.SQP

qpTest :: Assertion
qpTest =
  -- From example here:
  -- http://www.mathworks.com/help/optim/ug/quadprog.html#zmw57dd0e69051
  let a = matrix 2 [1, (-1), (-1), 2]
      b = vector [(-2), (-6)]
      e = matrix 2 [ (-1), (-1)
                   , 1, (-2)
                   , (-2), (-1)
                   , 1, 0
                   , 0, 1 ]
      f = vector [2, 2, 3, 0, 0]
      result = solveQuadProg (a, b) Nothing $ Just (e, f)
  in  result @?= Right (vector [0.6666666666666669,
                                1.3333333333333335],
                        -8.22222222222222)

-- sqpTest :: Assertion
-- sqpTest =
--   let a = matrix 2 [1, (-1), (-1), 2]
--       b = vector [(-2), (-6)]
--       e = matrix 2 [ (1), (1)
--                    , (-1), (2)
--                    , (2), (1)
--                    , (-1), 0
--                    , 0, (-1) ]
--       f = vector [(-2), (-2), (-3), 0, 0]
--       affineIneqs x = (e, f)
--       ineqs x = (e #> x) + f
--       cost x = (x `dot` ((0.5 * (a #> x)) + b)) + 0.0
--       problem = Problem
--                 { _cost = cost
--                 , _approxCost = const (a, b, 0)
--                 , _trueIneqs = ineqs
--                 , _approxAffineIneqs = affineIneqs
--                 , _numVariables = 2
--                 , _numIneqs = 5
--                 }
--       (xResult, fResult) = optimize problem (vector [0, 0])
--       (xTrue, fTrue) = (vector [0.6666666666666669,
--                              1.3333333333333335],
--                         -8.22222222222222)
--       xError = norm_2 $ xTrue - xResult
--       fError = abs $ fTrue - fResult
--   in  (xError < 1e-8, fError < 1e-8) @?= (True, True)

getQuadraticApproximation :: ([Double] -> Double)
                          -> ([Double] -> [Double])
                          -> ([Double] -> [[Double]])
                          -> Vector Double
                          -> (Matrix Double, Vector Double, Double)
getQuadraticApproximation fL fL' fL'' x =
  let xL = toList x
      f' = fromList $ fL' xL
      f'' = makeSymmetric $ fromLists $ fL'' xL
      c  = fL xL - (x `dot` f') + 0.5 * (x `dot` (f'' #> x))
  in  ( f'', f' - (f'' #> x), c)

getLinearApproximation :: ([Double] -> [Double])
                       -> ([Double] -> [[Double]])
                       -> Vector Double
                       -> (Matrix Double, Vector Double)
getLinearApproximation fL fL' x =
  let xL = toList x
      f' = fromLists $ fL' xL
  in  (f', vector (fL xL) - (f' #> x))

-- rosenbrockTest :: Assertion
-- rosenbrockTest =
--   -- From: http://www.mathworks.com/help/optim/ug/constrained-nonlinear-optimization-algorithms.html#f26633
--   let f [x, y] = (1 - x)^(2 :: Int) + 100*(y - x^(2 :: Int))^(2 :: Int)
--       cost x = let xL = toList x in f xL
--       convexCost = getQuadraticApproximation f
--       g [x, y] = x^(2 :: Int) + y^(2 :: Int) - 1.5
--       ineqs x = let xL = toList x in vector [g xL]
--       convexIneqs = getLinearApproximation f
--       problem = Problem
--                 { _cost = cost
--                 , _approxCost = convexCost
--                 , _trueIneqs = ineqs
--                 , _approxAffineIneqs = convexIneqs
--                 , _numVariables = 2
--                 , _numIneqs = 1
--                 }
--       (xResult, fResult) = optimize problem (vector [(-1.9), 2.0])
--       (xTrue, fTrue) = (vector [0.9072, 0.8228],
--                         0.00861632761)
--       xError = norm_2 $ xTrue - xResult
--       fError = abs $ fTrue - fResult
--   in  (xError < 1e-4, fError < 1e-4) @?= (True, True)

makeSymmetric :: Matrix Double -> Matrix Double
makeSymmetric m =
  let upperTriIxs = concat
        [ [(r, c) | c <- [(r + 1)..(cols m - 1)]] | r <- [0..(rows m - 1)] ]
      symElems (r, c) = [ ((r, c), m ! r ! c)
                        , ((c, r), m ! r ! c) ]
  in  (diag . takeDiag) m + assoc (size m) 0.0 (concatMap symElems upperTriIxs)

-- From: http://www.mathworks.com/help/optim/ug/nonlinear-equality-and-inequality-constraints.html
nonlinearEqTest :: Assertion
nonlinearEqTest =
  let f [x, y] = exp(x) * (4*x*x + 2*y*y + 4*x*y + 2*y + 1)
      g [x, y] = [negate (x * y) - 10]
      h [x, y] = [(x * x) + y - 1]

      cost = f . toList
      approxCost = getQuadraticApproximation f (grad f) (hessian f)
      ineq = vector . g . toList
      approxIneq = getLinearApproximation g (jacobian g)
      eq = vector . h . toList
      approxEq = getLinearApproximation h (jacobian h)
      problem = Problem
                { _cost = cost
                , _approxCost = approxCost
                , _trueIneqs = ineq
                , _approxAffineIneqs = approxIneq
                , _eqs = eq
                , _approxEqs = approxEq
                , _numVariables = 2
                , _numIneqs = 1
                , _numEqs = 1
                }
      (xResult, fResult) = optimize problem $ vector [(-1.0), 1.0]
      (xTrue, fTrue) = (vector [(-0.7529), 0.4332], 1.5093)
  in  (xResult, fResult) @?= (xTrue, fTrue)

main = defaultMain tests
  where tests = [ -- testCase "qp-test" qpTest
                -- , testCase "sqp-test" sqpTest
                -- , testCase "rosenbrock-test" rosenbrockTest
                  testCase "nonlineareq-test" nonlinearEqTest
                ]
