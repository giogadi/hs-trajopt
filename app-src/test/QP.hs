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

sqpTest :: Assertion
sqpTest =
  let a = matrix 2 [1, (-1), (-1), 2]
      b = vector [(-2), (-6)]
      e = matrix 2 [ (1), (1)
                   , (-1), (2)
                   , (2), (1)
                   , (-1), 0
                   , 0, (-1) ]
      f = vector [(-2), (-2), (-3), 0, 0]
      affineIneqs x = (e, f)
      ineqs x = (e #> x) + f
      cost x = (x `dot` ((0.5 * (a #> x)) + b)) + 0.0
      problem = Problem
                { _cost = cost
                , _approxCost = const (a, b, 0)
                , _trueIneqs = ineqs
                , _approxAffineIneqs = affineIneqs
                , _numVariables = 2
                , _numIneqs = 5
                }
      (xResult, fResult) = optimize problem (vector [0, 0])
      (xTrue, fTrue) = (vector [0.6666666666666669,
                             1.3333333333333335],
                        -8.22222222222222)
      xError = norm_2 $ xTrue - xResult
      fError = abs $ fTrue - fResult
  in  (xError < 1e-8, fError < 1e-8) @=? (True, True)

rosenbrockTest :: Assertion
rosenbrockTest =
  let f [x, y] = (1 - x)^(2 :: Int) + 100*(y - x^(2 :: Int))^(2 :: Int)
      cost x = let xL = toList x in f xL
      convexCost x = let xL = toList x
                         f' = fromList $ grad f xL
                         f'' = fromLists $ hessian f xL
                         c  = f xL - (x `dot` f') + (x `dot` (f'' #> x))
                     in  ( f'', f', c)
      g [x, y] = x^(2 :: Int) + y^(2 :: Int) - 1.5
      ineqs x = let xL = toList x in vector [g xL]
      convexIneqs x = let xL = toList x
                          g' = fromList $ grad g xL
                      in  (fromRows [g'], vector [g xL - (x `dot` g')])
      problem = Problem
                { _cost = cost
                , _approxCost = convexCost
                , _trueIneqs = ineqs
                , _approxAffineIneqs = convexIneqs
                , _numVariables = 2
                , _numIneqs = 1
                }
      (xResult, fResult) = optimize problem (vector [(-1.9), 2.0])
      (xTrue, fTrue) = (vector [0.9072, 0.8228],
                        267.62)
  in  (xResult, fResult) @=? (xTrue, fTrue)

main = defaultMain tests
  where tests = [ testCase "qp-test" qpTest
                , testCase "sqp-test" sqpTest
                , testCase "rosenbrock-test" rosenbrockTest
                ]
