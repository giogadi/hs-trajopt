import Test.HUnit
import Test.Framework
import Test.Framework.Providers.HUnit

import Numeric.LinearAlgebra.HMatrix
import Numeric.Minimization.QuadProgPP

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
      problem = Problem
                { _cost = (a, b)
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

main = defaultMain tests
  where tests = [ testCase "qp-test" qpTest
                , testCase "sqp-test" sqpTest
                ]
