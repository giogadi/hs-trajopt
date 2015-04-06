import Test.HUnit
import Test.Framework
import Test.Framework.Providers.HUnit

import Numeric.LinearAlgebra.HMatrix
import Numeric.Minimization.QuadProgPP

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

main = defaultMain tests
  where tests = [ testCase "qp-test" qpTest ]
