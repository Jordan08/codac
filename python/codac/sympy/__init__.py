from codac._core._sympy import *
from codac import AnalyticFunction, AnalyticFunction_Scalar, AnalyticFunction_Vector

def sympy_diff(f):
  if isinstance(f.f, (AnalyticFunction_Scalar)):
    return AnalyticFunction(sympy_diff_scalar(f.f))
  elif isinstance(f.f, (AnalyticFunction_Vector)):
    return AnalyticFunction(sympy_diff_vector(f.f))
  else:
    codac_error("sympy_diff: invalid function input")