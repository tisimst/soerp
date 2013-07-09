"""
Generalizes mathematical operators that work on numeric objects (from the math
module) compatible with objects with uncertainty distributions
"""
from soerp import _make_UF_compatible_object
import ad.admath as admath
#import sys

__author__ = 'Abraham Lee'

__all__ = admath.__all__

def abs(x):
    return _make_UF_compatible_object(admath.abs(x))

def acos(x):
    return _make_UF_compatible_object(admath.acos(x))

def acosh(x):
    return _make_UF_compatible_object(admath.acosh(x))

def acot(x):
    return _make_UF_compatible_object(admath.acot(x))

def acoth(x):
    return _make_UF_compatible_object(admath.acoth(x))

def acsc(x):
    return _make_UF_compatible_object(admath.acsc(x))

def acsch(x):
    return _make_UF_compatible_object(admath.acsch(x))

def asec(x):
    return _make_UF_compatible_object(admath.asec(x))

def asech(x):
    return _make_UF_compatible_object(admath.asech(x))

def asin(x):
    return _make_UF_compatible_object(admath.asin(x))

def asinh(x):
    return _make_UF_compatible_object(admath.asinh(x))

def atan(x):
    return _make_UF_compatible_object(admath.atan(x))

def atanh(x):
    return _make_UF_compatible_object(admath.atanh(x))

def ceil(x):
    return _make_UF_compatible_object(admath.ceil(x))

def cos(x):
    return _make_UF_compatible_object(admath.cos(x))

def cosh(x):
    return _make_UF_compatible_object(admath.cosh(x))

def cot(x):
    return _make_UF_compatible_object(admath.cot(x))

def coth(x):
    return _make_UF_compatible_object(admath.coth(x))

def csc(x):
    return _make_UF_compatible_object(admath.csc(x))

def csch(x):
    return _make_UF_compatible_object(admath.csch(x))

def degrees(x):
    return _make_UF_compatible_object(admath.degrees(x))

def erf(x):
    return _make_UF_compatible_object(admath.erf(x))

def erfc(x):
    return _make_UF_compatible_object(admath.erfc(x))

def exp(x):
    return _make_UF_compatible_object(admath.exp(x))

def expm1(x):
    return _make_UF_compatible_object(admath.expm1(x))

def fabs(x):
    return _make_UF_compatible_object(admath.fabs(x))

def factorial(x):
    return _make_UF_compatible_object(admath.factorial(x))

def floor(x):
    return _make_UF_compatible_object(admath.floor(x))

def gamma(x):
    return _make_UF_compatible_object(admath.gamma(x))

def lgamma(x):
    return _make_UF_compatible_object(admath.lgamma(x))

def hypot(x,y):
    return _make_UF_compatible_object(admath.hypot(x,y))

def ln(x):
    return _make_UF_compatible_object(admath.ln(x))

def log(x):
    return _make_UF_compatible_object(admath.log(x))

def log10(x):
    return _make_UF_compatible_object(admath.log10(x))

def log1p(x):
    return _make_UF_compatible_object(admath.log1p(x))

def pow(x):
    return _make_UF_compatible_object(admath.pow(x))

def radians(x):
    return _make_UF_compatible_object(admath.radians(x))

def sec(x):
    return _make_UF_compatible_object(admath.sec(x))

def sech(x):
    return _make_UF_compatible_object(admath.sech(x))

def sin(x):
    return _make_UF_compatible_object(admath.sin(x))

def sinh(x):
    return _make_UF_compatible_object(admath.sinh(x))

def sqrt(x):
    return _make_UF_compatible_object(admath.sqrt(x))

def tan(x):
    return _make_UF_compatible_object(admath.tan(x))

def tanh(x):
    return _make_UF_compatible_object(admath.tanh(x))

def trunc(x):
    return _make_UF_compatible_object(admath.trunc(x))


    
