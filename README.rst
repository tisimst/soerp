===============================
``soerp`` Package Documentation
===============================

Overview
========

``soerp`` is the Python implementation of the original Fortran code `SOERP` 
by N. D. Cox to apply a second-order analysis to `error propagation`_ (or 
uncertainty analysis). The ``soerp`` package allows you to **easily** and 
**transparently** track the effects of uncertainty through mathematical 
calculations. Advanced mathematical functions, similar to those in the standard 
math_ module can also be evaluated directly.

In order to correctly use ``soerp``, the **first eight statistical moments** 
of the underlying distribution are required. These are the *mean*, *variance*, 
and then the *standardized third through eighth moments*. These can be input 
manually in the form of an array, but they can also be **conveniently 
generated** using either the **nice constructors** or directly by using the 
distributions from the ``scipy.stats`` sub-module. See the examples below for 
usage examples of both input methods. The result of all calculations generates a 
*mean*, *variance*, and *standardized skewness and kurtosis* coefficients.


Required Packages
=================

- ad_ : For first- and second-order automatic differentiation (install this first).

Suggested Packages
==================

- NumPy_ : Numeric Python

- SciPy_ : Scientific Python (the nice distribution constructors require this)

- Matplotlib_ : Python plotting library

Basic examples
==============

Let's begin by importing all the available constructors::

    >>> from soerp import *   # uv, N, U, Exp, etc.

Now, we can see that there are several equivalent ways to specify a statistical distribution, say a Normal distribution with a mean value of 10 and a standard deviation of 1:

- Manually input the first 8 moments (mean, variance, and 3rd-8th standardized central moments)::

    >>> x = uv([10, 1, 0, 3, 0, 15, 0, 105])

- Use the ``rv`` kwarg to input a distribution from the ``scipy.stats`` module::

    >>> x = uv(rv=ss.norm(loc=10, scale=1))

- Use a built-in convenience constructor (typically the easiest if you can)::

    >>> x = N(10, 1)

A Simple Example
----------------

Now let's walk through an example of a three-part assembly stack-up::

    >>> x1 = N(24, 1)  # normally distributed
    >>> x2 = N(37, 4)  # normally distributed
    >>> x3 = Exp(2)  # exponentially distributed
    >>> Z = (x1*x2**2)/(15*(1.5 + x3))

We can now see the results of the calculations in two ways:

#. The usual ``print`` statement (or simply the object if in a terminal)::

    >>> Z  # "print" is optional at the command-line
    uv(1176.45, 99699.6822917, 0.708013052944, 6.16324345127)

#. The ``describe`` class method that explains briefly what the values are::

    >>> Z.describe()
    SOERP Uncertain Value:
     > Mean...................  1176.45
     > Variance...............  99699.6822917
     > Skewness Coefficient...  0.708013052944
     > Kurtosis Coefficient...  6.16324345127

Distribution Moments
--------------------

The eight moments of any input variable (and four of any output variable) can be accessed using the ``moments`` class method, as in::

    >>> x1.moments()
    [24.0, 1.0, 0.0, 3.0000000000000053, 0.0, 15.000000000000004, 0.0, 105.0]
    >>> Z.moments()
    [1176.45, 99699.6822917, 0.708013052944, 6.16324345127]

Correlations
------------

Statistical correlations are correctly handled, even after calculations have taken place::

    >>> x1 - x1
    0.0
    >>> square = x1**2
    >>> square - x1*x1
    0.0

Derivatives
-----------

Derivatives with respect to original variables are calculated via the ad_ package and are accessed using the **intuitive class methods**::

    >>> Z.d(x1)  # dZ/dx1
    45.63333333333333

    >>> Z.d2(x2)  # d^2Z/dx2^2
    1.6

    >>> Z.d2c(x1, x3)  # d^2Z/dx1dx3 (order doesn't matter)
    -22.816666666666666
    
When we need multiple derivatives at a time, we can use the ``gradient`` and ``hessian`` class methods::

    >>> Z.gradient([x1, x2, x3])
    [45.63333333333333, 59.199999999999996, -547.6]

    >>> Z.hessian([x1, x2, x3])
    [[0.0, 2.466666666666667, -22.816666666666666], [2.466666666666667, 1.6, -29.6], [-22.816666666666666, -29.6, 547.6]]

Error Components/Variance Contributions
---------------------------------------

Another useful feature is available through the ``error_components`` class method that has various ways of representing the first- and second-order variance components::

    >>> Z.error_components(pprint=True)
    COMPOSITE VARIABLE ERROR COMPONENTS
    uv(37.0, 16.0, 0.0, 3.0) = 58202.9155556 or 58.378236%
    uv(24.0, 1.0, 0.0, 3.0) = 2196.15170139 or 2.202767%
    uv(0.5, 0.25, 2.0, 9.0) = -35665.8249653 or 35.773258%

Advanced Example
----------------

Here's a *slightly* more advanced example, estimating the statistical properties of volumetric gas flow through an orifice meter::

    >>> from soerp.umath import *  # sin, exp, sqrt, etc.
    >>> H = N(64, 0.5)
    >>> M = N(16, 0.1)
    >>> P = N(361, 2)
    >>> t = N(165, 0.5)
    >>> C = 38.4
    >>> Q = C*umath.sqrt((520*H*P)/(M*(t + 460)))
    >>> Q.describe()
    SOERP Uncertain Value:
     > Mean...................  1330.99973939
     > Variance...............  58.210762839
     > Skewness Coefficient...  0.0109422068056
     > Kurtosis Coefficient...  3.00032693502
 
This seems to indicate that even though there are products, divisions, and the usage of ``sqrt``, the result resembles a normal distribution (i.e., Q ~ N(1331, 7.63), where the standard deviation = sqrt(58.2) = 7.63).

Main Features
=============

1. **Transparent calculations** with derivatives automatically calculated. 
   **No or little modification** to existing code required.

2. Basic `NumPy` support without modification. Vectorized calculations built-in  
   to the ``ad`` package.

3. Nearly all standard `math`_ module functions supported through the 
   ``soerp.umath`` sub-module. If you think a function is in there, it probably 
   is.

4. Nearly all derivatives calculated analytically using ``ad`` functionality.

5. **Easy continuous distribution constructors**: 

   - ``N(mu, sigma)`` : `Normal distribution`_

   - ``U(a, b)`` : `Uniform distribution`_

   - ``Exp(lamda, [mu])`` : `Exponential distribution`_

   - ``Gamma(k, theta)`` : `Gamma distribution`_

   - ``Beta(alpha, beta, [a, b])`` : `Beta distribution`_

   - ``LogN(mu, sigma)`` : `Log-normal distribution`_

   - ``X2(k)`` : `Chi-squared distribution`_

   - ``F(d1, d2)`` : `F-distribution`_

   - ``Tri(a, b, c)`` : `Triangular distribution`_

   - ``T(v)`` : `T-distribution`_

   - ``Weib(lamda, k)`` : `Weibull distribution`_

   The location, scale, and shape parameters follow the notation in the 
   respective Wikipedia articles. *Discrete distributions are not recommended 
   for use at this time. If you need discrete distributions, try the* mcerp_ 
   *python package instead.*

Installation
============

**Make sure you install the** `ad`_ **package first!**

You have several easy, convenient options to install the ``soerp`` package 
(administrative privileges may be required)

1. Download the package files below, unzip to any directory, and run::

    $ [sudo] python setup.py install
   
2. Simply copy the unzipped ``soerp-XYZ`` directory to any other location that 
   python can find it and rename it ``soerp``.
   
3. If ``setuptools`` is installed, run::

    $ easy_install --upgrade soerp
   
4. If ``pip`` is installed, run::

    $ pip install --upgrade soerp

Python 3
--------

To use this package with Python 3.x, you will need to run the ``2to3`` 
conversion tool at the command-line using the following syntax while in 
the unzipped ``soerp`` directory::

    $ 2to3 -w -f all *.py
    
This should take care of the main changes required. Then, run
``python3 setup.py install``. If bugs continue to pop up,
please email the author.
    
See Also
========

- uncertainties_ : First-order error propagation

- mcerp_ : Real-time latin-hypercube sampling-based Monte Carlo error propagation

Contact
=======

Please send **feature requests, bug reports, or feedback** to 
`Abraham Lee`_.

Acknowledgements
================

The author wishes to thank `Eric O. LEBIGOT`_ who first developed the
`uncertainties`_ python package (for first-order error propagation), 
from which many inspiring ideas (like maintaining object correlations, etc.) 
are re-used and/or have been slightly evolved. *If you don't need second
order functionality, his package is an excellent alternative since it is
optimized for first-order uncertainty analysis.*

References
==========

- N.D. Cox, 1979, *Tolerance Analysis by Computer*, Journal of Quality Technology, Vol. 11, No. 2, pp. 80-87



.. _error propagation: http://en.wikipedia.org/wiki/Propagation_of_uncertainty
.. _math: http://docs.python.org/library/math.html
.. _ad: http://pypi.python.org/pypi/ad
.. _mcerp: http://pypi.python.org/pypi/mcerp
.. _NumPy: http://www.numpy.org/
.. _SciPy: http://scipy.org
.. _Matplotlib: http://matplotlib.org/
.. _uncertainties: http://pypi.python.org/pypi/uncertainties
.. _Abraham Lee: mailto: tisimst@gmail.com
.. _Eric O. LEBIGOT: http://www.linkedin.com/pub/eric-lebigot/22/293/277
.. _PEP8: http://www.python.org/dev/peps/pep-0008
.. _Normal distribution: http://en.wikipedia.org/wiki/Normal_distribution
.. _Uniform distribution: http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
.. _Exponential distribution: http://en.wikipedia.org/wiki/Exponential_distribution
.. _Gamma distribution: http://en.wikipedia.org/wiki/Gamma_distribution
.. _Beta distribution: http://en.wikipedia.org/wiki/Beta_distribution
.. _Log-normal distribution: http://en.wikipedia.org/wiki/Log-normal_distribution
.. _Chi-squared distribution: http://en.wikipedia.org/wiki/Chi-squared_distribution
.. _F-distribution: http://en.wikipedia.org/wiki/F-distribution
.. _Triangular distribution: http://en.wikipedia.org/wiki/Triangular_distribution
.. _T-distribution: http://en.wikipedia.org/wiki/Student's_t-distribution
.. _Weibull distribution: http://en.wikipedia.org/wiki/Weibull_distribution
