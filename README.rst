``soerp`` Package Documentation
===============================

.. contents::

Overview
--------

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
-----------------

- ad_ : For automatic differentiation (install this first).

Suggested Packages
------------------

- NumPy_ : Numeric Python

- SciPy_ : Scientific Python (the nice distribution constructors require this)

- Matplotlib_ : Python plotting library

Basic examples
--------------
::

    >>> from soerp import *   # the constructors for uncertain variables

    # these are equivalent ways to specify the distribution
    >>> x = uv([10, 1, 0, 3, 0, 15, 0, 105])  # manually input moments
    >>> x = uv(rv=ss.norm(loc=10, scale=1))  # directly using the scipy.stats distributions
    >>> x = N(10, 1)  # a nice constructor (still requires scipy)

    >>> x1 = N(24, 1)  # normally distributed
    >>> x2 = N(37, 4)  # normally distributed
    >>> x3 = Exp(2)  # exponentially distributed

    >>> Z = (x1*x2**2)/(15*(1.5 + x3))
    >>> Z  # output compactly shows the mean, variance, and standardized skewness and kurtosis
    uv(1176.45, 99699.6822917, 0.708013052944, 6.16324345127)

    >>> Z.describe()  # use for more detailed output
    SOERP Uncertain Value:
     > Mean...................  1176.45
     > Variance...............  99699.6822917
     > Skewness Coefficient...  0.708013052944
     > Kurtosis Coefficient...  6.16324345127
        
    >>> x1.moments()  # the eight moments can be accessed at any time
    [24.0, 1.0, 0.0, 3.0000000000000053, 0.0, 15.000000000000004, 0.0, 105.0]
    
    >>> x1-x1  # correlations are correctly handled
    0.0
    
    # convenient access to derivatives
    >>> Z.d(x1)  # first derivative wrt x1 (returns all if no input, 0 if derivative doesn't exist)
    45.63333333333333

    >>> Z.d2(x2)  # second derivative wrt x2
    1.6

    >>> Z.d2c(x1, x3)  # second cross-derivative wrt x1 and x3 (order doesn't matter)
    -22.816666666666666
    
    >>> Z.gradient([x1, x2, x3])  # convenience function, useful for optimization
    [45.63333333333333, 59.199999999999996, -547.6]

    >>> Z.hessian([x1, x2, x3])   # another convenience function
    [[0.0, 2.466666666666667, -22.816666666666666], [2.466666666666667, 1.6, -29.6], [-22.816666666666666, -29.6, 547.6]]

    >>> Z.error_components(pprint=True)  # show how each variable is contributing errors
    COMPOSITE VARIABLE ERROR COMPONENTS
    uv(37.0, 16.0, 0.0, 3.0) = 58202.9155556 or 58.378236%
    uv(24.0, 1.0, 0.0, 3.0) = 2196.15170139 or 2.202767%
    uv(0.5, 0.25, 2.0, 9.0) = -35665.8249653 or 35.773258%

    # a more advanced example (volumetric gas flow through orifice meter)
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
 
Main Features
-------------

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
   respective Wikipedia articles.

Installation
------------

**Make sure you install the** `ad`_ **package first!**

You have several easy, convenient options to install the ``soerp`` package 
(administrative privileges may be required)

1. Download the package files below, unzip to any directory, and run 
   ``python setup.py install`` from the command-line.
   
2. Simply copy the unzipped ``soerp-XYZ`` directory to any other location that 
   python can find it and rename it ``soerp``.
   
3. If ``setuptools`` is installed, run ``easy_install --upgrade soerp`` from 
   the command-line.
   
4. If ``pip`` is installed, run ``pip --upgrade soerp`` from the command-line

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
--------

- uncertainties_ : First order error propagation

- mcerp_ : Real-time Monte Carlo, Latin-Hypercube Sampling-based, Error Propagation

Contact
-------

Please send **feature requests, bug reports, or feedback** to 
`Abraham Lee`_.

Acknowledgements
----------------

A lot of the credit goes to `Eric O. LEBIGOT`_ who first developed 
`uncertainties`_, a very nice first-order package for error propagation, from 
which many inspiring ideas (like correlating variables, etc.) are re-used or 
slightly evolved. If you **don't** need second order functionality, I recommend 
using his package.

References
----------

- N.D. Cox, 1979, Tolerance Analysis by Computer, Journal of Quality Technology, Vol. 11, No. 2, pp. 80-87



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
