python-fit
==========

The python-fit module is designed for people who need to fit data
frequently and quickly. The module is not designed for huge amounts of
control over the minimization process but rather tries to make fitting
data simple and painless. If you want to fit data several times a day,
every day, and you really just want to see if the fit you've made looks
good against your data, check out this software. If you want one very
statistically aware and neurotically controlled fit, you might consider
looking elsewhere.

With python-fit, you get work done.

Features
--------

-  Common fitting curves built-in
-  Default parameters for built-in functions intelligently calculated
   using your data.
-  Fit with user defined functions, too.
-  A ready-to-plot-fit always conviently returned.
-  Get fit parameters and associated errors.
-  Chi-squared residual.

Example Usage
-------------

Fit with built in functions:

.. code-block:: python

    from pylab import *
    ion()
    import fit
    from numpy import random, exp
    random.seed(0)

    # Create some data to fit
    x = arange(-10, 10, .2)
    # A gaussian of height 10, width 2, centered at zero. With noise.
    y = 10*exp(-x**2/8) + (random.rand(100) - 0.5)


    # No need to provide first guess at parameters for fit.gaus
    (xf, yf), params, err, chi = fit.fit(fit.gaus, x,y)

    print "N:    %.2f +/- %.3f" % (params[0], err[0])
    print "N:    %.2f +/- %.3f" % (params[1], err[1])
    print "N:    %.2f +/- %.3f" % (params[2], err[2])

    plot(x,y, 'bo', label='Data')
    plot(xf, yf, 'r-', label='Fit')
    legend() 

Functions available:

.. code-block:: python

    Gaussian curve
    Exponential curve
    Double exponential
    Polynomials for degrees 0-20
    Power-Law
    Crystal Ball
    ... and more!

Fit with user defined functions:

.. code-block:: python

    x = (4.2105303, 5.2631601, 6.2405997, 7.5187997, 8.7218, 9.7744402, 10.676691, 11.65414, 12.63158, 13.83459, 14.887219, 16.015039, 17.06767, 
    18.270679, 19.24812, 20.300751, 21.50376, 23.157888, 25.789471, 28.345871, 30.601501, 33.458643, 39.022559, 46.015039, 48.270679)
    y = (0.18942001, 0.2099, 0.23891001, 0.27816002, 0.31911, 0.35836001, 0.39932001, 0.43686003, 0.46416002, 0.49829001, 0.51536004, 0.52556, 0.51876995, 
    0.5, 0.47271, 0.44026, 0.39249001, 0.33106002, 0.24060, 0.17746, 0.13311001, 0.11262, 0.095566, 0.095566, 0.095566)

    x = array(x)
    y = array(y)

    # Guassian with quadratic background.
    def example_function(params, x):
        N,mu,sigma,a,b,c = params
        return N*exp(-0.5 * ((x-mu)/sigma)**2 ) + a*x**2 + b*x + c
      
    # It will still try to guess parameters, but they are dumb!
    (xf,yf),p,e,chi = fit.fit(example_function, x,y)
    plot(x,y, 'bo', label='Data')
    plot(xf,yf, 'r-', label='Fit')
    legend()

Even though ``example_function`` is defined by the user, python-fit will
guess parameters (the median value of the xdata for all parameters; it
works if x and y are on similar scales). If the fit fails, then provide
some decent parameters as a first guess:

.. code-block:: python

    results = fit.fit(example_function, x, y, default_pars = [1, 12, 10, 1, 1, 1])
    plot(results[0][0], results[0][1], 'r--')

Fit a sub-range:

.. code-block:: python

    clf()
    results = fit.fit(fit.gaus, x, y, data_range=[0, 23])
    plot(results[0][0], results[0][1], 'r-.')

Define your own weights to prevent outliers from wreaking havoc on your
fit:

.. code-block:: python

    # Create some outliers.
    y_outlier = y + (random.rand(len(y))**20)*3
    # I'll just make a cut and say outliers are above 0.55
    weights = 1. * (y_outlier < .55)
    results = fit.fit(example_function, x, y_outlier, we=weights)
    clf()
    plot(x,y_outlier, 'bo', label='Data w/ Outliers')
    plot(results[0][0], results[0][1], 'r-.', label='Fit around outliers')
