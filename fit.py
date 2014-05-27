"""
fit.py

A fitting package for python. Designed for ease of use for fast results.

List of built in functions that will try to return intelligent default
parameters for you when you use them:
  gaus
  expo
  double_exp
  line
  crystal_ball
  crystal_ball_norm

Built in functions that do not provide intelligent default values:
  fofoi - First over first inverse: a/(x+b) + c*x+d


================== Example Usage ===================================

from pylab import *
import fit
ion()


x = (4.2105303, 5.2631601, 6.2405997, 7.5187997, 8.7218, 9.7744402, 10.676691, 11.65414, 12.63158, 13.83459, 14.887219, 16.015039, 17.06767, 18.270679, 19.24812, 20.300751, 21.50376, 23.157888, 25.789471, 28.345871, 30.601501, 33.458643, 39.022559, 46.015039, 48.270679)

y = (0.18942001, 0.2099, 0.23891001, 0.27816002, 0.31911, 0.35836001, 0.39932001, 0.43686003, 0.46416002, 0.49829001, 0.51536004, 0.52556, 0.51876995, 0.5, 0.47271, 0.44026, 0.39249001, 0.33106002, 0.24060, 0.17746, 0.13311001, 0.11262, 0.095566, 0.095566, 0.095566)


x = array(x)
y = array(y)

# Guassian with linear background.
def example_function(params, x):
    N,mu,sigma,a,b,c = params
    return N*exp(-0.5 * ((x-mu)/sigma)**2 ) + a*x**2 + b*x + c

#f,p,e,chi = fit.fit(fit.expo, x,y)
(xf,yf),p,e,chi = fit.fit(example_function, x,y)
plot(x,y, 'b.')
plot(xf,yf)

=====================================================================

"""
from scipy.odr import odrpack as odr
from scipy.odr import models

from numpy import exp, zeros, linspace, array, diff, average, logical_and, argmax

def gaus(params, x):
  """ An un-normalized Gaussian curve.
  params: N, mu, sigma """
  N,mu,sigma = params
  return N*exp(-0.5 * ((x-mu)/sigma)**2 )

def expo(params, x):
  """ An exponential curve. *NB* The constant is in the exponent! 
  params: const, slope """
  const, slope = params
  return exp(const + slope*x)

def double_exp(params, x):
  """ Two exponential constants in one.
  params: const1, slope1, mu1, const2, slope2, mu2 """
  exp1 = params[0]*exp(params[1]*x + params[2])
  exp2 = params[3]*exp(params[4]*x + params[5])
  return exp1*exp2

def line(params, x):
  """ Just another name for pol1. """
  intercept, slope = params
  return slope*x + intercept

def crystal_ball(params, x):
  """ A Gaussian curve on one side and a power-law on the other side. Used in
  physics to model lossy processes.
  See http://en.wikipedia.org/wiki/Crystal_Ball_function
  Note that the definition used here differs slightly. At the time of this
  writing, the wiki article has some discrepancies in definitions/plots. This
  definition makes it easier to fit the function by using complex numbers
  and by negating any negative values for a and n.

  This version of the crystal ball is normalized by an additional parameter.
  params: N, a, n, xb, sig
  """
  x = x+0j # Prevent warnings...
  N, a, n, xb, sig = params
  if a < 0:
    a = -a
  if n < 0:
    n = -n
  aa = abs(a)
  A = (n/aa)**n * exp(- aa**2 / 2)
  B = n/aa - aa
  total = 0.*x
  total += ((x-xb)/sig  > -a) * N * exp(- (x-xb)**2/(2.*sig**2))
  total += ((x-xb)/sig <= -a) * N * A * (B - (x-xb)/sig)**(-n)
  try:
    return total.real
  except:
    return total
  return total

def crystal_ball_norm(params, x):
  """ A Gaussian curve on one side and a power-law on the other side. Used in
  physics to model lossy processes.
  See http://en.wikipedia.org/wiki/Crystal_Ball_function
  Note that the definition used here differs slightly. At the time of this
  writing, the wiki article has some discrepancies in definitions/plots. This
  definition makes it easier to fit the function by using complex numbers
  and by negating any negative values for a and n.

  This version of the crystal ball is normalized by an internal normalization
  process.
  params: a, n, xb, sig
  """
  x = x+0j # Prevent warnings...
  a, n, xb, sig = params
  if a < 0:
    a = -a
  if n < 0:
    n = -n
  aa = abs(a)
  A = (n/aa)**n * exp(- aa**2 / 2)
  B = n/aa - aa
  C = n/aa / (n-1.) * exp(-aa**2/2.)
  D = sqrt(pi/2.) * (1. + erf(aa/sqrt(2.)))
  N = 1. / (sig * (C+D))
  total = 0.*x
  total += ((x-xb)/sig  > -a) * N * exp(- (x-xb)**2/(2.*sig**2))
  total += ((x-xb)/sig <= -a) * N * A * (B - (x-xb)/sig)**(-n)
  try:
    return total.real
  except:
    return total
  return total

def pow_law(params, x):
  """ Power law curve.
  params: base,exponent
  """
  base,exponent = params
  return base*(x**exponent)

def fofoi(params, x):
  """ fofoi - First over first inverse: a/(x+b) + c*x+d
  This function is useful when exponential curves and power-law curves don't
  quite capture a dual slope function.
  params: a, b, c, d
  """
  a,b,c,d = params
  return a/(x+b) + c*x+d

def get_default_params(xdata, ydata, func):
  """ Some predefined algorithms that help get the everyday functions
  set up for proper bounds.
  """
  if func == gaus:
    mu = xdata[ydata.argmax()]
    sigma = xdata.std()
    N = max(ydata)
    return array([N,mu,sigma])
  elif func == expo:
    const = average(xdata)
    slope = average(diff(ydata)/xdata[1:])
    return array([const,slope])
  elif func == crystal_ball:
    N = ydata.max()
    xb = xdata[ydata.argmax()]
    sigma = sum(ydata > ydata.mean()*1.8)*1. / len(xdata) * (xdata.max() - xdata.min())
    n = 2.
    a = .5
    return array([N, a, n, xb, sigma])
  elif func == crystal_ball_norm:
    xb = xdata[ydata.argmax()]
    sigma = sum(ydata > ydata.mean()*1.8)*1. / len(xdata) * (xdata.max() - xdata.min())
    n = 2.
    a = .5
    return array([a, n, xb, sigma])
  elif func == pow_law:
    base = 2.6
    exponent = -.15
    return array([base, exponent])
  else:
    # Just put some values from the x's into the params.
    for i in xrange(21):
      try:
        args = [range(i), 0]
        result = func(*args)
      except ValueError:
        pass
      else:
        break
    params = zeros(i) + xdata[len(xdata)/2]
    return params

def fit(func, x,y, default_pars=None, data_range=None, we=None, verbose=False, itmax=200):
  ''' The meat of the fitting package. See docs of fit.py for more details.
  Functions available are gaus and expo and more.

  Error implementation provided via an example by: Tiago, 20071114

  Performs a least squares fit to the data, with errors!
  Uses scipy odrpack, but for least squares.
  
  IN: (func, x, y, verbose, itmax)
     func         - A function that accepts input in the form:
                    func(params, x)
     x,y (arrays) - data to fit
     default_pars - Optional default parameters to start minimization at.
     data_range   - Fit a subrange of (x,y). Provide a tuple of the form
                    (x_min, x_max).
     we           - Weighting for data points as delivered to ODR. You
                    probably want it the same length as x.
     verbose      - can be 0,1,2 for different levels of output
                    (False or True are the same as 0 or 1)
     itmax (int)  - optional maximum number of iterations
     
  OUT: (fit, params, err)
     fit    -  xfit,yfit arrays that can immediately plotted to see the
              results of your fit. xf and yf are defined as: 
                  xfit = linspace( min(x), max(x), len(x)*10)
                  yfit = func(xfit)
     params - the coefficients of your fit in the order the function takes.
     err    - standard error (1-sigma) on the coefficients

  '''

  # If this is a histogram output, correct it for the user.
  if len(x) == len(y) + 1:
    x = (x[1:] + x[:-1])/2.
  # Take a slice of data.
  if data_range:
    y = y[logical_and(x > data_range[0], x < data_range[1])]
    x = x[logical_and(x > data_range[0], x < data_range[1])]

  # http://www.scipy.org/doc/api_docs/SciPy.odr.odrpack.html
  # see models.py and use ready made models!!!!
  if default_pars != None:
    beta0 = array(default_pars)
  else:
    beta0 = get_default_params(x,y,func)
  model_func   = models.Model(func)
  
  mydata = odr.Data(x, y, we)
  myodr  = odr.ODR(mydata, model_func,maxit=itmax, beta0=beta0)

  # Set type of fit to least-squares:
  myodr.set_job(fit_type=2)
  if verbose == 2: myodr.set_iprint(final=2)
        
  fit = myodr.run()

  # Display results:
  if verbose: fit.pprint()

  if fit.stopreason[0] == 'Iteration limit reached':
    print '(WWW) poly_lsq: Iteration limit reached, result not reliable!'

  # Results and errors
  coeff = fit.beta
  err   = fit.sd_beta
  chi   = fit.sum_square

  # The resulting fit.
  xfit = linspace( min(x), max(x), len(x)*10)
  yfit = func(fit.beta, xfit)

  return array([xfit,yfit]),coeff,err, chi



###########################################################
# Some more function definitions.
###########################################################


# First, create a polynomial decorator.

def polN_dec(func):
  """ Polynomial decorator """
  degree = int(func.func_name.replace("pol",""))
  def polN(*args):
    params, x = args
    if len(params) != degree+1:
      raise ValueError
    return sum( [param*x**i for i,param in enumerate(params)] )
  polN.__doc__ = \
  """A polynomial function of degree {}. Parameters are provided in
  increasing order of degree, i.e. [0] + [1]*x + [2]*x^2 + ...
  Input: (params, x)
    params - An array of length {}
    x      - The data to evaluate given the previous parameters.
  """
  polN.__doc__ = polN.__doc__.format(degree, degree+1)
  polN.__name__ = "pol{}".format(degree)
  return polN

# Create polynomials from 0 to 20.
@polN_dec
def pol0(params,x): return
@polN_dec
def pol1(params,x): return
@polN_dec
def pol2(params,x): return
@polN_dec
def pol3(params,x): return
@polN_dec
def pol4(params,x): return
@polN_dec
def pol5(params,x): return
@polN_dec
def pol6(params,x): return
@polN_dec
def pol7(params,x): return
@polN_dec
def pol8(params,x): return
@polN_dec
def pol9(params,x): return
@polN_dec
def pol10(params,x): return
@polN_dec
def pol11(params,x): return
@polN_dec
def pol12(params,x): return
@polN_dec
def pol13(params,x): return
@polN_dec
def pol14(params,x): return
@polN_dec
def pol15(params,x): return
@polN_dec
def pol16(params,x): return
@polN_dec
def pol17(params,x): return
@polN_dec
def pol18(params,x): return
@polN_dec
def pol19(params,x): return
@polN_dec
def pol20(params,x): return
