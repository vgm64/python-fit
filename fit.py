"""
fit.py

A fitting package for python. Based on a user experience like that of root.



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

#f,p,e = fit.fit(fit.expo, x,y)
(xf,yf),p,e = fit.fit(example_function, x,y)
plot(x,y, 'b.')
plot(xf,yf)

=====================================================================

"""
from scipy.odr import odrpack as odr
from scipy.odr import models

from numpy import exp, zeros, linspace, array, diff, average

def gaus(params, x):
    N,mu,sigma = params
    return N*exp(-0.5 * ((x-mu)/sigma)**2 )

def expo(params, x):
    const, slope = params
    return exp(const + slope*x)

def double_exp(params, x):
    exp1 = params[0]*exp(params[1]*x + params[2])
    exp2 = params[3]*exp(params[4]*x + params[5])
    return exp1*exp2

def line(params, x):
    intercept, slope = params
    return slope*x + intercept

def test_func(params, x):
    N,mu,sigma = params
    return N*exp(-0.5 * ((x-mu)/sigma)**4 )

def get_default_params(xdata, ydata, func):
    """ Some predefined algorithms that help get the everyday functions
    set up for proper bounds.
    """
    if func == gaus:
        print "gaus"
        mu = xdata[ydata.argmax()]
        sigma = xdata.mean()
        N = max(ydata) - min(ydata)
        #print "Initial parameters: ", [N, mu,sigma]
        return array([N,mu,sigma])
    elif func == expo:
        print "expo"
        const = average(xdata)
        slope = average(diff(ydata)/xdata[1:])
        return array([const,slope])
    #elif func == 
    else:
        # I can try to loop from 1 to ~5 and try to f
        for i in xrange(10):
            try:
                args = [range(i), 0]
                #args = [0,range(i)]
                result = func(*args)
            except:
                pass
            else:
                break
        params = zeros(i) + xdata[len(xdata)/2]
        return params

def fit(func, x,y, default_pars=None, verbose=False,itmax=200):
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

    # http://www.scipy.org/doc/api_docs/SciPy.odr.odrpack.html
    # see models.py and use ready made models!!!!
    if default_pars:
        beta0 = array(default_pars)
    else:
        beta0 = get_default_params(x,y,func)
    model_func   = models.Model(func)
    mydata = odr.Data(x, y)
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

    # The resulting fit.
    xfit = linspace( min(x), max(x), len(x)*10)
    yfit = func(fit.beta, xfit)

    return array([xfit,yfit]),coeff,err



###########################################################
# Some more function definitions.
###########################################################


# First, polynomials from order 0 to 20.

def polN_dec(func):
    def polN(*args):
        params, x = args
        degree = int(func.func_name.replace("pol",""))
        if len(params) != degree+1:
            raise TypeError
        return sum( [param*x**i for i,param in enumerate(params)] )
    return polN

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
