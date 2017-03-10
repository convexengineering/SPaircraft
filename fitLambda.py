import pandas as pd 
import numpy as np 
from numpy import logspace, log, log10
import matplotlib.pyplot as plt
# you have to change it to this and run the file form gpfit
from gpfit.fit import fit

# Fitting for y = (1+lam+lam^2)/(1+lam)
lambdas = np.linspace(0.1,0.35,20)
p = 1.+2.*lambdas

y = (1.+lambdas+lambdas**2.)/(1.+lambdas)
lambda_max = np.amax(lambdas)
y_max = np.amax(y)
pmax = np.amax(p)
logy = log(y/y_max)
logp = log(p/pmax)
Type = 'SMA'
K = 2
cstrt,rms_error = fit(logp[1:-2],logy[1:-2],K,Type)


geny = (0.205*(p/pmax)**0.772 + 0.795*(p/pmax)**-0.125)**(1/.166)*y_max
genyorig = (0.86*(p)**(-2.38) + 0.14*(p)**0.56)**(1/3.94)


plt.plot(lambdas,y,'r')
plt.plot(lambdas, geny, 'b')
plt.plot(lambdas,genyorig,'c')
plt.show()



