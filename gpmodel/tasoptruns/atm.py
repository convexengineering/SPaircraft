"air properties"
import numpy as np
import matplotlib.pyplot as plt
from gpfit.fit import fit

def atm(altitude):
    Tsl = 288.15 # K
    Latm = 0.0065 #K/m
    Tatm = Tsl - Latm*altitude
    psl = 101325.0 # Pa
    patm = (Tatm/Tsl)**5.257*psl
    Ratm = 8.31447 #J/mol/K
    Matm = 0.0289644 # kg/mol
    rho = patm/(Ratm/Matm*Tatm)
    return rho, Tatm, patm

def plotfit():
    h = np.linspace(0, 11000, 100)
    r, T, p = atm(h)
    x = np.log(r/1.225)
    y = np.log(T/288.15)
    cn, rm = fit(x, y, 1, "MA")
    yfit = cn.evaluate(x)

    fig, ax = plt.subplots()
    ax.plot(1.225*np.exp(x), 288.15*np.exp(y), "o", mfc="None")
    ax.plot(1.225*np.exp(x), 288.15*np.exp(yfit))
    ax.grid()
    ax.set_xlabel("air density [kg/m^3]")
    ax.set_ylabel("air temperature [K]")
    fig.savefig("rhofit.pdf")
