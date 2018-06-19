# Plotting tools
import matplotlib.pyplot as plt
from gpkit.small_scripts import mag
import numpy as np

def gen_plots(sol):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    rng = np.cumsum(mag(sol('R_{segment}')))
    alt = mag(sol('hft'))
    plt.plot(rng, alt)
    plt.ylabel('Altitude [ft]', fontsize=18)
    plt.xlabel('Down Range Distance', fontsize=18)
    plt.title('Aircraft Altitude Profile')
    plt.show()

def gen_D82_plots(sol):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    rng = mag(sol('R_{segment}'))
    alt = mag(sol('hft'))
    rng = np.cumsum(rng)

    tasrng = [0, 11.51, 27.36, 52.64, 103.28, 103.28, 2825.41, 2825.41, 2869.08, 2912.76, 2956.43, 3000]
    tasalt = [0, 9619.5, 19239.0, 28858.5, 38478.0, 38478.0, 41681.3, 41681.3, 32129.3, 21998.5, 11288.7, 0]

    plt.plot(rng, alt)
    plt.plot(tasrng, tasalt)
    plt.legend(['SP Model', 'TASOPT'], loc=4, fontsize=18)
    plt.ylabel('Altitude [ft]', fontsize=22)
    plt.xlabel('Down Range Distance [nm]', fontsize=22)
    plt.title('D8.2 Altitude Profile',fontsize=22)
    plt.ylim([0, 46000])
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tick_params(axis='both', which='minor', labelsize=16)
    plt.savefig('D8_altitude_profile.eps', bbox_inches="tight")
    plt.show()

def gen_D8_737_plots(solD8, sol737):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    rngD8 = np.cumsum(mag(solD8('R_{segment}')))
    altD8 = mag(solD8('hft'))
    rng73 = np.cumsum(mag(sol737('R_{segment}')))
    alt73 = mag(sol737('hft'))

    tasrngD8 = [0, 11.51, 27.36, 52.64, 103.28, 103.28, 2825.41, 2825.41, 2869.08, 2912.76, 2956.43, 3000]
    tasaltD8 = [0, 9619.5, 19239.0, 28858.5, 38478.0, 38478.0, 41681.3, 41681.3, 32129.3, 21998.5, 11288.7, 0]

    tasrng73 = [0, 13.68, 31.34, 59.96, 115.05, 115.05, 2875.38, 2875.38, 2906.56, 2937.74, 2968.92, 3000]
    tasalt73 = [0, 8750, 17500, 26250, 35000, 35000, 39677.3, 39677.3, 29758., 19838.6, 9919.3, 0]

    plt.plot(rngD8, altD8)
    plt.plot(tasrngD8, tasaltD8)
    plt.plot(rng73, alt73)
    plt.plot(tasrng73, tasalt73)
    plt.legend(['D8 SP Model', 'D8 TASOPT', '737 SP Model', '737 TASOPT'], loc=4)
    plt.ylabel('Altitude [ft]', fontsize=18)
    plt.xlabel('Down Range Distance [nm]', fontsize=18)
    plt.title('737 and D8 Altitude Profile', fontsize=18)
    plt.savefig('737_D8_altitude_profile.eps', bbox_inches="tight")
    plt.show()

def gen_D8_D8_no_BLI_plots(solD8, solno_BLI):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    rngD8 = np.cumsum(mag(solD8('R_{segment}')))
    altD8 = mag(solD8('hft'))
    rngno_BLI = np.cumsum(mag(solno_BLI('R_{segment}')))
    altno_BLI = mag(solno_BLI('hft'))

    plt.plot(rngD8, altD8)
    plt.plot(rngno_BLI, altno_BLI)
    plt.legend(['D8', 'D8 w/out BLI (rear podded engines)'], loc=4)
    plt.ylabel('Altitude [ft]', fontsize=18)
    plt.xlabel('Down Range Distance [nm]', fontsize=18)
    plt.title('D8 Altitude Profile with and without BLI')
    plt.savefig('D8_D8_no_BLI_altitude_profile.pdf', bbox_inches="tight")
    plt.show()

def gen_737_plots(sol):
    """
    function to generate plots of interesting values
    """
    rng = np.cumsum(mag(sol('R_{segment}')))
    alt = mag(sol('hft'))

    #generate an altitude profile plot
    tasrng = [0, 13.68, 31.34, 59.96, 115.05, 115.05, 2875.38, 2875.38, 2906.56, 2937.74, 2968.92, 3000]
    tasalt = [0, 8750, 17500, 26250, 35000, 35000, 39677.3, 39677.3, 29758., 19838.6, 9919.3, 0]

    plt.plot(rng, alt)
    plt.plot(tasrng, tasalt)
    plt.legend(['SP Model', 'TASOPT'], loc=4, fontsize=18)
    plt.ylabel('Altitude [ft]', fontsize=22)
    plt.xlabel('Down Range Distance [nm]', fontsize=22)
    plt.title('737 Altitude Profile', fontsize=22)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tick_params(axis='both', which='minor', labelsize=16)
    plt.savefig('737_altitude_profile.eps', bbox_inches="tight")
    plt.show()

def gen_777_plots(sol):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    tasrng = [0, 15.6, 33.24, 60.40, 107.98, 107.98, 5850.37, 5850.37, 5887.82, 5925.28, 5962.74, 6000]
    tasalt = [0, 7994.2, 15988.5, 23982.8, 31977.0, 31977.0, 39723.4, 39723.4, 31282.2, 21847.9, 11420.5, 0]

    rng = np.cumsum(mag(sol('R_{segment}')))
    alt = mag(sol('hft'))

    plt.plot(rng, alt)
    plt.plot(tasrng, tasalt)
    plt.legend(['SP Model', 'TASOPT'], loc=4, fontsize=18)
    plt.ylabel('Altitude [ft]', fontsize=22)
    plt.xlabel('Down Range Distance [nm]', fontsize=22)
    plt.title('777 Altitude Profile', fontsize=22)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tick_params(axis='both', which='minor', labelsize=16)
    plt.savefig('777_altitude_profile.eps', bbox_inches="tight")
    plt.show()
