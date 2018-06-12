# Plotting tools
import matplotlib.pyplot as plt
from gpkit.small_scripts import mag
import numpy as np

def gen_plots(sol):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    rng = []
    alt = []
    for i in range(len(sol('R_{climb}'))):
           rng.append(mag(sol('R_{climb}')[i][0]))
    for i in range(len(sol('R_{cruise}'))):
           rng.append(mag(sol('R_{cruise}')[i][0]))
    for i in range(len(sol('hft')['hft_Mission/FlightState/Altitude'])):
           alt.append(mag(sol('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rng = np.cumsum(rng)
    plt.plot(rng, alt)
    plt.ylabel('Altitude [ft]', fontsize=18)
    plt.xlabel('Down Range Distance', fontsize=18)
    plt.title('Aircraft Altitude Profile')
#    plt.savefig('M08_D8_wing_profile_drag.pdf', bbox_inches="tight")
    plt.show()

def gen_D82_plots(sol):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    rng = [0]
    alt = [0]
    for i in range(len(sol('R_{climb}'))):
           rng.append(mag(sol('R_{climb}')[i][0]))
    for i in range(len(sol('R_{cruise}'))):
           rng.append(mag(sol('R_{cruise}')[i][0]))
    for i in range(len(sol('hft')['hft_Mission/FlightState/Altitude'])):
           alt.append(mag(sol('hft')['hft_Mission/FlightState/Altitude'][i][0]))
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
    rngD8 = []
    altD8 = []
    for i in range(len(solD8('R_{climb}'))):
           rngD8.append(mag(solD8('R_{climb}')[i][0]))
    for i in range(len(solD8('Rng'))):
           rngD8.append(mag(solD8('Rng')[i][0]))
    for i in range(len(solD8('hft')['hft_Mission/FlightState/Altitude'])):
           altD8.append(mag(solD8('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rngD8 = np.cumsum(rngD8)

    rng73 = []
    alt73 = []
    for i in range(len(sol737('R_{climb}'))):
           rng73.append(mag(sol737('R_{climb}')[i][0]))
    for i in range(len(sol737('Rng'))):
           rng73.append(mag(sol737('Rng')[i][0]))
    for i in range(len(sol737('hft')['hft_Mission/FlightState/Altitude'])):
           alt73.append(mag(sol737('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rng73 = np.cumsum(rng73)

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
    rngD8 = []
    altD8 = []
    for i in range(len(solD8('R_{climb}'))):
           rngD8.append(mag(solD8('R_{climb}')[i][0]))
    for i in range(len(solD8('Rng'))):
           rngD8.append(mag(solD8('Rng')[i][0]))
    for i in range(len(solD8('hft')['hft_Mission/FlightState/Altitude'])):
           altD8.append(mag(solD8('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rngD8 = np.cumsum(rngD8)

    rngno_BLI = []
    altno_BLI = []
    for i in range(len(solno_BLI('R_{climb}'))):
           rngno_BLI.append(mag(solno_BLI('R_{climb}')[i][0]))
    for i in range(len(solno_BLI('Rng'))):
           rngno_BLI.append(mag(solno_BLI('Rng')[i][0]))
    for i in range(len(solno_BLI('hft')['hft_Mission/FlightState/Altitude'])):
           altno_BLI.append(mag(solno_BLI('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rngno_BLI = np.cumsum(rngno_BLI)

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

    #generate an altitude profile plot
    tasrng = [0, 13.68, 31.34, 59.96, 115.05, 115.05, 2875.38, 2875.38, 2906.56, 2937.74, 2968.92, 3000]
    tasalt = [0, 8750, 17500, 26250, 35000, 35000, 39677.3, 39677.3, 29758., 19838.6, 9919.3, 0]

    rng = [0]
    alt = [0]

    for i in range(len(sol('R_{climb}'))):
           rng.append(mag(sol('R_{climb}')[i][0]))
    for i in range(len(sol('R_{cruise}'))):
           rng.append(mag(sol('R_{cruise}')[i][0]))
    for i in range(len(sol('hft')['hft_Mission/FlightState/Altitude'])):
           alt.append(mag(sol('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rng = np.cumsum(rng)
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

    rng = [0]
    alt = [0]

    for i in range(len(sol('R_{climb}'))):
           rng.append(mag(sol('R_{climb}')[i][0]))
    for i in range(len(sol('R_{cruise}'))):
           rng.append(mag(sol('R_{cruise}')[i][0]))
    for i in range(len(sol('hft')['hft_Mission/FlightState/Altitude'])):
           alt.append(mag(sol('hft')['hft_Mission/FlightState/Altitude'][i][0]))

    rng = np.cumsum(rng)
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
