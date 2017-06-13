"""
function to compute aircraft values after solving model
"""

from gpkit.small_scripts import mag

def post_compute(sol, Nclimb):
    eta_P_fan = 2/(1 + 1.94384*mag(sol('u_8')[Nclimb][0])/mag(sol('V')['V_Mission/FlightState'][Nclimb][0]))

    print "Fan Propulsive Efficiency in Cruise Segment 1"
    print "---------------------"
    print eta_P_fan
