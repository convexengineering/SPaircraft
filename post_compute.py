"""
function to compute aircraft values after solving model
"""

from gpkit.small_scripts import mag

def post_compute(sol, Nclimb):
    eta_P_fan = 2/(1 + 1.94384*mag(sol('u_{8}')[Nclimb])/mag(sol('V')[Nclimb]))

    print "Fan Propulsive Efficiency in Cruise Segment 1"
    print "---------------------"
    print eta_P_fan
