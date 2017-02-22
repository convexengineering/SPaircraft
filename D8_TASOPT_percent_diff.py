"""
Method to compute and print percent differences between SP D8 model and TASOPT
"""

from gpkit.small_scripts import mag

def percent_diff(sol, version):
    """
    Method to compute and print percent differences between SP D8 model and TASOPT

    INPUTS
    ------
    sol: solution from D8 SP model
    version: either 0, 1, or 2....corresponds to comparing to the D8.0, D8.1, or D8.2
    """
    if version == 2:
        #weights to compare
        print "WEIGHT DIFFERENCES"
        print "\n"
        print "Total Fuel Weight Percent Diff: %s" % compute_diff(mag(sol('W_{f_{total}}').to('lbf')), 26959.1)
        print "\n"
        print "Total Aircraft Weight Percent Diff: %s" % compute_diff(mag(sol('W_{total}').to('lbf')), 134758.2)
        print "\n"
        print "Engine Weight Percent Diff: %s" % compute_diff(mag(sol('W_{engine}').to('lbf')), 8168.2)
        print "\n"
        print "Fuselage Weight Percent Diff: %s" % compute_diff(mag(sol('W_{fuse}').to('lbf')), 30487.4)
        print "\n"
        print "Payload Weight Percent Diff: %s" % compute_diff(mag(sol('W_{payload}').to('lbf')), 38715.5)
        print "\n"
        print "VT Weight Percent Diff: %s" % compute_diff(mag(sol('W_{VT}').to('lbf')), 376.4)
        print "\n"
        print "HT Weight Percent Diff: %s" % compute_diff(mag(sol('W_{HT}').to('lbf')), 1578.5)
        print "\n"
        print "Wing Weight Percent Diff: %s" % compute_diff(mag(sol('W_{wing}').to('lbf')), 19737.5)

        #wing value to compare
        print "\n\n\n"
        print "WING DIFFERENCES"
        print "\n"
        print "Wing Span Percent Diff: %s" % compute_diff(mag(sol('b').to('ft')), 140.042)
        print "\n"
        print "Wing Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR')), 15.7490)
        print "\n"
        print "Wing Area Percent Diff: %s" % compute_diff(mag(sol('S').to('ft^2')), 1245.27)
        

        #HT valuees to compare
        print "\n\n\n"
        print "HORIZONTAL TAIL DIFFERENCES"
        print "\n"
##        print "HT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR_h')), 12.0000)
        print "\n"
        print "HT Area Percent Diff: %s" % compute_diff(mag(sol('S_h').to('ft^2')), 243.76  )
        print "\n"
        print "HT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{h}')), .913) 
        

        #VT values to compare
        print "\n\n\n"
        print "VERTICAL TAIL DIFFERENCES"
        print "\n"
##        print "VT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('A_{vt}')), 2.2000)
        print "\n"
        print "VT Span Percent Diff: %s" % compute_diff(mag(sol('b_{vt}').to('ft')), 11.59)
        print "\n"
        print "VT Area Percent Diff: %s" % compute_diff(mag(sol('S_{vt}').to('ft^2')), 122.03)
        print "\n"
        print "VT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{vt}')), .03) 

def compute_diff(sp, tasopt):
    """
    Method to actually compute the percent difference
    """
    diff = 100*(sp-tasopt)/sp
    
    return diff
