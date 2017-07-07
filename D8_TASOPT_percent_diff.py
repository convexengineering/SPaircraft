"""
Method to compute and print percent differences between SP D8 model and TASOPT
"""

from gpkit.small_scripts import mag

def percent_diff(sol,aircraft):
    """
    Method to compute and print percent differences between SP D8 model and TASOPT

    INPUTS
    ------
    sol: solution from D8 SP model
    aircraft: string corresponding to different aircraft
    b737800 compares to TASOPT 737-800
    optimal737 compares to TASOPT 737-800 w/physics base tail sizing
    b777300ER compares to TASOPT 777-300ER
    """
    if aircraft == 'optimal737':
       #weights to compare for TASOPT 737-800 run...sizing w/out specifying tail volume
        print "WEIGHT DIFFERENCES"
        print "\n"
        print "Total Fuel Weight Percent Diff: %s" % compute_diff(mag(sol('W_{f_{total}}').to('lbf')), 45057.0)
        print "\n"
        print "Total Aircraft Weight Percent Diff: %s" % compute_diff(mag(sol('W_{total}').to('lbf')), 166502.0)
        print "\n"
        print "Engine Weight Percent Diff: %s" % compute_diff(2.*mag(sol('W_{engsys}').to('lbf')), 11632.6)
        print "\n"
        print "Fuselage Weight Percent Diff: %s" % compute_diff(mag(sol('W_{fuse}').to('lbf')), 35641.6)
        print "\n"
        print "Payload Weight Percent Diff: %s" % compute_diff(mag(sol('W_{payload}').to('lbf')), 38715.5)
        print "\n"
        print "VT Weight Percent Diff: %s" % compute_diff(mag(sol('W_{VT}').to('lbf')), 1477.2)
        print "\n"
        print "HT Weight Percent Diff: %s" % compute_diff(mag(sol('W_{HT}').to('lbf')), 1463.8)
        print "\n"
        print "Wing Weight Percent Diff: %s" % compute_diff(mag(sol('W_{wing}').to('lbf')), 21691.8)

        #wing value to compare
        print "\n\n\n"
        print "WING DIFFERENCES"
        print "\n"
        print "Wing Span Percent Diff: %s" % compute_diff(mag(sol('b').to('ft')), 113.586)
##        print "\n"
##        print "Wing Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR')), 10.1)
        print "\n"
        print "Wing Area Percent Diff: %s" % compute_diff(mag(sol('S').to('ft^2')), 1277.41)
        

        #HT values to compare
        print "\n\n\n"
        print "HORIZONTAL TAIL DIFFERENCES"
##        print "\n"
##        print "HT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR_{ht}')), 6)
        print "\n"
        print "HT Area Percent Diff: %s" % compute_diff(mag(sol('S_{ht}').to('ft^2')), 313.86  )
##        print "\n"
##        print "HT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{h}')), 1.450) 
        

        #VT values to compare
        print "\n\n\n"
        print "VERTICAL TAIL DIFFERENCES"
        print "\n"
##        print "VT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('A_{vt}')), 2.0)
        print "\n"
        print "VT Span Percent Diff: %s" % compute_diff(mag(sol('b_{vt}').to('ft')), 22.99)
        print "\n"
        print "VT Area Percent Diff: %s" % compute_diff(mag(sol('S_{vt}').to('ft^2')), 264.29)
##        print "\n"
##        print "VT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{vt}')), .1)


        #drag values to compare
        print "\n\n\n"
        print "CRUISE SEGMENT 1 DRAG DIFFERENCES"
        print "\n"
        print "L/D Percent Diff: %s" % compute_diff(mag(sol('L/D')[0]), 17.252)
        print "\n"
        print "Overall Cd Percent Diff: %s" % compute_diff(mag(sol('C_D')[0]), 0.03292)
        print "\n"
        print "Nacelle Cd Percent Diff: %s" % compute_diff(mag(sol('C_{d_nacelle}')[0]),0.00191)
        print "\n"
        print "HT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_h}')[0]), 0.00175*sol('S')/sol('S_{ht}'))
        print "\n"
        print "Fuselage Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{fuse}}')[0]), 0.00801)
        print "\n"
        print "VT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{vis}}')[0]), 0.00153*sol('S')/sol('S_{vt}'))
        print "\n"
        print "Induced Drag Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{i_w}}')[0]), 0.01112)
        print "\n"
        print "Wing Profile Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{p_w}}')[0]), 0.00833)
        

        #Fuselage values to compare
        print "\n\n\n"
        print "FUSELAGE DIFFERENCES"
        print "\n"
        print "Weight of HB material: %s "  % compute_diff(mag(sol('W_{hbend}')), 2021.3)
        print "\n"
        print "Weight of VB material: %s "  % compute_diff(mag(sol('W_{vbend}')), 1210.5)
    
def compute_diff(sp, tasopt):
    """
    Method to actually compute the percent difference
    """
    diff = 100*(sp-tasopt)/tasopt
    
    return diff
