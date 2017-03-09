"""
Method to compute and print percent differences between SP D8 model and TASOPT
"""

from gpkit.small_scripts import mag

def percent_diff(sol, version, Nclimb):
    """
    Method to compute and print percent differences between SP D8 model and TASOPT

    INPUTS
    ------
    sol: solution from D8 SP model
    version: either 0, 1, or 2....corresponds to comparing to the D8.0, D8.1, or D8.2
    """

    if version == 800:
       #weights to compare for TASOPT 737-800 run
        print "WEIGHT DIFFERENCES"
        print "\n"
        print "Total Fuel Weight Percent Diff: %s" % compute_diff(mag(sol('W_{f_{total}}').to('lbf')), 47487.3)
        print "\n"
        print "Total Aircraft Weight Percent Diff: %s" % compute_diff(mag(sol('W_{total}').to('lbf')), 174979.1)
        print "\n"
        print "Engine Weight Percent Diff: %s" % compute_diff(mag(sol('W_{engine}').to('lbf')), 12280.0)
        print "\n"
        print "Fuselage Weight Percent Diff: %s" % compute_diff(mag(sol('W_{fuse}').to('lbf')), 37025.2)
        print "\n"
        print "Payload Weight Percent Diff: %s" % compute_diff(mag(sol('W_{payload}').to('lbf')), 38715.5)
        print "\n"
        print "VT Weight Percent Diff: %s" % compute_diff(mag(sol('W_{VT}').to('lbf')), 1764.1)
        print "\n"
        print "HT Weight Percent Diff: %s" % compute_diff(mag(sol('W_{HT}').to('lbf')), 2616.1)
        print "\n"
        print "Wing Weight Percent Diff: %s" % compute_diff(mag(sol('W_{wing}').to('lbf')), 23717.3)

        #wing value to compare
        print "\n\n\n"
        print "WING DIFFERENCES"
        print "\n"
        print "Wing Span Percent Diff: %s" % compute_diff(mag(sol('b').to('ft')), 116.427)
##        print "\n"
##        print "Wing Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR')), 10.1)
        print "\n"
        print "Wing Area Percent Diff: %s" % compute_diff(mag(sol('S').to('ft^2')), 1342.10)
        

        #HT values to compare
        print "\n\n\n"
        print "HORIZONTAL TAIL DIFFERENCES"
##        print "\n"
##        print "HT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR_h')), 6)
        print "\n"
        print "HT Area Percent Diff: %s" % compute_diff(mag(sol('S_h').to('ft^2')), 456.86  )
##        print "\n"
##        print "HT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{h}')), 1.450) 
        

        #VT values to compare
        print "\n\n\n"
        print "VERTICAL TAIL DIFFERENCES"
        print "\n"
##        print "VT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('A_{vt}')), 2.0)
        print "\n"
        print "VT Span Percent Diff: %s" % compute_diff(mag(sol('b_{vt}').to('ft')), 24.39)
        print "\n"
        print "VT Area Percent Diff: %s" % compute_diff(mag(sol('S_{vt}').to('ft^2')), 297.50)
##        print "\n"
##        print "VT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{vt}')), .1)


        #drag values to compare
        print "\n\n\n"
        print "CRUISE SEGMENT 1 DRAG DIFFERENCES"
        print "\n"
        print "Overall Cd Percent Diff: %s" % compute_diff(mag(sol('C_D')['C_D_Mission, CruiseSegment, CruiseP, AircraftP.1'][0]), 0.03304)
        print "\n"
        print "Nacelle Cd Percent Diff: %s" % compute_diff(mag(sol('C_{d_nacelle}')['C_{d_nacelle}_Mission, CruiseSegment, CruiseP, AircraftP.1'][0]), 0.00191)
        print "\n"
        print "HT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_h}')['C_{D_h}_Mission, CruiseSegment, CruiseP, AircraftP.1, HorizontalTailPerformance.1'][0]), 0.00239)
        print "\n"
        print "Fuselage Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{fuse}}')['C_{D_{fuse}}_Mission, CruiseSegment, CruiseP, AircraftP.1, FuselagePerformance.1'][0]), 0.00762)
        print "\n"
        print "VT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{vis}}')['C_{D_{vis}}_Mission, CruiseSegment, CruiseP, AircraftP.1, VerticalTailPerformance.1'][0]), 0.00163)
        print "\n"
        print "Wing Cd (w/Induced Drag) Percent Diff: %s" % compute_diff(mag(sol('C_{d_w}')['C_{d_w}_Mission, CruiseSegment, CruiseP, AircraftP.1, WingPerformance.1'][0]), 0.00833)
        print "\n"
        print "Wing Profile Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{p_w}}')['C_{D_{p_w}}_Mission, CruiseSegment, CruiseP, AircraftP.1, WingPerformance.1'][0]), 0.00833)
        
        #crusie TSFC
        print "\n\n\n"
        print "CRUISE SEGMENT 1 TSFC DIFFERENCES"
        print "\n"
        print "Initial Cruise TSFC Percent Diff: %s" % compute_diff(mag(sol('TSFC')[Nclimb]), 0.64009)
        

        #Fuselage values to compare
        print "\n\n\n"
        print "FUSELAGE DIFFERENCES"
        print "\n"
        print

    if version == 801:
       #weights to compare for TASOPT 737-800 run...sizing w/out specifying tail volume
        print "WEIGHT DIFFERENCES"
        print "\n"
        print "Total Fuel Weight Percent Diff: %s" % compute_diff(mag(sol('W_{f_{total}}').to('lbf')), 45057.0)
        print "\n"
        print "Total Aircraft Weight Percent Diff: %s" % compute_diff(mag(sol('W_{total}').to('lbf')), 166502.0)
        print "\n"
        print "Engine Weight Percent Diff: %s" % compute_diff(mag(sol('W_{engine}').to('lbf')), 7467.9)
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
##        print "HT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR_h')), 6)
        print "\n"
        print "HT Area Percent Diff: %s" % compute_diff(mag(sol('S_h').to('ft^2')), 313.86  )
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
        print "Overall Cd Percent Diff: %s" % compute_diff(mag(sol('C_D')['C_D_Mission, CruiseSegment, CruiseP, AircraftP.1'][0]), 0.03292)
        print "\n"
        print "Nacelle Cd Percent Diff: %s" % compute_diff(mag(sol('C_{d_nacelle}')['C_{d_nacelle}_Mission, CruiseSegment, CruiseP, AircraftP.1'][0]),0.00191)
        print "\n"
        print "HT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_h}')['C_{D_h}_Mission, CruiseSegment, CruiseP, AircraftP.1, HorizontalTailPerformance.1'][0]), 0.00175)
        print "\n"
        print "Fuselage Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{fuse}}')['C_{D_{fuse}}_Mission, CruiseSegment, CruiseP, AircraftP.1, FuselagePerformance.1'][0]), 0.00801)
        print "\n"
        print "VT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{vis}}')['C_{D_{vis}}_Mission, CruiseSegment, CruiseP, AircraftP.1, VerticalTailPerformance.1'][0]), 0.00153)
        print "\n"
        print "Wing Cd (w/Induced Drag) Percent Diff: %s" % compute_diff(mag(sol('C_{d_w}')['C_{d_w}_Mission, CruiseSegment, CruiseP, AircraftP.1, WingPerformance.1'][0]), 0.00861)
        print "\n"
        print "Wing Profile Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{p_w}}')['C_{D_{p_w}}_Mission, CruiseSegment, CruiseP, AircraftP.1, WingPerformance.1'][0]), 0.00861)
        
        #crusie TSFC
        print "\n\n\n"
        print "CRUISE SEGMENT 1 TSFC DIFFERENCES"
        print "\n"
        print "Initial Cruise TSFC Percent Diff: %s" % compute_diff(mag(sol('TSFC')[Nclimb]), 0.64030)
        

        #Fuselage values to compare
        print "\n\n\n"
        print "FUSELAGE DIFFERENCES"
        print "\n"
        print
    
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
##        print "\n"
##        print "Wing Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR')), 15.7490)
        print "\n"
        print "Wing Area Percent Diff: %s" % compute_diff(mag(sol('S').to('ft^2')), 1245.27)
        

        #HT values to compare
        print "\n\n\n"
        print "HORIZONTAL TAIL DIFFERENCES"
##        print "\n"
##        print "HT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR_h')), 12.0000)
        print "\n"
        print "HT Area Percent Diff: %s" % compute_diff(mag(sol('S_h').to('ft^2')), 243.76  )
##        print "\n"
##        print "HT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{h}')), .913) 
        

        #VT values to compare
        print "\n\n\n"
        print "VERTICAL TAIL DIFFERENCES"
        print "\n"
##        print "VT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('A_{vt}')), 2.2000)
        print "\n"
        print "VT Span Percent Diff: %s" % compute_diff(mag(sol('b_{vt}').to('ft')), 11.59)
        print "\n"
        print "VT Area Percent Diff: %s" % compute_diff(mag(sol('S_{vt}').to('ft^2')), 122.03)
##        print "\n"
##        print "VT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{vt}')), .03)


        #drag values to compare
        print "\n\n\n"
        print "CRUISE SEGMENT 1 DRAG DIFFERENCES"
        print "\n"
        print "Overall Cd Percent Diff: %s" % compute_diff(mag(sol('C_D')['C_D_Mission, CruiseSegment, CruiseP, AircraftP.1'][0]), 0.03212)
        print "\n"
        print "Nacelle Cd Percent Diff: %s" % compute_diff(mag(sol('C_{d_nacelle}')['C_{d_nacelle}_Mission, CruiseSegment, CruiseP, AircraftP.1'][0]), 0.00054)
        print "\n"
        print "HT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_h}')['C_{D_h}_Mission, CruiseSegment, CruiseP, AircraftP.1, HorizontalTailPerformance.1'][0]), 0.00227)
        print "\n"
        print "Fuselage Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{fuse}}')['C_{D_{fuse}}_Mission, CruiseSegment, CruiseP, AircraftP.1, FuselagePerformance.1'][0]), 0.00866)
        print "\n"
        print "VT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{vis}}')['C_{D_{vis}}_Mission, CruiseSegment, CruiseP, AircraftP.1, VerticalTailPerformance.1'][0]), 0.00089)

        #crusie TSFC
        print "\n\n\n"
        print "CRUISE SEGMENT 1 TSFC DIFFERENCES"
        print "\n"
        print "Initial Cruise TSFC Percent Diff: %s" % compute_diff(mag(sol('TSFC')[Nclimb]), 0.49971)
        

        #Fuselage values to compare
        print "\n\n\n"
        print "FUSELAGE DIFFERENCES"
        print "\n"
        print

def compute_diff(sp, tasopt):
    """
    Method to actually compute the percent difference
    """
    diff = 100*(sp-tasopt)/sp
    
    return diff
