"""
Method to compute and print percent differences between SP D8 model and TASOPT
"""

from gpkit.small_scripts import mag

def percent_diff(sol,aircraft,Nclimb):
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
    if aircraft == 'M072_737':
       #weights to compare for TASOPT 777-300ER run
        print "WEIGHT DIFFERENCES"
        print "\n"
        print "Total Fuel Weight Percent Diff: %s" % compute_diff(mag(sol('W_{f_{total}}').to('lbf')), 40530.2)
        print "\n"
        print "Total Aircraft Weight Percent Diff: %s" % compute_diff(mag(sol('W_{total}').to('lbf')), 155258.4)
        print "\n"
        print "Engine Weight Percent Diff: %s" % compute_diff(2.*mag(sol('W_{engsys}').to('lbf')), 11601.5)
        print "\n"
        print "Fuselage Weight Percent Diff: %s" % compute_diff(mag(sol('W_{fuse}').to('lbf')), 35360.9)
        print "\n"
        print "Payload Weight Percent Diff: %s" % compute_diff(mag(sol('W_{payload}').to('lbf')), 38715.5)
        print "\n"
        print "VT Weight Percent Diff: %s" % compute_diff(mag(sol('W_{VT}').to('lbf')), 1383.8)
        print "\n"
        print "HT Weight Percent Diff: %s" % compute_diff(mag(sol('W_{HT}').to('lbf')), 1424.0)
        print "\n"
        print "Wing Weight Percent Diff: %s" % compute_diff(mag(sol('W_{wing}').to('lbf')), 16150.7)

        #wing value to compare
        print "\n\n\n"
        print "WING DIFFERENCES"
        print "\n"
        print "Wing Span Percent Diff: %s" % compute_diff(mag(sol('b').to('ft')), 118.500)
##        print "\n"
##        print "Wing Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR')), 10.1)
        print "\n"
        print "Wing Area Percent Diff: %s" % compute_diff(mag(sol('S').to('ft^2')), 1208.24)
        

        #HT values to compare
        print "\n\n\n"
        print "HORIZONTAL TAIL DIFFERENCES"
##        print "\n"
##        print "HT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR_{ht}')), 6)
        print "\n"
        print "HT Area Percent Diff: %s" % compute_diff(mag(sol('S_{ht}').to('ft^2')), 367.39)
##        print "\n"
##        print "HT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{h}')), 1.450) 
        

        #VT values to compare
        print "\n\n\n"
        print "VERTICAL TAIL DIFFERENCES"
        print "\n"
##        print "VT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('A_{vt}')), 2.0)
        print "\n"
        print "VT Span Percent Diff: %s" % compute_diff(mag(sol('b_{vt}').to('ft')), 22.50)
        print "\n"
        print "VT Area Percent Diff: %s" % compute_diff(mag(sol('S_{vt}').to('ft^2')), 253.03)
##        print "\n"
##        print "VT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{vt}')), .1)


        #drag values to compare
        print "\n\n\n"
        print "CRUISE SEGMENT 1 DRAG DIFFERENCES"
        print "\n"
        print "Overall Cd Percent Diff: %s" % compute_diff(mag(sol('C_D')['C_D_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 0.03581)
        print "\n"
        # print "L/D Percent Diff: %s" % compute_diff(mag(sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 17.252)
        # print "\n"
        print "Nacelle Cd Percent Diff: %s" % compute_diff(mag(sol('C_{d_nacelle}')['C_{d_nacelle}_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 0.00181)
        print "\n"
        print "HT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_h}')['C_{D_h}_Mission/CruiseSegment/CruiseP/AircraftP/HorizontalTailPerformance'][0]), 0.00154*sol('S')/sol('S_{ht}'))
        print "\n"
        print "Fuselage Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{fuse}}')['C_{D_{fuse}}_Mission/CruiseSegment/CruiseP/AircraftP/FuselagePerformance'][0]), 0.00234)
        print "\n"
        print "VT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{vis}}')['C_{D_{vis}}_Mission/CruiseSegment/CruiseP/AircraftP/VerticalTailPerformance'][0]), 0.00157*sol('S')/sol('S_{vt}'))
        print "\n"
        print "Induced Drag Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{i_w}}')['C_{D_{i_w}}_Mission/CruiseSegment/CruiseP/AircraftP/WingPerformance'][0]), 0.01242)
        print "\n"
        print "Wing Profile Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{p_w}}')['C_{D_{p_w}}_Mission/CruiseSegment/CruiseP/AircraftP/WingPerformance'][0]),0.00915)
        
        #crusie TSFC
        print "\n\n\n"
        print "CRUISE SEGMENT 1 TSFC DIFFERENCES"
        print "\n"
        print "Initial Cruise TSFC Percent Diff: %s" % compute_diff(mag(sol('TSFC')[Nclimb]), 0.58399)
        

        #Fuselage values to compare
        print "\n\n\n"
        print "FUSELAGE DIFFERENCES"
        print "\n"
        # print "Weight of HB material: %s "  % compute_diff(mag(sol('W_{hbend}')), 2021.3)
        # print "\n"
        # print "Weight of VB material: %s "  % compute_diff(mag(sol('W_{vbend}')), 1210.5)


    if aircraft == 'D12':
       #weights to compare for TASOPT 777-300ER run
        print "WEIGHT DIFFERENCES"
        print "\n"
        print "Total Fuel Weight Percent Diff: %s" % compute_diff(mag(sol('W_{f_{total}}').to('lbf')), 152030.4)
        print "\n"
        print "Total Aircraft Weight Percent Diff: %s" % compute_diff(mag(sol('W_{total}').to('lbf')), 528633.7)
        print "\n"
        print "Engine Weight Percent Diff: %s" % compute_diff(2.*mag(sol('W_{engsys}').to('lbf')), 35537.6)
        print "\n"
        print "Fuselage Weight Percent Diff: %s" % compute_diff(mag(sol('W_{fuse}').to('lbf')), 107426.3)
        print "\n"
        print "Payload Weight Percent Diff: %s" % compute_diff(mag(sol('W_{payload}').to('lbf')), 115046.0)
        print "\n"
        print "VT Weight Percent Diff: %s" % compute_diff(mag(sol('W_{VT}').to('lbf')), 2368.3)
        print "\n"
        print "HT Weight Percent Diff: %s" % compute_diff(mag(sol('W_{HT}').to('lbf')), 9058.3)
        print "\n"
        print "Wing Weight Percent Diff: %s" % compute_diff(mag(sol('W_{wing}').to('lbf')), 91144.3)

        #wing value to compare
        print "\n\n\n"
        print "WING DIFFERENCES"
        print "\n"
        print "Wing Span Percent Diff: %s" % compute_diff(mag(sol('b').to('ft')), 199.748)
##        print "\n"
##        print "Wing Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR')), 10.1)
        print "\n"
        print "Wing Area Percent Diff: %s" % compute_diff(mag(sol('S').to('ft^2')), 3989.93)
        

        #HT values to compare
        print "\n\n\n"
        print "HORIZONTAL TAIL DIFFERENCES"
##        print "\n"
##        print "HT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR_{ht}')), 6)
        print "\n"
        print "HT Area Percent Diff: %s" % compute_diff(mag(sol('S_{ht}').to('ft^2')), 689.07  )
##        print "\n"
##        print "HT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{h}')), 1.450) 
        

        #VT values to compare
        print "\n\n\n"
        print "VERTICAL TAIL DIFFERENCES"
        print "\n"
##        print "VT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('A_{vt}')), 2.0)
        print "\n"
        print "VT Span Percent Diff: %s" % compute_diff(mag(sol('b_{vt}').to('ft')), 21.36)
        print "\n"
        print "VT Area Percent Diff: %s" % compute_diff(mag(sol('S_{vt}').to('ft^2')), 456.15)
##        print "\n"
##        print "VT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{vt}')), .1)


        #drag values to compare
        print "\n\n\n"
        print "CRUISE SEGMENT 1 DRAG DIFFERENCES"
        print "\n"
        print "Overall Cd Percent Diff: %s" % compute_diff(mag(sol('C_D')['C_D_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 0.02428)
        print "\n"
        # print "L/D Percent Diff: %s" % compute_diff(mag(sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 17.252)
        # print "\n"
        print "Nacelle Cd Percent Diff: %s" % compute_diff(mag(sol('C_{d_nacelle}')['C_{d_nacelle}_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 0.00062)
        print "\n"
        print "HT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_h}')['C_{D_h}_Mission/CruiseSegment/CruiseP/AircraftP/HorizontalTailPerformance'][0]), 0.00154*sol('S')/sol('S_{ht}'))
        print "\n"
        print "Fuselage Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{fuse}}')['C_{D_{fuse}}_Mission/CruiseSegment/CruiseP/AircraftP/FuselagePerformance'][0]), 0.00521)
        print "\n"
        print "VT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{vis}}')['C_{D_{vis}}_Mission/CruiseSegment/CruiseP/AircraftP/VerticalTailPerformance'][0]), 0.00093*sol('S')/sol('S_{vt}'))
        print "\n"
        print "Induced Drag Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{i_w}}')['C_{D_{i_w}}_Mission/CruiseSegment/CruiseP/AircraftP/WingPerformance'][0]), 0.00899)
        print "\n"
        print "Wing Profile Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{p_w}}')['C_{D_{p_w}}_Mission/CruiseSegment/CruiseP/AircraftP/WingPerformance'][0]), 0.00723)
        
        #crusie TSFC
        print "\n\n\n"
        print "CRUISE SEGMENT 1 TSFC DIFFERENCES"
        print "\n"
        print "Initial Cruise TSFC Percent Diff: %s" % compute_diff(mag(sol('TSFC')[Nclimb]), 0.54649)
        

        #Fuselage values to compare
        print "\n\n\n"
        print "FUSELAGE DIFFERENCES"
        print "\n"
        # print "Weight of HB material: %s "  % compute_diff(mag(sol('W_{hbend}')), 2021.3)
        # print "\n"
        # print "Weight of VB material: %s "  % compute_diff(mag(sol('W_{vbend}')), 1210.5)

    if aircraft == 'b777300ER':
       #weights to compare for TASOPT 777-300ER run
        print "WEIGHT DIFFERENCES"
        print "\n"
        print "Total Fuel Weight Percent Diff: %s" % compute_diff(mag(sol('W_{f_{total}}').to('lbf')), 212236.0)
        print "\n"
        print "Total Aircraft Weight Percent Diff: %s" % compute_diff(mag(sol('W_{total}').to('lbf')), 625008.3)
        print "\n"
        print "Engine Weight Percent Diff: %s" % compute_diff(2.*mag(sol('W_{engsys}').to('lbf')), 53717.8)
        print "\n"
        print "Fuselage Weight Percent Diff: %s" % compute_diff(mag(sol('W_{fuse}').to('lbf')), 107426.3)
        print "\n"
        print "Payload Weight Percent Diff: %s" % compute_diff(mag(sol('W_{payload}').to('lbf')), 103541.4)
        print "\n"
        print "VT Weight Percent Diff: %s" % compute_diff(mag(sol('W_{VT}').to('lbf')), 4872.7)
        print "\n"
        print "HT Weight Percent Diff: %s" % compute_diff(mag(sol('W_{HT}').to('lbf')), 5844.6)
        print "\n"
        print "Wing Weight Percent Diff: %s" % compute_diff(mag(sol('W_{wing}').to('lbf')), 99869.0)

        #wing value to compare
        print "\n\n\n"
        print "WING DIFFERENCES"
        print "\n"
        print "Wing Span Percent Diff: %s" % compute_diff(mag(sol('b').to('ft')), 197.733)
##        print "\n"
##        print "Wing Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR')), 10.1)
        print "\n"
        print "Wing Area Percent Diff: %s" % compute_diff(mag(sol('S').to('ft^2')), 4624.26)
        

        #HT values to compare
        print "\n\n\n"
        print "HORIZONTAL TAIL DIFFERENCES"
##        print "\n"
##        print "HT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR_{ht}')), 6)
        print "\n"
        print "HT Area Percent Diff: %s" % compute_diff(mag(sol('S_{ht}').to('ft^2')), 969.73  )
##        print "\n"
##        print "HT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{h}')), 1.450) 
        

        #VT values to compare
        print "\n\n\n"
        print "VERTICAL TAIL DIFFERENCES"
        print "\n"
##        print "VT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('A_{vt}')), 2.0)
        print "\n"
        print "VT Span Percent Diff: %s" % compute_diff(mag(sol('b_{vt}').to('ft')), 38.04)
        print "\n"
        print "VT Area Percent Diff: %s" % compute_diff(mag(sol('S_{vt}').to('ft^2')), 615.92)
##        print "\n"
##        print "VT Volume Coefficient Percent Diff: %s" % compute_diff(mag(sol('V_{vt}')), .1)


        #drag values to compare
        print "\n\n\n"
        print "CRUISE SEGMENT 1 DRAG DIFFERENCES"
        print "\n"
        print "Overall Cd Percent Diff: %s" % compute_diff(mag(sol('C_D')['C_D_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 0.02539)
        print "\n"
        # print "L/D Percent Diff: %s" % compute_diff(mag(sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 17.252)
        # print "\n"
        print "Nacelle Cd Percent Diff: %s" % compute_diff(mag(sol('C_{d_nacelle}')['C_{d_nacelle}_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 0.00151)
        print "\n"
        print "HT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_h}')['C_{D_h}_Mission/CruiseSegment/CruiseP/AircraftP/HorizontalTailPerformance'][0]), 0.00126*sol('S')/sol('S_{ht}'))
        print "\n"
        print "Fuselage Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{fuse}}')['C_{D_{fuse}}_Mission/CruiseSegment/CruiseP/AircraftP/FuselagePerformance'][0]), 0.00532)
        print "\n"
        print "VT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{vis}}')['C_{D_{vis}}_Mission/CruiseSegment/CruiseP/AircraftP/VerticalTailPerformance'][0]), 0.00093*sol('S')/sol('S_{vt}'))
        print "\n"
        print "Induced Drag Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{i_w}}')['C_{D_{i_w}}_Mission/CruiseSegment/CruiseP/AircraftP/WingPerformance'][0]), 0.00942)
        print "\n"
        print "Wing Profile Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{p_w}}')['C_{D_{p_w}}_Mission/CruiseSegment/CruiseP/AircraftP/WingPerformance'][0]), 0.00695)
        
        #crusie TSFC
        print "\n\n\n"
        print "CRUISE SEGMENT 1 TSFC DIFFERENCES"
        print "\n"
        print "Initial Cruise TSFC Percent Diff: %s" % compute_diff(mag(sol('TSFC')[Nclimb]), 0.56581)
        

        #Fuselage values to compare
        print "\n\n\n"
        print "FUSELAGE DIFFERENCES"
        print "\n"
        # print "Weight of HB material: %s "  % compute_diff(mag(sol('W_{hbend}')), 2021.3)
        # print "\n"
        # print "Weight of VB material: %s "  % compute_diff(mag(sol('W_{vbend}')), 1210.5)

    if aircraft == 'b737800':
       #weights to compare for TASOPT 737-800 run
        print "WEIGHT DIFFERENCES"
        print "\n"
        print "Total Fuel Weight Percent Diff: %s" % compute_diff(mag(sol('W_{f_{total}}').to('lbf')), 45057.0)
        print "\n"
        print "Total Aircraft Weight Percent Diff: %s" % compute_diff(mag(sol('W_{total}').to('lbf')), 174979.1)
        print "\n"
        print "Engine Weight Percent Diff: %s" % compute_diff(2.*mag(sol('W_{engsys}').to('lbf')), 12280.0)
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
##        print "HT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR_{ht}')), 6)
        print "\n"
        print "HT Area Percent Diff: %s" % compute_diff(mag(sol('S_{ht}').to('ft^2')), 456.86  )
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
        print "Overall Cd Percent Diff: %s" % compute_diff(mag(sol('C_D')['C_D_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 0.03304)
        print "\n"
        # print "L/D Percent Diff: %s" % compute_diff(mag(sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 17.252)
        # print "\n"
        print "Nacelle Cd Percent Diff: %s" % compute_diff(mag(sol('C_{d_nacelle}')['C_{d_nacelle}_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 0.00191)
        print "\n"
        print "HT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_h}')['C_{D_h}_Mission/CruiseSegment/CruiseP/AircraftP/HorizontalTailPerformance'][0]), 0.00239*sol('S')/sol('S_{ht}'))
        print "\n"
        print "Fuselage Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{fuse}}')['C_{D_{fuse}}_Mission/CruiseSegment/CruiseP/AircraftP/FuselagePerformance'][0]), 0.00762)
        print "\n"
        print "VT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{vis}}')['C_{D_{vis}}_Mission/CruiseSegment/CruiseP/AircraftP/VerticalTailPerformance'][0]), 0.00163*sol('S')/sol('S_{vt}'))
        print "\n"
        print "Induced Drag Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{i_w}}')['C_{D_{i_w}}_Mission/CruiseSegment/CruiseP/AircraftP/WingPerformance'][0]), 0.01117)
        print "\n"
        print "Wing Profile Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{p_w}}')['C_{D_{p_w}}_Mission/CruiseSegment/CruiseP/AircraftP/WingPerformance'][0]), 0.00833)
        
        #crusie TSFC
        print "\n\n\n"
        print "CRUISE SEGMENT 1 TSFC DIFFERENCES"
        print "\n"
        print "Initial Cruise TSFC Percent Diff: %s" % compute_diff(mag(sol('TSFC')[Nclimb]), 0.64009)
        

        #Fuselage values to compare
        print "\n\n\n"
        print "FUSELAGE DIFFERENCES"
        print "\n"
        # print "Weight of HB material: %s "  % compute_diff(mag(sol('W_{hbend}')), 2021.3)
        # print "\n"
        # print "Weight of VB material: %s "  % compute_diff(mag(sol('W_{vbend}')), 1210.5)

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
        print "L/D Percent Diff: %s" % compute_diff(mag(sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 17.252)
        print "\n"
        print "Overall Cd Percent Diff: %s" % compute_diff(mag(sol('C_D')['C_D_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 0.03292)
        print "\n"
        print "Nacelle Cd Percent Diff: %s" % compute_diff(mag(sol('C_{d_nacelle}')['C_{d_nacelle}_Mission/CruiseSegment/CruiseP/AircraftP'][0]),0.00191)
        print "\n"
        print "HT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_h}')['C_{D_h}_Mission/CruiseSegment/CruiseP/AircraftP/HorizontalTailPerformance'][0]), 0.00175*sol('S')/sol('S_{ht}'))
        print "\n"
        print "Fuselage Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{fuse}}')['C_{D_{fuse}}_Mission/CruiseSegment/CruiseP/AircraftP/FuselagePerformance'][0]), 0.00801)
        print "\n"
        print "VT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{vis}}')['C_{D_{vis}}_Mission/CruiseSegment/CruiseP/AircraftP/VerticalTailPerformance'][0]), 0.00153*sol('S')/sol('S_{vt}'))
        print "\n"
        print "Induced Drag Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{i_w}}')['C_{D_{i_w}}_Mission/CruiseSegment/CruiseP/AircraftP/WingPerformance'][0]), 0.01112)
        print "\n"
        print "Wing Profile Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{p_w}}')['C_{D_{p_w}}_Mission/CruiseSegment/CruiseP/AircraftP/WingPerformance'][0]), 0.00861)
        
        #crusie TSFC
        print "\n\n\n"
        print "CRUISE SEGMENT 1 TSFC DIFFERENCES"
        print "\n"
        print "Initial Cruise TSFC Percent Diff: %s" % compute_diff(mag(sol('TSFC')[Nclimb]), 0.64030)
        

        #Fuselage values to compare
        print "\n\n\n"
        print "FUSELAGE DIFFERENCES"
        print "\n"
        print "Weight of HB material: %s "  % compute_diff(mag(sol('W_{hbend}')), 2021.3)
        print "\n"
        print "Weight of VB material: %s "  % compute_diff(mag(sol('W_{vbend}')), 1210.5)
    
    if aircraft == 'D82':
        #weights to compare for D8.2
        print "WEIGHT DIFFERENCES"
        print "\n"
        print "Total Fuel Weight Percent Diff: %s" % compute_diff(mag(sol('W_{f_{total}}').to('lbf')), 26959.1)
        print "\n"
        print "Total Aircraft Weight Percent Diff: %s" % compute_diff(mag(sol('W_{total}').to('lbf')), 134758.2)
        print "\n"
        print "Engine Weight Percent Diff: %s" % compute_diff(2.*mag(sol('W_{engsys}').to('lbf')), 10502.8)
        print "\n"
        print "Fuselage Weight Percent Diff: %s" % compute_diff(mag(sol('W_{fuse}').to('lbf')), 30487.4)
        print "\n"
        print "Payload Weight Percent Diff: %s" % compute_diff(mag(sol('W_{payload}').to('lbf')), 38715.5)
        print "\n"
        print "VT Weight Percent Diff: %s" % compute_diff(2.*mag(sol('W_{VT}').to('lbf')), 376.4)
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
##        print "HT Aspect Ratio Percent Diff: %s" % compute_diff(mag(sol('AR_{ht}')), 12.0000)
        print "\n"
        print "HT Area Percent Diff: %s" % compute_diff(mag(sol('S_{ht}').to('ft^2')), 243.76  )
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
        print "Overall Cd Percent Diff: %s" % compute_diff(mag(sol('C_D')['C_D_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 0.03212)
        print "\n"
        print "L/D Percent Diff: %s" % compute_diff(mag(sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 21.477)
        print "\n"
        print "Nacelle Cd Percent Diff: %s" % compute_diff(mag(sol('C_{d_nacelle}')['C_{d_nacelle}_Mission/CruiseSegment/CruiseP/AircraftP'][0]), 0.00054)
        print "\n"
        print "HT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_h}')['C_{D_h}_Mission/CruiseSegment/CruiseP/AircraftP/HorizontalTailPerformance'][0]), 0.00227*sol('S')/sol('S_{ht}'))
        print "\n"
        print "Fuselage Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{fuse}}')['C_{D_{fuse}}_Mission/CruiseSegment/CruiseP/AircraftP/FuselagePerformance'][0]), 0.00866)
        print "\n"
        print "VT Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{vis}}')['C_{D_{vis}}_Mission/CruiseSegment/CruiseP/AircraftP/VerticalTailPerformance'][0]), 0.00089*sol('S')/sol('S_{vt}'))
        print "\n"
        print "Wing Cd percent Diff: %s" % compute_diff(mag(sol('C_{D_{p_w}}')['C_{D_{p_w}}_Mission/CruiseSegment/CruiseP/AircraftP/WingPerformance'][0]),0.00947)
        print "Induced Drag Cd Percent Diff: %s" % compute_diff(mag(sol('C_{D_{i_w}}')['C_{D_{i_w}}_Mission/CruiseSegment/CruiseP/AircraftP/WingPerformance'][0]), 0.01063)
        
        #crusie TSFC
        print "\n\n\n"
        print "CRUISE SEGMENT 1 TSFC DIFFERENCES"
        print "\n"
        print "Initial Cruise TSFC Percent Diff: %s" % compute_diff(mag(sol('TSFC')[Nclimb]), 0.49971)
        

        #Fuselage values to compare
        print "\n\n\n"
        print "FUSELAGE DIFFERENCES"
        print "\n"
        print "\n"
        print "Weight of HB material: %s "  % compute_diff(mag(sol('W_{hbend}')), 359.7)
        print "\n"
        print "Weight of VB material: %s "  % compute_diff(mag(sol('W_{vbend}')), 5.5)

def compute_diff(sp, tasopt):
    """
    Method to actually compute the percent difference
    """
    diff = 100*(sp-tasopt)/tasopt
    
    return diff
