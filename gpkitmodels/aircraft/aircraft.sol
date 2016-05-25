Beginning signomial solve.
Allowing negative dual variable(s) as small as -1.23498217014e-11.
Solving took 6 GP solves and 3.27 seconds.

Cost
----
 5.164e+04 [N] 

Free Variables
--------------
                         Aircraft | 
                             AR_w : 9.298 Wing aspect ratio
                              C_D : 0.02911 Drag coefficient
                              C_L : 0.5616 Lift coefficient
                          C_{L_w} : 0.5616 Lift coefficient, wing
                       C_{L_{aw}} : 5.616 Lift curve slope, wing
                                D : 2.909e+04  [N] Total aircraft drag (cruise)
                         D_{fuse} : 1.038e+04  [N] Fuselage drag
                           D_{ht} : 1011  [N] Horizontal tail drag
                           D_{vt} : 944.6  [N] Vertical tail drag
                         D_{wing} : 1.675e+04  [N] Wing drag
                              L_w : 5.612e+05  [N] Wing lift
                      L_{v_{max}} : 6.399e+05  [N] Maximum load for structural sizing
                              S_w : 98  [m**2] Wing reference area
                           V_{TO} : 73.38  [m/s] Takeoff speed
                       V_{\infty} : 231.7  [m/s] Cruise velocity
                                W : 5.612e+05  [N] Total aircraft weight
                         W_{fuel} : 5.164e+04  [N] Fuel weight
                         W_{fuse} : 1.616e+05  [N] Fuselage weight
                           W_{ht} : 5196  [N] Horizontal tail weight
                           W_{lg} : 1.212e+04  [N] Weight of landing gear
                          W_{pay} : 1.616e+05  [N] Payload weight
                           W_{vt} : 1963  [N] Vertical tail weight
                         W_{wing} : 1.204e+05  [N] Wing weight
                           W_{zf} : 4.729e+05  [N] Zero fuel weight
                        \bar{c}_w : 3.728  [m] Mean aerodynamic chord (wing)
                      \frac{L}{D} : 19.29 Lift/drag ratio
                              \xi : 0.07295 Takeoff parameter
                           b_{vt} : 3.707  [m] Vertical tail span
                           c_{vt} : 9.146  [m] Vertical tail root chord
                         h_{hold} : 0.7714  [m] Hold height
                         l_{fuse} : 52.47  [m] Fuselage length
                    p_{\lambda_v} : 1.6 1 + 2*Tail taper ratio
                         w_{fuse} : 3.71  [m] Fuselage width
                              x_w : 19.6  [m] Position of wing aerodynamic center
                     x_{CG_{eng}} : 19.6  [m] x-location of engine CG
                      x_{CG_{fu}} : 17.67  [m] x-location of fuselage CG
                      x_{CG_{ht}} : 51.32  [m] Horizontal tail CG location
                      x_{CG_{lg}} : 19.3  [m] x-location of landing gear CG
                      x_{CG_{vt}} : 47.9  [m] x-location of tail CG
                    x_{CG_{wing}} : 19.6  [m] x-location of wing CG
                           x_{CG} : 17.6  [m] x-location of CG
                           x_{TO} : 1524  [m] Takeoff distance
                           x_{up} : 29.61  [m] Fuselage upsweep point
                                y : 0.1318 Takeoff parameter
                          z_{bre} : 0.1036 Breguet parameter
                                    
               Fuselage, Aircraft | 
                        A_{floor} : 0.05549  [m**2] Floor beam x-sectional area
                         A_{fuse} : 10.81  [m**2] Fuselage x-sectional area
                         A_{hold} : 1.68  [m**2] Cargo hold x-sectional area
                         A_{skin} : 0.01087  [m**2] Skin cross sectional area
                     D_{friction} : 9646  [N] Friction drag
                      D_{upsweep} : 736.7  [N] Drag due to fuse upsweep
                               FF : 1.057 Fuselage form factor
                        M_{floor} : 4.442e+05  [N*m] Max bending moment in floor beams
                        P_{floor} : 1.137e+06  [N] Distributed floor load
                         R_{fuse} : 1.855  [m] Fuselage radius
                         S_{bulk} : 21.62  [m**2] Bulkhead surface area
                        S_{floor} : 5.686e+05  [N] Maximum shear in floor beams
                         S_{nose} : 49.82  [m**2] Nose surface area
                         V_{bulk} : 0.02016  [m**3] Bulkhead skin volume
                        V_{cabin} : 315  [m**3] Cabin volume
                        V_{cargo} : 6.796  [m**3] Cargo volume
                         V_{cone} : 0.1657  [m**3] Cone skin volume
                          V_{cyl} : 0.2653  [m**3] Cylinder skin volume
                        V_{floor} : 0.1734  [m**3] Floor volume
                         V_{hold} : 41.02  [m**3] Hold volume
                         V_{lugg} : 18.24  [m**3] Luggage volume
                         V_{nose} : 0.04645  [m**3] Nose skin volume
                          W_{apu} : 5657  [N] APU weight
                         W_{buoy} : 1517  [N] Buoyancy weight
                         W_{cone} : 7898  [N] Cone weight
                        W_{floor} : 9865  [N] Floor weight
                        W_{insul} : 4307  [N] Insulation material weight
                         W_{lugg} : 1.79e+04  [N] Passenger luggage weight
                         W_{padd} : 6.465e+04  [N] Misc weights (galley, toilets, doors etc.)
                         W_{pass} : 1.337e+05  [N] Passenger weight
                         W_{seat} : 2.79e+04  [N] Seating weight
                        W_{shell} : 1.582e+04  [N] Shell weight
                         W_{skin} : 8791  [N] Skin weight
                       W_{window} : 1.062e+04  [N] Window weight
                   \lambda_{cone} : 0.4 Tailcone radius taper ratio (xshell2->xtail)
                             \phi : 0.07883 Upsweep angle
                     \rho_{cabin} : 0.8711  [kg/m**3] Air density in cabin
                         \sigma_x : 3.831e+07  [N/m**2] Axial stress in skin
                  \sigma_{\theta} : 1.034e+08  [N/m**2] Skin hoop stress
                      \tau_{cone} : 1.034e+08  [N/m**2] Shear stress in cone
                                f : 14.14 Fineness ratio
                        h_{floor} : 0.08359  [m] Floor beam height
                         l_{cone} : 22.86  [m] Cone length
                        l_{floor} : 28.12  [m] Floor length
                         l_{nose} : 5.2  [m] Nose length
                        l_{shell} : 24.41  [m] Shell length
                         n_{pass} : 167 Number of passengers
                         n_{rows} : 31 Number of rows
                        t_{shell} : 0.001259  [m] Shell thickness
                         t_{skin} : 0.0009324  [m] Skin thickness
                        w_{floor} : 3.125  [m] Floor width
                           xVbulk : 0.5969  [m**4] Volume moment of bulkhead
                            xVcyl : 4.617  [m**4] Volume moment of cylinder
                           xVnose : 0.1208  [m**4] Volume moment of nose
                            xWapu : 2.069e+05  [N*m] Moment of APU
                           xWcone : 3.242e+05  [N*m] Moment of cone
                            xWfix : 2.802e+04  [N*m] Moment of fixed weights
                          xWfloor : 1.717e+05  [N*m] Moment of floor weight
                           xWfuse : 2.856e+06  [N*m] Fuselage moment
                          xWinsul : 7.496e+04  [N*m] Moment of insulation material
                           xWpadd : 1.125e+06  [N*m] Moment of misc weights
                           xWseat : 4.856e+05  [N*m] Moment of seats
                          xWshell : 2.543e+05  [N*m] Mass moment of shell
                           xWskin : 1.413e+05  [N*m] Mass moment of skin
                         xWwindow : 1.848e+05  [N*m] Mass moment of windows
                       x_{shell1} : 5.2  [m] Start of cylinder section
                       x_{shell2} : 29.61  [m] End of cylinder section
                                    
         HorizontalTail, Aircraft | 
                             AR_h : 7.707 Horizontal tail aspect ratio
                          C_{D_h} : 0.006697 Horizontal tail drag coefficient
                      C_{D_{0_h}} : 0.005348 Horizontal tail parasitic drag coefficient
                          C_{L_h} : 0.1785 Lift coefficient (htail)
                     C_{L_{ah_0}} : 5.348 Isolated lift curve slope (htail)
                       C_{L_{ah}} : 2.963 Lift curve slope (htail)
                              K_f : 0.4434 Empirical factor for fuselage-wing interference
                              L_h : 2.695e+04  [N] Horizontal tail downforce
                      L_{{max}_h} : 4.702e+05  [N] Maximum load
                         Re_{c_h} : 1.001e+07 Cruise Reynolds number (Horizontal tail)
                             S.M. : 0.05 Stability margin
                              S_h : 14.81  [m**2] Horizontal tail area
              \Delta x_{{lead}_h} : 32.57  [m] Distance from CG to horizontal tail leading edge
             \Delta x_{{trail}_h} : 34.88  [m] Distance from CG to horizontal tail trailing edge
                           \alpha : 0.06026 Horizontal tail angle of attack
                     \bar{c}_{ht} : 1.591  [m] Mean aerodynamic chord (ht)
                        \lambda_h : 0.2 Horizontal tail taper ratio
                           \tau_h : 0.15 Horizontal tail thickness/chord ratio
                           b_{ht} : 10.68  [m] Horizontal tail span
                       c_{root_h} : 2.31  [m] Horizontal tail root chord
                        c_{tip_h} : 0.462  [m] Horizontal tail tip chord
                              e_h : 0.9752 Oswald efficiency factor
                     f(\lambda_h) : 0.0033 Empirical efficiency function of taper
                           l_{ht} : 34.73  [m] Horizontal tail moment arm
                           p_{ht} : 1.4 Substituted variable = 1 + 2*taper
                           q_{ht} : 1.2 Substituted variable = 1 + taper
                 y_{\bar{c}_{ht}} : 3.052  [m] Spanwise location of mean aerodynamic chord
                                    
            LandingGear, Aircraft | 
                                B : 15.75  [m] Landing gear base
                         E_{land} : 2.657e+05  [J] Max KE to be absorbed in landing
                          F_{w_m} : 4458 Weight factor (main)
                          F_{w_n} : 400.9 Weight factor (nose)
                              I_m : 5.007e-06  [m**4] Area moment of inertia (main strut)
                              I_n : 7.204e-07  [m**4] Area moment of inertia (nose strut)
                              L_m : 4.489e+05  [N] Max static load through main gear
                              L_n : 1.122e+05  [N] Min static load through nose gear
                      L_{n_{dyn}} : 4.834e+04  [N] Dyn. braking load, nose gear
                          L_{w_m} : 1.122e+05  [N] Static load per wheel (main)
                          L_{w_n} : 5.612e+04  [N] Static load per wheel (nose)
                           S_{sa} : 0.2959  [m] Stroke of the shock absorber
                                T : 5.662  [m] Main landing gear track
                           W_{mg} : 1.1e+04  [N] Weight of main gear
                           W_{ms} : 1048  [N] Weight of main struts
                           W_{mw} : 1781  [N] Weight of main wheels (per strut)
                           W_{ng} : 1113  [N] Weight of nose gear
                           W_{ns} : 85.63  [N] Weight of nose strut
                           W_{nw} : 410.8  [N] Weight of nose wheels (total)
                         W_{wa,m} : 200.2  [lbf] Wheel assembly weight for single main gear wheel
                         W_{wa,n} : 46.18  [lbf] Wheel assembly weight for single nose gear wheel
                       \Delta x_m : 3.149  [m] Distance b/w main gear and CG
                       \Delta x_n : 12.6  [m] Distance b/w nose gear and CG
                       \tan(\phi) : 0.2679 Angle b/w main gear and CG
                       \tan(\psi) : 1.963 Tip over angle
                      d_{nacelle} : 2.05  [m] Nacelle diameter
                         d_{oleo} : 0.3119  [m] Diameter of oleo shock absorber
                          d_{t_m} : 39.72  [in] Diameter of main gear tires
                          d_{t_n} : 31.78  [in] Diameter of nose gear tires
                              l_m : 2.375  [m] Length of main gear
                              l_n : 1.627  [m] Length of nose gear
                         l_{oleo} : 0.7399  [m] Length of oleo shock absorber
                              r_m : 0.04181  [m] Radius of main gear struts
                              r_n : 0.04592  [m] Radius of nose gear struts
                              t_m : 0.02181  [m] Thickness of main gear strut wall
                              t_n : 0.002368  [m] Thickness of nose gear strut wall
                          w_{t_m} : 0.3439  [m] Width of main tires
                          w_{t_n} : 0.2751  [m] Width of nose tires
                              x_m : 20.75  [m] x-location of main gear
                              x_n : 5  [m] x-location of nose gear
                              y_m : 2.831  [m] y-location of main gear (symmetric)
                                    
           VerticalTail, Aircraft | 
                          A_{fan} : 2.405  [m**2] Engine reference area
                           A_{vt} : 0.7093 Vertical tail aspect ratio
                      C_{D_{vis}} : 0.004781 Viscous drag coefficient
                       C_{L_{vt}} : 0.3905 Vertical tail lift coefficient
                           D_{wm} : 3609  [N] Engine out windmill drag
                     L_{max_{vt}} : 1.28e+06  [N] Maximum load for structural sizing
                           L_{vt} : 2.271e+04  [N] Vertical tail lift in engine out
                          Re_{vt} : 4.099e+07 Vertical tail reynolds number, cruise
                                S : 38.76  [m**2] Vertical tail reference area (full)
                           S_{vt} : 19.38  [m**2] Vertical tail reference area (half)
                       W_{struct} : 3927  [N] Full span weight
                \Delta x_{lead_v} : 25.73  [m] Distance from CG to vertical tail leading edge
               \Delta x_{trail_v} : 34.88  [m] Distance from CG to vertical tail trailing edge
                     \bar{c}_{vt} : 6.519  [m] Vertical tail mean aero chord
                     \lambda_{vt} : 0.3 Vertical tail taper ratio
                        \tau_{vt} : 0.09958 Vertical tail thickness/chord ratio
                                b : 7.415  [m] Vertical tail full span
                    c_{root_{vt}} : 9.146  [m] Vertical tail root chord
                     c_{tip_{vt}} : 2.744  [m] Vertical tail tip chord
                           l_{vt} : 28.2  [m] Vertical tail moment arm
                           p_{vt} : 1.6 Substituted variable = 1 + 2*taper
                           q_{vt} : 1.3 Substituted variable = 1 + taper
                 z_{\bar{c}_{vt}} : 1.004  [m] Vertical location of mean aerodynamic chord
                                    
                   Wing, Aircraft | 
                          C_{D_w} : 0.01677 Drag coefficient, wing
                      C_{D_{p_w}} : 0.005637 Wing parasitic drag coefficient
                      L_{max_{w}} : 3.112e+06  [N] Maximum load
                             Re_w : 2.344e+07 Cruise Reynolds number (wing)
                    V_{fuel, max} : 93.38  [m**3] Available fuel volume
                         \alpha_w : 0.1 Wing angle of attack
              \bar{A}_{fuel, max} : 0.069 Non-dim. fuel area
                          \lambda : 0.2 Wing taper ratio
                           \tau_w : 0.15 Wing thickness/chord ratio
                              b_w : 30.19  [m] Wing span
                         c_{root} : 5.411  [m] Wing root chord
                          c_{tip} : 1.082  [m] Wing tip chord
                                e : 0.9702 Oswald efficiency factor
                     f(\lambda_w) : 0.0033 Empirical efficiency function of taper
                              p_w : 1.4 Substituted variable = 1 + 2*taper
                              q_w : 1.2 Substituted variable = 1 + taper
                    y_{\bar{c}_w} : 8.624  [m] Spanwise location of mean aerodynamic chord
                                    
WingBox, HorizontalTail, Aircraft | 
                          I_{cap} : 2.377e-05 Non-dim spar cap area moment of inertia
                              M_r : 2.114e+05  [N] Root moment per root chord
                          W_{cap} : 3382  [N] Weight of spar caps
                          W_{web} : 329.3  [N] Weight of shear web
                              \nu : 0.8612 Dummy variable = $(t^2 + t + 1)/(t+1)$
                          t_{cap} : 0.005417 Non-dim. spar cap thickness
                          t_{web} : 0.002345 Non-dim. shear web thickness
                                    
  WingBox, VerticalTail, Aircraft | 
                                A : 1.419 Aspect ratio
                          I_{cap} : 7.456e-07 Non-dim spar cap area moment of inertia
                              M_r : 1.21e+05  [N] Root moment per root chord
                          W_{cap} : 2107  [N] Weight of spar caps
                          W_{web} : 697.4  [N] Weight of shear web
                              \nu : 0.8225 Dummy variable = $(t^2 + t + 1)/(t+1)$
                          t_{cap} : 0.0003581 Non-dim. spar cap thickness
                          t_{web} : 0.0007934 Non-dim. shear web thickness
                                    
          WingBox, Wing, Aircraft | 
                          I_{cap} : 3.458e-05 Non-dim spar cap area moment of inertia
                              M_r : 1.688e+06  [N] Root moment per root chord
                          W_{cap} : 7.984e+04  [N] Weight of spar caps
                          W_{web} : 6158  [N] Weight of shear web
                              \nu : 0.8612 Dummy variable = $(t^2 + t + 1)/(t+1)$
                          t_{cap} : 0.008251 Non-dim. spar cap thickness
                          t_{web} : 0.002828 Non-dim. shear web thickness

Constants
---------
                Aircraft | 
             C_{L_{max}} : 2.5 Max lift coefficient
                       M : 0.78 Cruise Mach number
                N_{lift} : 2 Wing loading multiplier
                   Range : 3000  [nmi] Range
                     T_e : 1.29e+05  [N] Engine thrust at takeoff
                  V_{ne} : 144  [m/s] Never exceed velocity
                 W_{eng} : 1e+04  [N] Engine weight
                     \mu : 1.4e-05  [N*s/m**2] Dynamic viscosity (35,000 ft)
                    \rho : 0.38  [kg/m**3] Air density (35,000 ft)
                  \rho_0 : 1.225  [kg/m**3] Air density (0 ft)
              \rho_{cap} : 2700  [kg/m**3] Density of spar cap material
              \rho_{web} : 2700  [kg/m**3] Density of shear web material
      \sigma_{max,shear} : 1.67e+08  [Pa] Allowable shear stress
            \sigma_{max} : 2.5e+08  [Pa] Allowable tensile stress
                       a : 297  [m/s] Speed of sound (35,000 ft)
                     c_T : 0.3  [lb/hr/lbf] Thrust specific fuel consumption
                 d_{fan} : 1.75  [m] Fan diameter
               f_{w,add} : 0.4 Wing added weight fraction
                       g : 9.81  [m/s**2] Gravitational acceleration
                     l_r : 5000  [ft] Runway length
                     r_h : 0.75 Fractional wing thickness at spar web
                       w : 0.5 Wingbox-width-to-chord ratio
                 y_{eng} : 4.83  [m] Engine moment arm
                           
      Fuselage, Aircraft | 
                      LF : 0.898 Load factor
                N_{land} : 6 Emergency landing load factor
                       R : 287  [J/K/kg] Universal gas constant
                     SPR : 6 Number of seats per row
               T_{cabin} : 300  [K] Cabin temperature
             W''_{floor} : 60  [N/m**2] Floor weight/area density
             W''_{insul} : 22  [N/m**2] Weight/area density of insulation material
               W'_{seat} : 150  [N] Weight per seat
             W'_{window} : 435  [N/m] Weight/length density of windows
           W_{avg. pass} : 180  [lbf] Average passenger weight
               W_{cargo} : 1e+04  [N] Cargo weight
            W_{carry on} : 15  [lbf] Ave. carry-on weight
             W_{checked} : 40  [lbf] Ave. checked bag weight
                 W_{fix} : 3000  [lbf] Fixed weights (pilots, cockpit seats, navcom)
                \Delta h : 1  [m] Distance from floor to widest part of fuselage
                \Delta p : 5.2e+04  [Pa] Pressure difference across fuselage skin
           \rho_{\infty} : 0.38  [kg/m**3] Air density (35,000ft)
             \rho_{bend} : 2700  [kg/m**3] Stringer density
            \rho_{cargo} : 150  [kg/m**3] Cargo density
             \rho_{cone} : 2700  [kg/m**3] Cone material density
            \rho_{floor} : 2700  [kg/m**3] Floor material density
             \rho_{lugg} : 100  [kg/m**3] Luggage density
             \rho_{skin} : 2700  [kg/m**3] Skin density
          \sigma_{floor} : 2.069e+08  [N/m**2] Max allowable cap stress
           \sigma_{skin} : 1.034e+08  [N/m**2] Max allowable skin stress
            \tau_{floor} : 2.069e+08  [N/m**2] Max allowable shear web stress
                 f_{apu} : 0.035 APU weight as fraction of payload weight
                f_{fadd} : 0.2 Fractional added weight of local reinforcements
               f_{frame} : 0.25 Fractional frame weight
              f_{lugg,1} : 0.4 Proportion of passengers with one suitcase
              f_{lugg,2} : 0.1 Proportion of passengers with two suitcases
                f_{padd} : 0.4 Other misc weight as fraction of payload weight
              f_{string} : 0.35 Fractional weight of stringers
                n_{seat} : 186  Number of seats
                     p_s : 31  [in] Seat pitch
               p_{cabin} : 7.5e+04  [Pa] Cabin air pressure (8,000ft)
                     r_E : 1 Ratio of stringer/skin moduli
               w_{aisle} : 0.51  [m] Aisle width
                w_{seat} : 0.5  [m] Seat width
                 w_{sys} : 0.1  [m] Width between cabin and skin for systems
                    xapu : 120  [ft] x-location of APU
                    xfix : 2.1  [m] x-location of fixed weight
                           
HorizontalTail, Aircraft | 
            C_{L_{hmax}} : 2.5 Max lift coefficient
            C_{m_{fuse}} : 0.05 Moment coefficient (fuselage)
              S.M._{min} : 0.05 Minimum stability margin
              \Delta x_w : 2  [m] Distance from aerodynamic centre to CG
          \alpha_{max,h} : 0.1 Max angle of attack, htail
                  \eta_h : 0.97 Lift efficiency (diff between sectional and actual lift)
               \eta_{ht} : 0.9 Tail efficiency
      \tan(\Lambda_{ht}) : 0.5774 tangent of horizontal tail sweep
            |C_{m_{ac}}| : 0.1 Moment coefficient about aerodynamic centre (wing)
                           
   LandingGear, Aircraft | 
                       E : 205  [GPa] Modulus of elasticity, 4340 steel
                       K : 2 Column effective length factor
                     N_s : 2 Factor of safety
                  \eta_s : 0.8 Shock absorber efficiency
            \lambda_{LG} : 2.5 Ratio of max to static load
               \rho_{st} : 7850  [kg/m**3] Density of 4340 Steel
            \sigma_{y_c} : 4.7e+08  [Pa] Compressive yield strength 4340 steel
            \tan(\gamma) : 0.08749 Dihedral angle
        \tan(\phi_{min}) : 0.2679 Lower bound on phi
        \tan(\psi_{max}) : 1.963 Upper bound on psi
      \tan(\theta_{max}) : 0.2679 Max rotation angle
               f_{add,m} : 1.5 Proportional added weight, main
               f_{add,n} : 1.5 Proportional added weight, nose
             h_{nacelle} : 0.5  [m] Min. nacelle clearance
                  n_{mg} : 2 Number of main gear struts
                 n_{wps} : 2 Number of wheels per strut
                p_{oleo} : 1800  [lbf/in**2] Oleo pressure
             t_{nacelle} : 0.15  [m] Nacelle thickness
                 w_{ult} : 10  [ft/s] Ultimate velocity of descent
                  z_{CG} : 2  [m] CG height relative to bottom of fuselage
                z_{wing} : 0.5  [m] Height of wing relative to base of fuselage
                           
  VerticalTail, Aircraft | 
              C_{D_{wm}} : 0.5 Windmill drag coefficient
            C_{L_{vmax}} : 2.6 Max lift coefficient
                     V_1 : 70  [m/s] Minimum takeoff velocity
                  \rho_c : 0.38  [kg/m**3] Air density (35,000ft)
               \rho_{TO} : 1.225  [kg/m**3] Air density (SL))
      \tan(\Lambda_{vt}) : 0.8391 Tangent of leading edge sweep (40 deg)
              c_{l_{vt}} : 0.5 Sectional lift force coefficient (engine out)
                     e_v : 0.8 Span efficiency of vertical tail
                           
          Wing, Aircraft | 
            C_{L_{wmax}} : 2.5 Max lift coefficient, wing
          \alpha_{max,w} : 0.1 Max angle of attack
           \cos(\Lambda) : 0.866 Cosine of quarter-chord sweep angle
                  \eta_w : 0.97 Lift efficiency (diff b/w sectional, actual lift)
             \rho_{fuel} : 817  [kg/m**3] Density of fuel
           \tan(\Lambda) : 0.5774 Tangent of quarter-chord sweep angle

Sensitivities
-------------
          Wing, Aircraft | 
            C_{L_{wmax}} : 0.3192 Max lift coefficient, wing
           \tan(\Lambda) : 0.06735 Tangent of quarter-chord sweep angle
                  \eta_w : -0.1465 Lift efficiency (diff b/w sectional, actual lift)
          \alpha_{max,w} : -0.2513 Max angle of attack
                           
  VerticalTail, Aircraft | 
                  \rho_c : 0.03472 Air density (35,000ft)
            C_{L_{vmax}} : 0.0266 Max lift coefficient
                     e_v : -0.01605 Span efficiency of vertical tail
               \rho_{TO} : -0.04466 Air density (SL))
              c_{l_{vt}} : -0.05721 Sectional lift force coefficient (engine out)
                     V_1 : -0.1425 Minimum takeoff velocity
                           
   LandingGear, Aircraft | 
               f_{add,m} : 0.01242 Proportional added weight, main
      \tan(\theta_{max}) : 0.01065 Max rotation angle
                  n_{mg} : -0.02991 Number of main gear struts
                 n_{wps} : -0.0339 Number of wheels per strut
                           
HorizontalTail, Aircraft | 
              \Delta x_w : 0.05261 Distance from aerodynamic centre to CG
            C_{L_{hmax}} : 0.01468 Max lift coefficient
      \tan(\Lambda_{ht}) : 0.01088 tangent of horizontal tail sweep
                  \eta_h : -0.02993 Lift efficiency (diff between sectional and actual lift)
               \eta_{ht} : -0.04177 Tail efficiency
                           
      Fuselage, Aircraft | 
                n_{seat} : 0.7349  Number of seats
                      LF : 0.5094 Load factor
           W_{avg. pass} : 0.4494 Average passenger weight
                w_{seat} : 0.3223 Seat width
           \rho_{\infty} : 0.3064 Air density (35,000ft)
                     SPR : 0.1623 Number of seats per row
                     p_s : 0.16 Seat pitch
                f_{padd} : 0.1481 Other misc weight as fraction of payload weight
                \Delta h : 0.1033 Distance from floor to widest part of fuselage
               W'_{seat} : 0.06547 Weight per seat
             W_{checked} : 0.05992 Ave. checked bag weight
               w_{aisle} : 0.0548 Aisle width
              f_{lugg,1} : 0.03995 Proportion of passengers with one suitcase
             \rho_{skin} : 0.03587 Skin density
                \Delta p : 0.03587 Pressure difference across fuselage skin
               W_{cargo} : 0.03361 Cargo weight
                 W_{fix} : 0.02683 Fixed weights (pilots, cockpit seats, navcom)
             W'_{window} : 0.02433 Weight/length density of windows
             \rho_{cone} : 0.02152 Cone material density
                 w_{sys} : 0.02149 Width between cabin and skin for systems
              f_{lugg,2} : 0.01997 Proportion of passengers with two suitcases
                 f_{apu} : 0.01495 APU weight as fraction of payload weight
             W''_{floor} : 0.01208 Floor weight/area density
              f_{string} : 0.01116 Fractional weight of stringers
                N_{land} : 0.01052 Emergency landing load factor
            \rho_{floor} : 0.01052 Floor material density
           \sigma_{skin} : -0.05739 Max allowable skin stress
                           
                Aircraft | 
                       g : 1.435 Gravitational acceleration
                     c_T : 1.062 Thrust specific fuel consumption
                   Range : 1.062 Range
                  V_{ne} : 0.721 Never exceed velocity
                N_{lift} : 0.339 Wing loading multiplier
                  \rho_0 : 0.3339 Air density (0 ft)
              \rho_{cap} : 0.2748 Density of spar cap material
               f_{w,add} : 0.0849 Wing added weight fraction
                 y_{eng} : 0.07245 Engine moment arm
                     T_e : 0.07126 Engine thrust at takeoff
                     \mu : 0.04465 Dynamic viscosity (35,000 ft)
                 W_{eng} : 0.02314 Engine weight
                     r_h : 0.02241 Fractional wing thickness at spar web
              \rho_{web} : 0.02241 Density of shear web material
      \sigma_{max,shear} : -0.02241 Allowable shear stress
                       w : -0.04182 Wingbox-width-to-chord ratio
            \sigma_{max} : -0.3166 Allowable tensile stress
                    \rho : -0.4166 Air density (35,000 ft)
                       a : -1.164 Speed of sound (35,000 ft)
                       M : -1.312 Cruise Mach number

