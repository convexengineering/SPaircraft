Beginning signomial solve.
Solving took 5 GP solves and 154 seconds.

Cost
----
 1.7e+05 [N] 

Free Variables
--------------
                       Aircraft |                                                                       
                              D : 2.665e+04  [N]       Total aircraft drag (cruise)                     
                       D_{fuse} : 1.256e+04  [N]       Fuselage drag                                    
                         D_{ht} : 2520       [N]       Horizontal tail drag                             
                         D_{vt} : 1568       [N]       Vertical tail drag                               
                        I_{cap} : 4.555e-06            Non-dim spar cap area moment of inertia          
                    L_{v_{max}} : 9.648e+05  [N]       Maximum load for structural sizing               
                            M_r : 3.675e+05  [N]       Root moment per root chord                       
                              W : 2.867e+05  [N]       Total aircraft weight                            
                        W_{cap} : 7450       [N]       Weight of spar caps                              
                       W_{fuse} : 2.109e+05  [N]       Fuselage weight                                  
                         W_{ht} : 1.3e+04    [N]       Horizontal tail weight                           
                         W_{lg} : 1.627e+04  [N]       Weight of landing gear                           
                         W_{vt} : 6498       [N]       Vertical tail weight                             
                        W_{web} : 1832       [N]       Weight of shear web                              
                            \nu : 0.8225               Dummy variable = $(t^2 + t + 1)/(t+1)$           
                         b_{vt} : 6.46       [m]       Vertical tail span                               
                         c_{vt} : 6.958      [m]       Vertical tail root chord                         
                       h_{hold} : 0.9556     [m]       Hold height                                      
                       l_{fuse} : 47         [m]       Fuselage length                                  
                  p_{\lambda_v} : 1.6                  1 + 2*Tail taper ratio                           
                        t_{cap} : 0.0009704            Non-dim. spar cap thickness                      
                        t_{web} : 0.001061             Non-dim. shear web thickness                     
                       w_{fuse} : 4.148      [m]       Fuselage width                                   
                         x_{up} : 29.61      [m]       Fuselage upsweep point                           
                                                                                                        
             Fuselage, Aircraft |                                                                       
                      A_{floor} : 0.04625    [m**2]    Floor beam x-sectional area                      
                       A_{fuse} : 13.52      [m**2]    Fuselage x-sectional area                        
                       A_{hold} : 2.435      [m**2]    Cargo hold x-sectional area                      
                       A_{skin} : 0.01359    [m**2]    Skin cross sectional area                        
                   D_{friction} : 1.018e+04  [N]       Friction drag                                    
                    D_{upsweep} : 2377       [N]       Drag due to fuse upsweep                         
                             FF : 1.07                 Fuselage form factor                             
                      M_{floor} : 5.166e+05  [N*m]     Max bending moment in floor beams                
                      P_{floor} : 1.137e+06  [N]       Distributed floor load                           
                       R_{fuse} : 2.074      [m]       Fuselage radius                                  
                       S_{bulk} : 27.03      [m**2]    Bulkhead surface area                            
                      S_{floor} : 5.686e+05  [N]       Maximum shear in floor beams                     
                       S_{nose} : 56.3       [m**2]    Nose surface area                                
                       V_{bulk} : 0.02818    [m**3]    Bulkhead skin volume                             
                      V_{cabin} : 960.4      [m**3]    Cabin volume                                     
                      V_{cargo} : 6.796      [m**3]    Cargo volume                                     
                       V_{cone} : 0.2962     [m**3]    Cone skin volume                                 
                        V_{cyl} : 0.3317     [m**3]    Cylinder skin volume                             
                      V_{floor} : 0.1681     [m**3]    Floor volume                                     
                       V_{hold} : 59.45      [m**3]    Hold volume                                      
                       V_{lugg} : 18.24      [m**3]    Luggage volume                                   
                       V_{nose} : 0.0587     [m**3]    Nose skin volume                                 
                        W_{apu} : 5657       [N]       APU weight                                       
                       W_{buoy} : 4627       [N]       Buoyancy weight                                  
                       W_{cone} : 1.412e+04  [N]       Cone weight                                      
                      W_{floor} : 1.068e+04  [N]       Floor weight                                     
                      W_{insul} : 4857       [N]       Insulation material weight                       
                       W_{lugg} : 1.79e+04   [N]       Passenger luggage weight                         
                       W_{padd} : 6.465e+04  [N]       Misc weights (galley, toilets, doors etc.)       
                       W_{pass} : 1.337e+05  [N]       Passenger weight                                 
                        W_{pay} : 1.616e+05  [N]       Payload weight                                   
                       W_{seat} : 2.79e+04   [N]       Seating weight                                   
                      W_{shell} : 2.634e+04  [N]       Shell weight                                     
                       W_{skin} : 1.464e+04  [N]       Skin weight                                      
                     W_{window} : 1.062e+04  [N]       Window weight                                    
                 \lambda_{cone} : 0.4                  Tailcone radius taper ratio (xshell2->xtail)     
                           \phi : 0.1143               Upsweep angle                                    
                   \rho_{cabin} : 0.8711     [kg/m**3] Air density in cabin                             
                       \sigma_x : 3.831e+07  [N/m**2]  Axial stress in skin                             
                \sigma_{\theta} : 1.034e+08  [N/m**2]  Skin hoop stress                                 
                    \tau_{cone} : 1.034e+08  [N/m**2]  Shear stress in cone                             
                              f : 11.33                Fineness ratio                                   
                      h_{floor} : 0.1185     [m]       Floor I-beam height                              
                       l_{cone} : 17.39      [m]       Cone length                                      
                      l_{floor} : 28.56      [m]       Floor length                                     
                       l_{nose} : 5.2        [m]       Nose length                                      
                      l_{shell} : 24.41      [m]       Shell length                                     
                       n_{pass} : 167                  Number of passengers                             
                       n_{rows} : 31                   Number of rows                                   
                      t_{shell} : 0.001408   [m]       Shell thickness                                  
                       t_{skin} : 0.001043   [m]       Skin thickness                                   
                      w_{floor} : 3.634      [m]       Floor width                                      
                         xVbulk : 0.8345     [m**4]    Volume moment of bulkhead                        
                          xVcyl : 5.773      [m**4]    Volume moment of cylinder                        
                         xVnose : 0.1526     [m**4]    Volume moment of nose                            
                          xWapu : 2.069e+05  [N*m]     Moment of APU                                    
                         xWcone : 5.409e+05  [N*m]     Moment of cone                                   
                          xWfix : 2.802e+04  [N*m]     Moment of fixed weights                          
                        xWfloor : 1.859e+05  [N*m]     Moment of floor weight                           
                         xWfuse : 3.164e+06  [N*m]     Fuselage moment                                  
                        xWinsul : 8.454e+04  [N*m]     Moment of insulation material                    
                         xWpadd : 1.125e+06  [N*m]     Moment of misc weights                           
                         xWseat : 4.856e+05  [N*m]     Moment of seats                                  
                        xWshell : 3.223e+05  [N*m]     Mass moment of shell                             
                         xWskin : 1.79e+05   [N*m]     Mass moment of skin                              
                       xWwindow : 1.848e+05  [N*m]     Mass moment of windows                           
                     x_{shell1} : 5.2        [m]       Start of cylinder section                        
                     x_{shell2} : 29.61      [m]       End of cylinder section                          
                                                                                                        
       HorizontalTail, Aircraft |                                                                       
                           AR_h : 1.701                Horizontal tail aspect ratio                     
                        C_{D_h} : 0.006961             Horizontal tail drag coefficient                 
                    C_{D_{0_h}} : 0.005157             Horizontal tail parasitic drag coefficient       
                        C_{L_h} : 0.098                Lift coefficient (htail)                         
                     C_{L_{ah}} : 2.206                Lift curve slope (htail)                         
                            K_f : 0.424                Empirical factor for fuselage-wing interference  
                    L_{{max}_h} : 1.149e+06  [N]       Maximum load                                     
                       Re_{c_h} : 3.15e+07             Cruise Reynolds number (Horizontal tail)         
                           S.M. : 0.05                 Stability margin                                 
                            S_h : 34.8       [m**2]    Horizontal tail area                             
            \Delta x_{{lead}_h} : 25.05      [m]       Distance from CG to horizontal tail leading edge 
           \Delta x_{{trail}_h} : 32         [m]       Distance from CG to horizontal tail trailing edge
                         \alpha : 0.04443              Horizontal tail angle of attack                  
                   \bar{c}_{ht} : 4.96       [m]       Mean aerodynamic chord (ht)                      
                      \lambda_h : 0.3                  Horizontal tail taper ratio                      
                         \tau_h : 0.15                 Horizontal tail thickness/chord ratio            
                         b_{ht} : 7.695      [m]       Horizontal tail span                             
                   c_{{root}_h} : 6.958      [m]       Horizontal tail root chord                       
                    c_{{tip}_h} : 2.087      [m]       Horizontal tail tip chord                        
                            e_h : 0.9962               Oswald efficiency factor                         
                   f(\lambda_h) : 0.002259             Empirical efficiency function of taper           
                            l_h : 27.49      [m]       Horizontal tail moment arm                       
                         p_{ht} : 1.6                  Substituted variable = 1 + 2*taper               
                         q_{ht} : 1.3                  Substituted variable = 1 + taper                 
                            x_w : 17         [m]       Position of wing aerodynamic center              
                    y_{\bar{c}} : 2.084      [m]       Vertical location of mean aerodynamic chord      
                                                                                                        
          LandingGear, Aircraft |                                                                       
                              B : 12.5       [m]       Landing gear base                                
                       E_{land} : 3.809e+05  [J]       Max KE to be absorbed in landing                 
                        F_{w_m} : 7159                 Weight factor (main)                             
                        F_{w_n} : 643.7                Weight factor (nose)                             
                            I_m : 1.339e-05  [m**4]    Area moment of inertia (main strut)              
                            I_n : 1.715e-06  [m**4]    Area moment of inertia (nose strut)              
                            L_m : 6.435e+05  [N]       Max static load through main gear                
                            L_n : 1.609e+05  [N]       Min static load through nose gear                
                    L_{n_{dyn}} : 1.096e+05  [N]       Dyn. braking load, nose gear                     
                        L_{w_m} : 1.609e+05  [N]       Static load per wheel (main)                     
                        L_{w_n} : 8.044e+04  [N]       Static load per wheel (nose)                     
                           S_sa : 0.2959     [m]       Stroke of the shock absorber                     
                              T : 9.66       [m]       Main landing gear track                          
                         W_{mg} : 1.469e+04  [N]       Weight of main gear                              
                         W_{ms} : 1405       [N]       Weight of main struts                            
                         W_{mw} : 2377       [N]       Weight of main wheels (per strut)                
                         W_{ng} : 1580       [N]       Weight of nose gear                              
                         W_{ns} : 209.7      [N]       Weight of nose strut                             
                         W_{nw} : 548.1      [N]       Weight of nose wheels (total)                    
                       W_{wa,m} : 267.2      [lbf]     Wheel assembly weight for single main gear wheel 
                       W_{wa,n} : 61.61      [lbf]     Wheel assembly weight for single nose gear wheel 
                     \Delta x_m : 2.5        [m]       Distance b/w main gear and CG                    
                     \Delta x_n : 10         [m]       Distance b/w nose gear and CG                    
                     \tan(\phi) : 0.2679               Angle b/w main gear and CG                       
                     \tan(\psi) : 1.963                Tip over angles                                  
                    d_{nacelle} : 2.05       [m]       Nacelle diameter                                 
                       d_{oleo} : 0.3735     [m]       Diameter of oleo shock absorber                  
                        d_{t_m} : 44.5       [in]      Diameter of main gear tires                      
                        d_{t_n} : 35.6       [in]      Diameter of nose gear tires                      
                            l_m : 3.245      [m]       Length of main gear                              
                            l_n : 2.322      [m]       Length of nose gear                              
                       l_{oleo} : 0.7399     [m]       Length of oleo shock absorber                    
                            r_m : 0.06902    [m]       Radius of main gear struts                       
                            r_n : 0.05408    [m]       Radius of nose gear struts                       
                            t_m : 0.01297    [m]       Thickness of main gear strut wall                
                            t_n : 0.003451   [m]       Thickness of nose gear strut wall                
                        w_{t_m} : 0.4088     [m]       Width of main tires                              
                        w_{t_n} : 0.327      [m]       Width of nose tires                              
                            x_m : 17.5       [m]       x-location of main gear                          
                            x_n : 5          [m]       x-location of nose gear                          
                            y_m : 4.83       [m]       y-location of main gear (symmetric)              
                                                                                                        
         VerticalTail, Aircraft |                                                                       
                        A_{fan} : 2.405      [m**2]    Engine reference area                            
                         A_{vt} : 1.428                Vertical tail aspect ratio                       
                    C_{D_{vis}} : 0.005157             Viscous drag coefficient                         
                     C_{L_{vt}} : 0.4389               Vertical tail lift coefficient                   
                         D_{wm} : 3112       [N]       Engine out windmill drag                         
                        L_{max} : 1.93e+06   [N]       Maximum wing load                                
                         L_{vt} : 3.318e+04  [N]       Vertical tail lift in engine out                 
                           Re_c : 3.15e+07             Vertical tail reynolds number, cruise            
                              S : 58.44      [m**2]    Reference area                                   
                         S_{vt} : 29.22      [m**2]    Vertical tail ref. area (half)                   
                     W_{struct} : 1.3e+04    [N]       Full span weight                                 
                \Delta x_{lead} : 16.52      [m]       Distance from CG to vertical tail leading edge   
               \Delta x_{trail} : 23.48      [m]       Distance from CG to vertical tail trailing edge  
                        \bar{c} : 4.96       [m]       Vertical tail mean aero chord                    
                   \lambda_{vt} : 0.3                  Vertical tail taper ratio                        
                           \tau : 0.15                 Thickness to chord ratio                         
                              b : 12.92      [m]       Span                                             
                       c_{root} : 6.958      [m]       Vertical tail root chord                         
                        c_{tip} : 2.087      [m]       Vertical tail tip chord                          
                         l_{vt} : 19.23      [m]       Vertical tail moment arm                         
                              p : 1.6                  Substituted variable = 1 + 2*taper               
                              q : 1.3                  Substituted variable = 1 + taper                 
                    z_{\bar{c}} : 1.75       [m]       Vertical location of mean aerodynamic chord      
                                                                                                        
WingBox, VerticalTail, Aircraft |                                                                       
                              A : 2.857                Aspect ratio                                     

Constants
---------
                Aircraft |                                                                                
                D_{wing} : 1e+04      [N]         Wing drag                                               
                N_{lift} : 2                      Wing loading multiplier                                 
              V_{\infty} : 234        [m/s]       Freestream velocity                                     
                  V_{ne} : 144        [m/s]       Never exceed velocity                                   
                 W_{eng} : 1e+04      [N]         Engine weight                                           
                W_{wing} : 3e+04      [N]         Wing weight                                             
                     \mu : 1.4e-05    [N*s/m**2]  Dynamic viscosity (35,000ft)                            
              \rho_{cap} : 2700       [kg/m**3]   Density of spar cap material                            
              \rho_{web} : 2700       [kg/m**3]   Density of shear web material                           
      \sigma_{max,shear} : 1.67e+08   [Pa]        Allowable shear stress                                  
            \sigma_{max} : 2.5e+08    [Pa]        Allowable tensile stress                                
                 d_{fan} : 1.75       [m]         Fan diameter                                            
               f_{w,add} : 0.4                    Wing added weight fraction                              
                       g : 9.81       [m/s**2]    Gravitational acceleration                              
                     r_h : 0.75                   Fractional wing thickness at spar web                   
                       w : 0.5                    Wingbox-width-to-chord ratio                            
                  x_{CG} : 15         [m]         CG location                                             
                 y_{eng} : 4.83       [m]         Engine moment arm                                       
                                                                                                          
      Fuselage, Aircraft |                                                                                
                      LF : 0.898                  Load factor                                             
                N_{land} : 6                      Emergency landing load factor                           
                       R : 287        [J/K/kg]    Universal gas constant                                  
                     SPR : 6                      Number of seats per row                                 
               T_{cabin} : 300        [K]         Cabin temperature                                       
             W''_{floor} : 60         [N/m**2]    Floor weight/area density                               
             W''_{insul} : 22         [N/m**2]    Weight/area density of insulation material              
               W'_{seat} : 150        [N]         Weight per seat                                         
             W'_{window} : 435        [N/m]       Weight/length density of windows                        
           W_{avg. pass} : 180        [lbf]       Average passenger weight                                
               W_{cargo} : 1e+04      [N]         Cargo weight                                            
            W_{carry on} : 15         [lbf]       Ave. carry-on weight                                    
             W_{checked} : 40         [lbf]       Ave. checked bag weight                                 
                 W_{fix} : 3000       [lbf]       Fixed weights (pilots, cockpit seats, navcom)           
                \Delta h : 1          [m]         Distance from floor to widest part of fuselage          
                \Delta p : 5.2e+04    [Pa]        Pressure difference across fuselage skin                
           \rho_{\infty} : 0.38       [kg/m**3]   Air density (35,000ft)                                  
             \rho_{bend} : 2700       [kg/m**3]   Stringer density                                        
            \rho_{cargo} : 150        [kg/m**3]   Cargo density                                           
             \rho_{cone} : 2700       [kg/m**3]   Cone material density                                   
            \rho_{floor} : 2700       [kg/m**3]   Floor material density                                  
             \rho_{lugg} : 100        [kg/m**3]   Luggage density                                         
             \rho_{skin} : 2700       [kg/m**3]   Skin density                                            
          \sigma_{floor} : 2.069e+08  [N/m**2]    Max allowable cap stress                                
           \sigma_{skin} : 1.034e+08  [N/m**2]    Max allowable skin stress                               
            \tau_{floor} : 2.069e+08  [N/m**2]    Max allowable shear web stress                          
                 f_{apu} : 0.035                  APU weight as fraction of payload weight                
                f_{fadd} : 0.2                    Fractional added weight of local reinforcements         
               f_{frame} : 0.25                   Fractional frame weight                                 
              f_{lugg,1} : 0.4                    Proportion of passengers with one suitcase              
              f_{lugg,2} : 0.1                    Proportion of passengers with two suitcases             
                f_{padd} : 0.4                    Other misc weight as fraction of payload weight         
              f_{string} : 0.35                   Fractional weight of stringers                          
                n_{seat} : 186                     Number of seats                                        
                     p_s : 31         [in]        Seat pitch                                              
               p_{cabin} : 7.5e+04    [Pa]        Cabin air pressure (8,000ft)                            
                     r_E : 1                      Ratio of stringer/skin moduli                           
               w_{aisle} : 0.51       [m]         Aisle width                                             
                w_{seat} : 0.5        [m]         Seat width                                              
                 w_{sys} : 0.1        [m]         Width between cabin and skin for systems                
             x_{CG_{fu}} : 15         [m]         x-location of fuselage CG                               
                    xapu : 120        [ft]        x-location of APU                                       
                    xfix : 2.1        [m]         x-location of fixed weight                              
                                                                                                          
HorizontalTail, Aircraft |                                                                                
                 C_{L_w} : 0.5                    Lift coefficient (wing)                                 
              C_{L_{aw}} : 6.283                  Lift curve slope (wing)                                 
            C_{L_{hmax}} : 2.6                    Max lift coefficient                                    
            C_{m_{fuse}} : 0.05                   Moment coefficient (fuselage)                           
              S.M._{min} : 0.05                   Minimum stability margin                                
                     S_w : 125        [m**2]      Wing area                                               
              \Delta x_w : 2          [m]         Distance from aerodynamic centre to CG                  
            \alpha_{max} : 0.1                    Max angle of attack (htail)                             
          \bar{c}_{wing} : 5          [m]         Mean aerodynamic chord (wing)                           
                  \eta_h : 0.97                   Lift efficiency (diff between sectional and actual lift)
                    \rho : 0.38       [kg/m**3]   Air density (35,000 ft)                                 
                  \rho_0 : 1.225      [kg/m**3]   Air density (0 ft)                                      
         \tan(\Lambda_h) : 0.5774                 tangent of horizontal tail sweep                        
            |C_{m_{ac}}| : 0.1                    Moment coefficient about aerodynamic centre (wing)      
                                                                                                          
   LandingGear, Aircraft |                                                                                
                       E : 205        [GPa]       Modulus of elasticity, 4340 steel                       
                       K : 2                      Column effective length factor                          
                     N_s : 2                      Factor of safety                                        
                     W_0 : 8.044e+05  [N]         Weight of aircraft excluding landing gear               
                  \eta_s : 0.8                    Shock absorber efficiency                               
            \lambda_{LG} : 2.5                    Ratio of max to static load                             
               \rho_{st} : 7850       [kg/m**3]   Density of 4340 Steel                                   
            \sigma_{y_c} : 4.7e+08    [Pa]        Compressive yield strength 4340 steel                   
            \tan(\gamma) : 0.08749                Tangent, dihedral angle                                 
        \tan(\phi_{min}) : 0.2679                 Lower bound on phi                                      
        \tan(\psi_{max}) : 1.963                  Upper bound on psi                                      
       \tan(\theta_{TO}) : 0.2679                 Takeoff pitch angle                                     
               f_{add,m} : 1.5                    Proportional added weight, main                         
               f_{add,n} : 1.5                    Proportional added weight, nose                         
             h_{nacelle} : 0.5        [m]         Min. nacelle clearance                                  
                  n_{mg} : 2                      Number of main gear struts                              
                 n_{wps} : 2                      Number of wheels per strut                              
                p_{oleo} : 1800       [lbf/in**2] Oleo pressure                                           
             t_{nacelle} : 0.15       [m]         Nacelle thickness                                       
                 w_{ult} : 10         [ft/s]      Ultimate velocity of descent                            
                x_{CG_0} : 18         [m]         x-location of CG excl. LG                               
                  z_{CG} : 2          [m]         CG height relative to bottom of fuselage                
                z_{wing} : 0.5        [m]         Height of wing relative to base of fuselage             
                                                                                                          
  VerticalTail, Aircraft |                                                                                
              C_{D_{wm}} : 0.5                    Windmill drag coefficient                               
            C_{L_{vmax}} : 2.6                    Max lift coefficient                                    
                     T_e : 1.29e+05   [N]         Thrust per engine at takeoff                            
                     V_1 : 65         [m/s]       Minimum takeoff velocity                                
                     V_c : 234        [m/s]       Cruise velocity                                         
                  \rho_c : 0.38       [kg/m**3]   Air density (35,000ft)                                  
               \rho_{TO} : 1.225      [kg/m**3]   Air density (SL))                                       
      \tan(\Lambda_{LE}) : 0.8391                 Tangent of leading edge sweep (40 deg)                  
              c_{l_{vt}} : 0.5                    Sectional lift force coefficient (engine out)           
                       e : 0.8                    Span efficiency of vertical tail                        
             x_{CG_{vt}} : 35         [m]         x-location of tail CG                                   

Sensitivities
-------------
  VerticalTail, Aircraft |                                                         
                     T_e : 0.3123   Thrust per engine at takeoff                   
            C_{L_{vmax}} : 0.1641   Max lift coefficient                           
                     V_c : 0.01826  Cruise velocity                                
      \tan(\Lambda_{LE}) : -0.02442 Tangent of leading edge sweep (40 deg)         
                       e : -0.0391  Span efficiency of vertical tail               
               \rho_{TO} : -0.1482  Air density (SL))                              
              c_{l_{vt}} : -0.2807  Sectional lift force coefficient (engine out)  
             x_{CG_{vt}} : -0.3326  x-location of tail CG                          
                     V_1 : -0.6246  Minimum takeoff velocity                       
                                                                                   
   LandingGear, Aircraft |                                                         
                     W_0 : 0.08551  Weight of aircraft excluding landing gear      
       \tan(\theta_{TO}) : 0.02623  Takeoff pitch angle                            
               f_{add,m} : 0.02097  Proportional added weight, main                
                       K : 0.01694  Column effective length factor                 
                  n_{mg} : -0.03137 Number of main gear struts                     
                 n_{wps} : -0.03806 Number of wheels per strut                     
                                                                                   
HorizontalTail, Aircraft |                                                         
                    \rho : 0.0146   Air density (35,000 ft)                        
              \Delta x_w : 0.01156  Distance from aerodynamic centre to CG         
                 C_{L_w} : 0.01024  Lift coefficient (wing)                        
                                                                                   
      Fuselage, Aircraft |                                                         
                n_{seat} : 0.9386    Number of seats                               
                     p_s : 0.584    Seat pitch                                     
                      LF : 0.2571   Load factor                                    
           W_{avg. pass} : 0.2268   Average passenger weight                       
                f_{padd} : 0.2206   Other misc weight as fraction of payload weight
             \rho_{cone} : 0.1061   Cone material density                          
                \Delta h : 0.1016   Distance from floor to widest part of fuselage 
               W'_{seat} : 0.09745  Weight per seat                                
             \rho_{skin} : 0.06319  Skin density                                   
                \Delta p : 0.06319  Pressure difference across fuselage skin       
           \rho_{\infty} : 0.06189  Air density (35,000ft)                         
                 f_{apu} : 0.04057  APU weight as fraction of payload weight       
                    xapu : 0.04057  x-location of APU                              
             W'_{window} : 0.03624  Weight/length density of windows               
              f_{string} : 0.03291  Fractional weight of stringers                 
             W_{checked} : 0.03025  Ave. checked bag weight                        
               f_{frame} : 0.02351  Fractional frame weight                        
             W''_{floor} : 0.02125  Floor weight/area density                      
              f_{lugg,1} : 0.02016  Proportion of passengers with one suitcase     
                f_{fadd} : 0.01881  Fractional added weight of local reinforcements
               W_{cargo} : 0.01696  Cargo weight                                   
             W''_{insul} : 0.01658  Weight/area density of insulation material     
                N_{land} : 0.01519  Emergency landing load factor                  
            \rho_{floor} : 0.01519  Floor material density                         
              f_{lugg,2} : 0.01008  Proportion of passengers with two suitcases    
          \sigma_{floor} : -0.01384 Max allowable cap stress                       
           \sigma_{skin} : -0.1692  Max allowable skin stress                      
                     SPR : -0.584   Number of seats per row                        
             x_{CG_{fu}} : -0.6204  x-location of fuselage CG                      
                                                                                   
                Aircraft |                                                         
                  V_{ne} : 0.3281   Never exceed velocity                          
                 y_{eng} : 0.3196   Engine moment arm                              
                       g : 0.2507   Gravitational acceleration                     
              V_{\infty} : 0.1652   Freestream velocity                            
                  x_{CG} : 0.1022   CG location                                    
                W_{wing} : 0.08823  Wing weight                                    
                D_{wing} : 0.05882  Wing drag                                      
                N_{lift} : 0.058    Wing loading multiplier                        
              \rho_{cap} : 0.04601  Density of spar cap material                   
                 W_{eng} : 0.02941  Engine weight                                  
               f_{w,add} : 0.01638  Wing added weight fraction                     
                 d_{fan} : 0.01507  Fan diameter                                   
                     \mu : 0.01238  Dynamic viscosity (35,000ft)                   
                     r_h : 0.01132  Fractional wing thickness at spar web          
              \rho_{web} : 0.01132  Density of shear web material                  
      \sigma_{max,shear} : -0.01132 Allowable shear stress                         
            \sigma_{max} : -0.04668 Allowable tensile stress                       

