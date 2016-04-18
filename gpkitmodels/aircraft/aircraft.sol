Beginning signomial solve.
Solving took 5 GP solves and 173 seconds.

Cost
----
 2.036e+05 [N] 

Free Variables
--------------
                       Aircraft |                                                                       
                        C_{L_w} : 0.1686               Lift coefficient (wing)                          
                     C_{L_{aw}} : 1.686                Lift curve slope (wing)                          
                              D : 9.395e+04  [N]       Total aircraft drag (cruise)                     
                       D_{fuse} : 1.6e+04    [N]       Fuselage drag                                    
                         D_{ht} : 1898       [N]       Horizontal tail drag                             
                         D_{vt} : 1145       [N]       Vertical tail drag                               
                       D_{wing} : 7.491e+04  [N]       Wing drag                                        
                        I_{cap} : 5.718e-07            Non-dim spar cap area moment of inertia          
                    L_{v_{max}} : 7.133e+05  [N]       Maximum load for structural sizing               
                            M_r : 2.768e+05  [N]       Root moment per root chord                       
                              W : 2.193e+05  [N]       Total aircraft weight                            
                        W_{cap} : 4677       [N]       Weight of spar caps                              
                       W_{fuse} : 1.663e+05  [N]       Fuselage weight                                  
                         W_{ht} : 1.1e+04    [N]       Horizontal tail weight                           
                         W_{lg} : 1.553e+04  [N]       Weight of landing gear                           
                         W_{vt} : 5501       [N]       Vertical tail weight                             
                        W_{web} : 3181       [N]       Weight of shear web                              
                       W_{wing} : 1.1e+04    [N]       Wing weight                                      
                            \nu : 0.8612               Dummy variable = $(t^2 + t + 1)/(t+1)$           
                         b_{vt} : 1.95       [m]       Vertical tail span                               
                         c_{vt} : 17.04      [m]       Vertical tail root chord                         
                       h_{hold} : 0.9165     [m]       Hold height                                      
                       l_{fuse} : 85.06      [m]       Fuselage length                                  
                  p_{\lambda_v} : 1.6                  1 + 2*Tail taper ratio                           
                        t_{cap} : 0.0001203            Non-dim. spar cap thickness                      
                        t_{web} : 0.0003637            Non-dim. shear web thickness                     
                       w_{fuse} : 4.035      [m]       Fuselage width                                   
                    x_{CG_{fu}} : 18.4       [m]       x-location of fuselage CG                        
                    x_{CG_{ht}} : 76.54      [m]       Horizontal tail CG location                      
                    x_{CG_{lg}} : 24.3       [m]       x-location of landing gear CG                    
                    x_{CG_{vt}} : 76.54      [m]       x-location of tail CG                            
                         x_{CG} : 22.87      [m]       x-location of CG                                 
                         x_{up} : 29.61      [m]       Fuselage upsweep point                           
                                                                                                        
             Fuselage, Aircraft |                                                                       
                      A_{floor} : 0.05188    [m**2]    Floor beam x-sectional area                      
                       A_{fuse} : 12.79      [m**2]    Fuselage x-sectional area                        
                       A_{hold} : 2.251      [m**2]    Cargo hold x-sectional area                      
                       A_{skin} : 0.01285    [m**2]    Skin cross sectional area                        
                   D_{friction} : 1.576e+04  [N]       Friction drag                                    
                    D_{upsweep} : 242.9      [N]       Drag due to fuse upsweep                         
                             FF : 1.059                Fuselage form factor                             
                      M_{floor} : 4.981e+05  [N*m]     Max bending moment in floor beams                
                      P_{floor} : 1.137e+06  [N]       Distributed floor load                           
                       R_{fuse} : 2.017      [m]       Fuselage radius                                  
                       S_{bulk} : 25.57      [m**2]    Bulkhead surface area                            
                      S_{floor} : 5.686e+05  [N]       Maximum shear in floor beams                     
                       S_{nose} : 54.6       [m**2]    Nose surface area                                
                       V_{bulk} : 0.02593    [m**3]    Bulkhead skin volume                             
                      V_{cabin} : 373.9      [m**3]    Cabin volume                                     
                      V_{cargo} : 6.796      [m**3]    Cargo volume                                     
                       V_{cone} : 0.1664     [m**3]    Cone skin volume                                 
                        V_{cyl} : 0.3138     [m**3]    Cylinder skin volume                             
                      V_{floor} : 0.1818     [m**3]    Floor volume                                     
                       V_{hold} : 54.95      [m**3]    Hold volume                                      
                       V_{lugg} : 18.24      [m**3]    Luggage volume                                   
                       V_{nose} : 0.05537    [m**3]    Nose skin volume                                 
                        W_{apu} : 5657       [N]       APU weight                                       
                       W_{buoy} : 1801       [N]       Buoyancy weight                                  
                       W_{cone} : 7935       [N]       Cone weight                                      
                      W_{floor} : 1.08e+04   [N]       Floor weight                                     
                      W_{insul} : 4714       [N]       Insulation material weight                       
                       W_{lugg} : 1.79e+04   [N]       Passenger luggage weight                         
                       W_{padd} : 6.465e+04  [N]       Misc weights (galley, toilets, doors etc.)       
                       W_{pass} : 1.337e+05  [N]       Passenger weight                                 
                        W_{pay} : 1.616e+05  [N]       Payload weight                                   
                       W_{seat} : 2.79e+04   [N]       Seating weight                                   
                      W_{shell} : 1.883e+04  [N]       Shell weight                                     
                       W_{skin} : 1.046e+04  [N]       Skin weight                                      
                     W_{window} : 1.062e+04  [N]       Window weight                                    
                 \lambda_{cone} : 0.4                  Tailcone radius taper ratio (xshell2->xtail)     
                           \phi : 0.04692              Upsweep angle                                    
                   \rho_{cabin} : 0.8711     [kg/m**3] Air density in cabin                             
                       \sigma_x : 3.831e+07  [N/m**2]  Axial stress in skin                             
                \sigma_{\theta} : 1.034e+08  [N/m**2]  Skin hoop stress                                 
                    \tau_{cone} : 1.034e+08  [N/m**2]  Shear stress in cone                             
                              f : 21.08                Fineness ratio                                   
                      h_{floor} : 0.1008     [m]       Floor I-beam height                              
                       l_{cone} : 42.61      [m]       Cone length                                      
                      l_{floor} : 28.44      [m]       Floor length                                     
                       l_{nose} : 5.2        [m]       Nose length                                      
                      l_{shell} : 24.41      [m]       Shell length                                     
                       n_{pass} : 167                  Number of passengers                             
                       n_{rows} : 31                   Number of rows                                   
                      t_{shell} : 0.001369   [m]       Shell thickness                                  
                       t_{skin} : 0.001014   [m]       Skin thickness                                   
                      w_{floor} : 3.504      [m]       Floor width                                      
                         xVbulk : 0.7678     [m**4]    Volume moment of bulkhead                        
                          xVcyl : 5.461      [m**4]    Volume moment of cylinder                        
                         xVnose : 0.144      [m**4]    Volume moment of nose                            
                          xWapu : 2.069e+05  [N*m]     Moment of APU                                    
                         xWcone : 4.55e+05   [N*m]     Moment of cone                                   
                          xWfix : 2.802e+04  [N*m]     Moment of fixed weights                          
                        xWfloor : 1.879e+05  [N*m]     Moment of floor weight                           
                         xWfuse : 3.059e+06  [N*m]     Fuselage moment                                  
                        xWinsul : 8.204e+04  [N*m]     Moment of insulation material                    
                         xWpadd : 1.125e+06  [N*m]     Moment of misc weights                           
                         xWseat : 4.856e+05  [N*m]     Moment of seats                                  
                        xWshell : 3.038e+05  [N*m]     Mass moment of shell                             
                         xWskin : 1.688e+05  [N*m]     Mass moment of skin                              
                       xWwindow : 1.848e+05  [N*m]     Mass moment of windows                           
                     x_{shell1} : 5.2        [m]       Start of cylinder section                        
                     x_{shell2} : 29.61      [m]       End of cylinder section                          
                                                                                                        
       HorizontalTail, Aircraft |                                                                       
                           AR_h : 0.3413               Horizontal tail aspect ratio                     
                        C_{D_h} : 0.005113             Horizontal tail drag coefficient                 
                    C_{D_{0_h}} : 0.005096             Horizontal tail parasitic drag coefficient       
                        C_{L_h} : 0.004318             Lift coefficient (htail)                         
                     C_{L_{ah}} : 0.5306               Lift curve slope (htail)                         
                            K_f : 0.3187               Empirical factor for fuselage-wing interference  
                    L_{{max}_h} : 1.178e+06  [N]       Maximum load                                     
                       Re_{c_h} : 7.457e+07            Cruise Reynolds number (Horizontal tail)         
                           S.M. : 0.05                 Stability margin                                 
                            S_h : 35.68      [m**2]    Horizontal tail area                             
            \Delta x_{{lead}_h} : 45.15      [m]       Distance from CG to horizontal tail leading edge 
           \Delta x_{{trail}_h} : 62.19      [m]       Distance from CG to horizontal tail trailing edge
                         \alpha : 0.008137             Horizontal tail angle of attack                  
                   \bar{c}_{ht} : 11.74      [m]       Mean aerodynamic chord (ht)                      
                 \bar{c}_{wing} : 5.546      [m]       Mean aerodynamic chord (wing)                    
                      \lambda_h : 0.2                  Horizontal tail taper ratio                      
                         \tau_h : 0.15                 Horizontal tail thickness/chord ratio            
                         b_{ht} : 3.49       [m]       Horizontal tail span                             
                     c_{root_h} : 17.04      [m]       Horizontal tail root chord                       
                      c_{tip_h} : 3.408      [m]       Horizontal tail tip chord                        
                            e_h : 0.9989               Oswald efficiency factor                         
                   f(\lambda_h) : 0.0033               Empirical efficiency function of taper           
                            l_h : 48.66      [m]       Horizontal tail moment arm                       
                         p_{ht} : 1.4                  Substituted variable = 1 + 2*taper               
                         q_{ht} : 1.2                  Substituted variable = 1 + taper                 
                            x_w : 24.87      [m]       Position of wing aerodynamic center              
               y_{\bar{c}_{ht}} : 0.997      [m]       Vertical location of mean aerodynamic chord      
                                                                                                        
          LandingGear, Aircraft |                                                                       
                              B : 13.8       [m]       Landing gear base                                
                       E_{land} : 3.809e+05  [J]       Max KE to be absorbed in landing                 
                        F_{w_m} : 7159                 Weight factor (main)                             
                        F_{w_n} : 643.7                Weight factor (nose)                             
                            I_m : 7.185e-06  [m**4]    Area moment of inertia (main strut)              
                            I_n : 1.114e-06  [m**4]    Area moment of inertia (nose strut)              
                            L_m : 6.435e+05  [N]       Max static load through main gear                
                            L_n : 1.609e+05  [N]       Min static load through nose gear                
                    L_{n_{dyn}} : 7.906e+04  [N]       Dyn. braking load, nose gear                     
                        L_{w_m} : 1.609e+05  [N]       Static load per wheel (main)                     
                        L_{w_n} : 8.044e+04  [N]       Static load per wheel (nose)                     
                           S_sa : 0.2959     [m]       Stroke of the shock absorber                     
                              T : 5.692      [m]       Main landing gear track                          
                         W_{mg} : 1.403e+04  [N]       Weight of main gear                              
                         W_{ms} : 1074       [N]       Weight of main struts                            
                         W_{mw} : 2377       [N]       Weight of main wheels (per strut)                
                         W_{ng} : 1498       [N]       Weight of nose gear                              
                         W_{ns} : 128        [N]       Weight of nose strut                             
                         W_{nw} : 548.1      [N]       Weight of nose wheels (total)                    
                       W_{wa,m} : 267.2      [lbf]     Wheel assembly weight for single main gear wheel 
                       W_{wa,n} : 61.61      [lbf]     Wheel assembly weight for single nose gear wheel 
                     \Delta x_m : 2.761      [m]       Distance b/w main gear and CG                    
                     \Delta x_n : 11.04      [m]       Distance b/w nose gear and CG                    
                     \tan(\phi) : 0.2679               Angle b/w main gear and CG                       
                     \tan(\psi) : 1.963                Tip over angles                                  
                    d_{nacelle} : 2.05       [m]       Nacelle diameter                                 
                       d_{oleo} : 0.3735     [m]       Diameter of oleo shock absorber                  
                        d_{t_m} : 44.5       [in]      Diameter of main gear tires                      
                        d_{t_n} : 35.6       [in]      Diameter of nose gear tires                      
                            l_m : 2.376      [m]       Length of main gear                              
                            l_n : 1.627      [m]       Length of nose gear                              
                       l_{oleo} : 0.7399     [m]       Length of oleo shock absorber                    
                            r_m : 0.04948    [m]       Radius of main gear struts                       
                            r_n : 0.04672    [m]       Radius of nose gear struts                       
                            t_m : 0.01888    [m]       Thickness of main gear strut wall                
                            t_n : 0.003479   [m]       Thickness of nose gear strut wall                
                        w_{t_m} : 0.4088     [m]       Width of main tires                              
                        w_{t_n} : 0.327      [m]       Width of nose tires                              
                            x_m : 25.63      [m]       x-location of main gear                          
                            x_n : 11.83      [m]       x-location of nose gear                          
                            y_m : 2.846      [m]       y-location of main gear (symmetric)              
                                                                                                        
         VerticalTail, Aircraft |                                                                       
                        A_{fan} : 2.405      [m**2]    Engine reference area                            
                         A_{vt} : 0.176                Vertical tail aspect ratio                       
                    C_{D_{vis}} : 0.005094             Viscous drag coefficient                         
                     C_{L_{vt}} : 0.2347               Vertical tail lift coefficient                   
                         D_{wm} : 3112       [N]       Engine out windmill drag                         
                   L_{max_{vt}} : 1.427e+06  [N]       Maximum load for structural sizing               
                         L_{vt} : 1.312e+04  [N]       Vertical tail lift in engine out                 
                        Re_{vt} : 7.716e+07            Vertical tail reynolds number, cruise            
                              S : 43.2       [m**2]    Vertical tail reference area (full)              
                         S_{vt} : 21.6       [m**2]    Vertical tail ref. area (half)                   
                     W_{struct} : 1.1e+04    [N]       Full span weight                                 
                \Delta x_{lead} : 45.15      [m]       Distance from CG to vertical tail leading edge   
               \Delta x_{trail} : 62.19      [m]       Distance from CG to vertical tail trailing edge  
                   \bar{c}_{vt} : 12.15      [m]       Vertical tail mean aero chord                    
                   \lambda_{vt} : 0.3                  Vertical tail taper ratio                        
                      \tau_{vt} : 0.15                 Vertical tail thickness/chord ratio              
                              b : 3.9        [m]       Span                                             
                  c_{root_{vt}} : 17.04      [m]       Vertical tail root chord                         
                   c_{tip_{vt}} : 5.113      [m]       Vertical tail tip chord                          
                         l_{vt} : 48.63      [m]       Vertical tail moment arm                         
                         p_{vt} : 1.6                  Substituted variable = 1 + 2*taper               
                         q_{vt} : 1.3                  Substituted variable = 1 + taper                 
               z_{\bar{c}_{vt}} : 0.5281     [m]       Vertical location of mean aerodynamic chord      
                                                                                                        
                 Wing, Aircraft |                                                                       
                             AR : 1.196                Wing aspect ratio                                
                        C_{D_w} : 0.0576               Drag coefficient                                 
                    L_{max_{w}} : 3.969e+06  [N]       Maximum load                                     
                       \alpha_w : 0.1                  Wing angle of attack                             
                        \lambda : 0.2                  Wing taper ratio                                 
                         \tau_w : 0.15                 Wing thickness/chord ratio                       
                            b_w : 12.22      [m]       Wing span                                        
                       c_{root} : 17.04      [m]       Wing root chord                                  
                        c_{tip} : 3.408      [m]       Wing tip chord                                   
                            e_w : 0.9961               Oswald efficiency factor                         
                   f(\lambda_w) : 0.0033               Empirical efficiency function of taper           
                            p_w : 1.4                  Substituted variable = 1 + 2*taper               
                            q_w : 1.2                  Substituted variable = 1 + taper                 
                                                                                                        
WingBox, VerticalTail, Aircraft |                                                                       
                              A : 0.3521               Aspect ratio                                     

Constants
---------
                Aircraft |                                                                                
                N_{lift} : 2                      Wing loading multiplier                                 
                     S_w : 125        [m**2]      Wing area                                               
              V_{\infty} : 234        [m/s]       Freestream velocity                                     
                  V_{ne} : 144        [m/s]       Never exceed velocity                                   
                 W_{eng} : 1e+04      [N]         Engine weight                                           
                     \mu : 1.4e-05    [N*s/m**2]  Dynamic viscosity (35,000ft)                            
                    \rho : 0.38       [kg/m**3]   Air density (35,000 ft)                                 
                  \rho_0 : 1.225      [kg/m**3]   Air density (0 ft)                                      
              \rho_{cap} : 2700       [kg/m**3]   Density of spar cap material                            
              \rho_{web} : 2700       [kg/m**3]   Density of shear web material                           
      \sigma_{max,shear} : 1.67e+08   [Pa]        Allowable shear stress                                  
            \sigma_{max} : 2.5e+08    [Pa]        Allowable tensile stress                                
                 d_{fan} : 1.75       [m]         Fan diameter                                            
               f_{w,add} : 0.4                    Wing added weight fraction                              
                       g : 9.81       [m/s**2]    Gravitational acceleration                              
                     r_h : 0.75                   Fractional wing thickness at spar web                   
                       w : 0.5                    Wingbox-width-to-chord ratio                            
            x_{CG_{eng}} : 15         [m]         x-location of engine CG                                 
           x_{CG_{wing}} : 15         [m]         x-location of wing CG                                   
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
                    xapu : 120        [ft]        x-location of APU                                       
                    xfix : 2.1        [m]         x-location of fixed weight                              
                                                                                                          
HorizontalTail, Aircraft |                                                                                
            C_{L_{hmax}} : 2.6                    Max lift coefficient                                    
            C_{m_{fuse}} : 0.05                   Moment coefficient (fuselage)                           
              S.M._{min} : 0.05                   Minimum stability margin                                
              \Delta x_w : 2          [m]         Distance from aerodynamic centre to CG                  
          \alpha_{max,h} : 0.1                    Max angle of attack (htail)                             
                  \eta_h : 0.97                   Lift efficiency (diff between sectional and actual lift)
         \tan(\Lambda_h) : 0.5774                 tangent of horizontal tail sweep                        
            |C_{m_{ac}}| : 0.1                    Moment coefficient about aerodynamic centre (wing)      
                                                                                                          
   LandingGear, Aircraft |                                                                                
                       E : 205        [GPa]       Modulus of elasticity, 4340 steel                       
                       K : 2                      Column effective length factor                          
                     N_s : 2                      Factor of safety                                        
              W_{0_{lg}} : 8.044e+05  [N]         Weight of aircraft excluding landing gear               
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
                                                                                                          
          Wing, Aircraft |                                                                                
             C_{D_{0_w}} : 0.05                   Wing parasitic drag coefficient                         
            C_{L_{wmax}} : 2.5                    Lift coefficient (wing)                                 
          \alpha_{max,w} : 0.1                    Max angle of attack                                     
                  \eta_w : 0.97                   Lift efficiency (diff between sectional and actual lift)
           \tan(\Lambda) : 0.5774                 tangent of wing sweep                                   

Sensitivities
-------------
        Wing, Aircraft |                                                         
           C_{D_{0_w}} : 0.3194   Wing parasitic drag coefficient                
          C_{L_{wmax}} : 0.08678  Lift coefficient (wing)                        
        \alpha_{max,w} : -0.02523 Max angle of attack                            
                                                                                 
VerticalTail, Aircraft |                                                         
                   T_e : 0.03554  Thrust per engine at takeoff                   
          C_{L_{vmax}} : 0.02521  Max lift coefficient                           
                   V_c : 0.0112   Cruise velocity                                
             \rho_{TO} : -0.01033 Air density (SL))                              
            c_{l_{vt}} : -0.01709 Sectional lift force coefficient (engine out)  
                     e : -0.01931 Span efficiency of vertical tail               
                   V_1 : -0.07108 Minimum takeoff velocity                       
                                                                                 
 LandingGear, Aircraft |                                                         
            W_{0_{lg}} : 0.09054  Weight of aircraft excluding landing gear      
             f_{add,m} : 0.02159  Proportional added weight, main                
                     K : 0.01301  Column effective length factor                 
                n_{mg} : -0.0387  Number of main gear struts                     
               n_{wps} : -0.04361 Number of wheels per strut                     
                                                                                 
    Fuselage, Aircraft |                                                         
              n_{seat} : 0.4103    Number of seats                               
                    LF : 0.2091   Load factor                                    
              f_{padd} : 0.1933   Other misc weight as fraction of payload weight
         W_{avg. pass} : 0.1845   Average passenger weight                       
              \Delta h : 0.1166   Distance from floor to widest part of fuselage 
                   p_s : 0.1156   Seat pitch                                     
             W'_{seat} : 0.08553  Weight per seat                                
         \rho_{\infty} : 0.05906  Air density (35,000ft)                         
           \rho_{skin} : 0.0562   Skin density                                   
              \Delta p : 0.0562   Pressure difference across fuselage skin       
               W_{fix} : 0.03894  Fixed weights (pilots, cockpit seats, navcom)  
           W'_{window} : 0.03174  Weight/length density of windows               
           \rho_{cone} : 0.02521  Cone material density                          
           W_{checked} : 0.0246   Ave. checked bag weight                        
           W''_{floor} : 0.01788  Floor weight/area density                      
               f_{apu} : 0.01742  APU weight as fraction of payload weight       
            f_{lugg,1} : 0.0164   Proportion of passengers with one suitcase     
            f_{string} : 0.01583  Fractional weight of stringers                 
              N_{land} : 0.0144   Emergency landing load factor                  
          \rho_{floor} : 0.0144   Floor material density                         
           W''_{insul} : 0.01409  Weight/area density of insulation material     
             W_{cargo} : 0.0138   Cargo weight                                   
             f_{frame} : 0.01131  Fractional frame weight                        
        \sigma_{floor} : -0.01325 Max allowable cap stress                       
         \sigma_{skin} : -0.08141 Max allowable skin stress                      
                   SPR : -0.1156  Number of seats per row                        
                                                                                 
              Aircraft |                                                         
            V_{\infty} : 0.6507   Freestream velocity                            
                   S_w : 0.4059   Wing area                                      
                  \rho : 0.2545   Air density (35,000 ft)                        
                V_{ne} : 0.224    Never exceed velocity                          
                     g : 0.1946   Gravitational acceleration                     
              N_{lift} : 0.08678  Wing loading multiplier                        
                \rho_0 : 0.08678  Air density (0 ft)                             
            \rho_{cap} : 0.05159  Density of spar cap material                   
                   r_h : 0.03509  Fractional wing thickness at spar web          
            \rho_{web} : 0.03509  Density of shear web material                  
               y_{eng} : 0.0326   Engine moment arm                              
               W_{eng} : 0.02978  Engine weight                                  
             f_{w,add} : 0.02477  Wing added weight fraction                     
               d_{fan} : 0.01745  Fan diameter                                   
                   \mu : 0.01561  Dynamic viscosity (35,000ft)                   
    \sigma_{max,shear} : -0.03509 Allowable shear stress                         
          \sigma_{max} : -0.05169 Allowable tensile stress                       

