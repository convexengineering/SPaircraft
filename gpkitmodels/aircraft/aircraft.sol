Beginning signomial solve.
Solving took 22 GP solves and 15.2 seconds.
Warning: Constraint [x_{CG_{vt}}_Aircraft >= 0.5*\D...] is not tight because the left hand side evaluated to 59.084882457 meter but the right hand side evaluated to 55.6692839631 meter (Allowable error: 0.0001%, Actual error: 6.1%)

Warning: Constraint [W_{lg}_Aircraft*x_{CG_{lg}}_Ai...] is not tight because the left hand side evaluated to 330011.289125 meter * newton but the right hand side evaluated to 319089.973057 meter * newton (Allowable error: 1.0%, Actual error: 3.4%)

Warning: Constraint [x_{CG_{ht}}_Aircraft >= 0.5*\D...] is not tight because the left hand side evaluated to 59.9522967927 meter but the right hand side evaluated to 58.9397092401 meter (Allowable error: 0.0001%, Actual error: 1.7%)


Cost
----
 1.91e+05 [N] 

Free Variables
--------------
                         Aircraft |                                                                       
                              C_D : 0.06741              Drag coefficient                                 
                              C_L : 0.3983               Lift coefficient                                 
                          C_{L_w} : 0.3983               Lift coefficient (wing)                          
                       C_{L_{aw}} : 3.983                Lift curve slope (wing)                          
                                D : 1.137e+05  [N]       Total aircraft drag (cruise)                     
                         D_{fuse} : 1.272e+04  [N]       Fuselage drag                                    
                           D_{ht} : 904.3      [N]       Horizontal tail drag                             
                           D_{vt} : 1000       [N]       Vertical tail drag                               
                         D_{wing} : 9.906e+04  [N]       Wing drag                                        
                              L_w : 6.718e+05  [N]       Wing lift                                        
                      L_{v_{max}} : 6.392e+05  [N]       Maximum load for structural sizing               
                              S_w : 162.1      [m**2]    Wing reference area                              
                           V_{TO} : 62.43      [m/s]     Takeoff speed                                    
                                W : 6.718e+05  [N]       Total aircraft weight                            
                         W_{fuel} : 1.91e+05   [N]       Fuel weight                                      
                         W_{fuse} : 1.647e+05  [N]       Fuselage weight                                  
                           W_{ht} : 1928       [N]       Horizontal tail weight                           
                           W_{lg} : 1.622e+04  [N]       Weight of landing gear                           
                          W_{pay} : 1.616e+05  [N]       Payload weight                                   
                           W_{vt} : 1160       [N]       Vertical tail weight                             
                         W_{wing} : 1.251e+05  [N]       Wing weight                                      
                           W_{zf} : 4.807e+05  [N]       Zero fuel weight                                 
                   \bar{c}_{wing} : 6.083      [m]       Mean aerodynamic chord (wing)                    
                      \frac{L}{D} : 5.909                Lift/drag ratio                                  
                              \xi : 0.2022               Takeoff parameter                                
                           b_{vt} : 3.218      [m]       Vertical tail span                               
                           c_{vt} : 9.254      [m]       Vertical tail root chord                         
                         h_{hold} : 0.9042     [m]       Hold height                                      
                         l_{fuse} : 60.3       [m]       Fuselage length                                  
                    p_{\lambda_v} : 1.6                  1 + 2*Tail taper ratio                           
                         w_{fuse} : 3.971      [m]       Fuselage width                                   
                              x_w : 20.47      [m]       Position of wing aerodynamic center              
                     x_{CG_{eng}} : 20.47      [m]       x-location of engine CG                          
                      x_{CG_{fu}} : 54.77      [m]       x-location of fuselage CG                        
                      x_{CG_{ht}} : 59.95      [m]       Horizontal tail CG location                      
                      x_{CG_{lg}} : 20.35      [m]       x-location of landing gear CG                    
                      x_{CG_{vt}} : 59.08      [m]       x-location of tail CG                            
                    x_{CG_{wing}} : 20.47      [m]       x-location of wing CG                            
                           x_{CG} : 18.47      [m]       x-location of CG                                 
                           x_{TO} : 1524       [m]       Takeoff distance                                 
                           x_{up} : 29.61      [m]       Fuselage upsweep point                           
                                y : 0.3721               Takeoff parameter                                
                          z_{bre} : 0.335                Breguet parameter                                
                                                                                                          
               Fuselage, Aircraft |                                                                       
                        A_{floor} : 0.06228    [m**2]    Floor beam x-sectional area                      
                         A_{fuse} : 12.38      [m**2]    Fuselage x-sectional area                        
                         A_{hold} : 2.175      [m**2]    Cargo hold x-sectional area                      
                         A_{skin} : 0.01245    [m**2]    Skin cross sectional area                        
                     D_{friction} : 1.173e+04  [N]       Friction drag                                    
                      D_{upsweep} : 985.5      [N]       Drag due to fuse upsweep                         
                               FF : 1.055                Fuselage form factor                             
                        M_{floor} : 4.876e+05  [N*m]     Max bending moment in floor beams                
                        P_{floor} : 1.137e+06  [N]       Distributed floor load                           
                         R_{fuse} : 1.985      [m]       Fuselage radius                                  
                         S_{bulk} : 24.76      [m**2]    Bulkhead surface area                            
                        S_{floor} : 5.686e+05  [N]       Maximum shear in floor beams                     
                         S_{nose} : 53.65      [m**2]    Nose surface area                                
                         V_{bulk} : 0.02471    [m**3]    Bulkhead skin volume                             
                        V_{cabin} : 361.8      [m**3]    Cabin volume                                     
                        V_{cargo} : 6.796      [m**3]    Cargo volume                                     
                         V_{cone} : 0.1358     [m**3]    Cone skin volume                                 
                          V_{cyl} : 0.3038     [m**3]    Cylinder skin volume                             
                        V_{floor} : 0.2136     [m**3]    Floor volume                                     
                         V_{hold} : 53.1       [m**3]    Hold volume                                      
                         V_{lugg} : 18.24      [m**3]    Luggage volume                                   
                         V_{nose} : 0.05354    [m**3]    Nose skin volume                                 
                          W_{apu} : 5657       [N]       APU weight                                       
                         W_{buoy} : 1743       [N]       Buoyancy weight                                  
                         W_{cone} : 6474       [N]       Cone weight                                      
                        W_{floor} : 1.15e+04   [N]       Floor weight                                     
                        W_{insul} : 4633       [N]       Insulation material weight                       
                         W_{lugg} : 1.79e+04   [N]       Passenger luggage weight                         
                         W_{padd} : 6.465e+04  [N]       Misc weights (galley, toilets, doors etc.)       
                         W_{pass} : 1.337e+05  [N]       Passenger weight                                 
                         W_{seat} : 2.79e+04   [N]       Seating weight                                   
                        W_{shell} : 1.822e+04  [N]       Shell weight                                     
                         W_{skin} : 1.012e+04  [N]       Skin weight                                      
                       W_{window} : 1.062e+04  [N]       Window weight                                    
                   \lambda_{cone} : 0.4                  Tailcone radius taper ratio (xshell2->xtail)     
                             \phi : 0.08321              Upsweep angle                                    
                     \rho_{cabin} : 0.8711     [kg/m**3] Air density in cabin                             
                         \sigma_x : 3.831e+07  [N/m**2]  Axial stress in skin                             
                  \sigma_{\theta} : 1.034e+08  [N/m**2]  Skin hoop stress                                 
                      \tau_{cone} : 1.034e+08  [N/m**2]  Shear stress in cone                             
                                f : 15.19                Fineness ratio                                   
                        h_{floor} : 0.08104    [m]       Floor I-beam height                              
                         l_{cone} : 23.14      [m]       Cone length                                      
                        l_{floor} : 28.38      [m]       Floor length                                     
                         l_{nose} : 5.2        [m]       Nose length                                      
                        l_{shell} : 24.41      [m]       Shell length                                     
                         n_{pass} : 167                  Number of passengers                             
                         n_{rows} : 31                   Number of rows                                   
                        t_{shell} : 0.001347   [m]       Shell thickness                                  
                         t_{skin} : 0.0009979  [m]       Skin thickness                                   
                        w_{floor} : 3.43       [m]       Floor width                                      
                           xVbulk : 0.7317     [m**4]    Volume moment of bulkhead                        
                            xVcyl : 10.39      [m**4]    Volume moment of cylinder                        
                           xVnose : 0.1392     [m**4]    Volume moment of nose                            
                            xWapu : 2.069e+05  [N*m]     Moment of APU                                    
                           xWcone : 6.611e+05  [N*m]     Moment of cone                                   
                            xWfix : 2.802e+04  [N*m]     Moment of fixed weights                          
                          xWfloor : 4.966e+05  [N*m]     Moment of floor weight                           
                           xWfuse : 9.022e+06  [N*m]     Fuselage moment                                  
                          xWinsul : 2.759e+05  [N*m]     Moment of insulation material                    
                           xWpadd : 1.688e+06  [N*m]     Moment of misc weights                           
                           xWseat : 9.095e+05  [N*m]     Moment of seats                                  
                          xWshell : 1.029e+06  [N*m]     Mass moment of shell                             
                           xWskin : 5.717e+05  [N*m]     Mass moment of skin                              
                         xWwindow : 4.711e+05  [N*m]     Mass moment of windows                           
                       x_{shell1} : 5.2        [m]       Start of cylinder section                        
                       x_{shell2} : 29.61      [m]       End of cylinder section                          
                                                                                                          
         HorizontalTail, Aircraft |                                                                       
                             AR_h : 4.637                Horizontal tail aspect ratio                     
                          C_{D_h} : 0.007073             Horizontal tail drag coefficient                 
                      C_{D_{0_h}} : 0.005308             Horizontal tail parasitic drag coefficient       
                          C_{L_h} : 0.1591               Lift coefficient (htail)                         
                       C_{L_{ah}} : 3.702                Lift curve slope (htail)                         
                              K_f : 0.3887               Empirical factor for fuselage-wing interference  
                              L_h : 2.034e+04  [N]       Horizontal tail downforce                        
                      L_{{max}_h} : 4.058e+05  [N]       Maximum load                                     
                         Re_{c_h} : 1.187e+07            Cruise Reynolds number (Horizontal tail)         
                             S.M. : 0.05                 Stability margin                                 
                              S_h : 12.29      [m**2]    Horizontal tail area                             
              \Delta x_{{lead}_h} : 39.11      [m]       Distance from CG to horizontal tail leading edge 
             \Delta x_{{trail}_h} : 41.83      [m]       Distance from CG to horizontal tail trailing edge
                           \alpha : 0.04298              Horizontal tail angle of attack                  
                     \bar{c}_{ht} : 1.869      [m]       Mean aerodynamic chord (ht)                      
                        \lambda_h : 0.2                  Horizontal tail taper ratio                      
                           \tau_h : 0.15                 Horizontal tail thickness/chord ratio            
                           b_{ht} : 7.549      [m]       Horizontal tail span                             
                       c_{root_h} : 2.713      [m]       Horizontal tail root chord                       
                        c_{tip_h} : 0.5426     [m]       Horizontal tail tip chord                        
                              e_h : 0.9849               Oswald efficiency factor                         
                     f(\lambda_h) : 0.0033               Empirical efficiency function of taper           
                              l_h : 40.83      [m]       Horizontal tail moment arm                       
                           p_{ht} : 1.4                  Substituted variable = 1 + 2*taper               
                           q_{ht} : 1.2                  Substituted variable = 1 + taper                 
                 y_{\bar{c}_{ht}} : 2.157      [m]       Vertical location of mean aerodynamic chord      
                                                                                                          
            LandingGear, Aircraft |                                                                       
                                B : 11.31      [m]       Landing gear base                                
                         E_{land} : 3.809e+05  [J]       Max KE to be absorbed in landing                 
                          F_{w_m} : 7161                 Weight factor (main)                             
                          F_{w_n} : 643.8                Weight factor (nose)                             
                              I_m : 7.204e-06  [m**4]    Area moment of inertia (main strut)              
                              I_n : 1.178e-06  [m**4]    Area moment of inertia (nose strut)              
                              L_m : 6.437e+05  [N]       Max static load through main gear                
                              L_n : 1.609e+05  [N]       Min static load through nose gear                
                      L_{n_{dyn}} : 9.654e+04  [N]       Dyn. braking load, nose gear                     
                          L_{w_m} : 1.609e+05  [N]       Static load per wheel (main)                     
                          L_{w_n} : 8.044e+04  [N]       Static load per wheel (nose)                     
                             S_sa : 0.2959     [m]       Stroke of the shock absorber                     
                                T : 5.754      [m]       Main landing gear track                          
                           W_{mg} : 1.471e+04  [N]       Weight of main gear                              
                           W_{ms} : 1410       [N]       Weight of main struts                            
                           W_{mw} : 2377       [N]       Weight of main wheels (per strut)                
                           W_{ng} : 1508       [N]       Weight of nose gear                              
                           W_{ns} : 137.3      [N]       Weight of nose strut                             
                           W_{nw} : 548.2      [N]       Weight of nose wheels (total)                    
                         W_{wa,m} : 267.2      [lbf]     Wheel assembly weight for single main gear wheel 
                         W_{wa,n} : 61.62      [lbf]     Wheel assembly weight for single nose gear wheel 
                       \Delta x_m : 2.262      [m]       Distance b/w main gear and CG                    
                       \Delta x_n : 9.052      [m]       Distance b/w nose gear and CG                    
                       \tan(\phi) : 0.2679               Angle b/w main gear and CG                       
                       \tan(\psi) : 1.963                Tip over angles                                  
                      d_{nacelle} : 2.05       [m]       Nacelle diameter                                 
                         d_{oleo} : 0.3735     [m]       Diameter of oleo shock absorber                  
                          d_{t_m} : 44.5       [in]      Diameter of main gear tires                      
                          d_{t_n} : 35.6       [in]      Diameter of nose gear tires                      
                              l_m : 2.379      [m]       Length of main gear                              
                              l_n : 1.627      [m]       Length of nose gear                              
                         l_{oleo} : 0.7396     [m]       Length of oleo shock absorber                    
                              r_m : 0.04326    [m]       Radius of main gear struts                       
                              r_n : 0.04638    [m]       Radius of nose gear struts                       
                              t_m : 0.02832    [m]       Thickness of main gear strut wall                
                              t_n : 0.003759   [m]       Thickness of nose gear strut wall                
                          w_{t_m} : 0.4088     [m]       Width of main tires                              
                          w_{t_n} : 0.3271     [m]       Width of nose tires                              
                              x_m : 20.73      [m]       x-location of main gear                          
                              x_n : 9.419      [m]       x-location of nose gear                          
                              y_m : 2.877      [m]       y-location of main gear (symmetric)              
                                                                                                          
           VerticalTail, Aircraft |                                                                       
                          A_{fan} : 2.405      [m**2]    Engine reference area                            
                           A_{vt} : 0.535                Vertical tail aspect ratio                       
                      C_{D_{vis}} : 0.004966             Viscous drag coefficient                         
                       C_{L_{vt}} : 0.3645               Vertical tail lift coefficient                   
                           D_{wm} : 3112       [N]       Engine out windmill drag                         
                     L_{max_{vt}} : 1.278e+06  [N]       Maximum load for structural sizing               
                           L_{vt} : 1.826e+04  [N]       Vertical tail lift in engine out                 
                          Re_{vt} : 4.19e+07             Vertical tail reynolds number, cruise            
                                S : 38.71      [m**2]    Vertical tail reference area (full)              
                           S_{vt} : 19.36      [m**2]    Vertical tail ref. area (half)                   
                       W_{struct} : 2320       [N]       Full span weight                                 
                  \Delta x_{lead} : 32.57      [m]       Distance from CG to vertical tail leading edge   
                 \Delta x_{trail} : 41.83      [m]       Distance from CG to vertical tail trailing edge  
                     \bar{c}_{vt} : 6.597      [m]       Vertical tail mean aero chord                    
                     \lambda_{vt} : 0.3                  Vertical tail taper ratio                        
                        \tau_{vt} : 0.1297               Vertical tail thickness/chord ratio              
                                b : 6.436      [m]       Span                                             
                    c_{root_{vt}} : 9.254      [m]       Vertical tail root chord                         
                     c_{tip_{vt}} : 2.776      [m]       Vertical tail tip chord                          
                           l_{vt} : 34.95      [m]       Vertical tail moment arm                         
                           p_{vt} : 1.6                  Substituted variable = 1 + 2*taper               
                           q_{vt} : 1.3                  Substituted variable = 1 + taper                 
                 z_{\bar{c}_{vt}} : 0.8715     [m]       Vertical location of mean aerodynamic chord      
                                                                                                          
                   Wing, Aircraft |                                                                       
                               AR : 5.89                 Wing aspect ratio                                
                          C_{D_w} : 0.05874              Drag coefficient                                 
                      L_{max_{w}} : 5.147e+06  [N]       Maximum load                                     
                         \alpha_w : 0.1                  Wing angle of attack                             
                          \lambda : 0.2                  Wing taper ratio                                 
                           \tau_w : 0.15                 Wing thickness/chord ratio                       
                              b_w : 30.9       [m]       Wing span                                        
                         c_{root} : 8.831      [m]       Wing root chord                                  
                          c_{tip} : 1.766      [m]       Wing tip chord                                   
                              e_w : 0.9809               Oswald efficiency factor                         
                     f(\lambda_w) : 0.0033               Empirical efficiency function of taper           
                              p_w : 1.4                  Substituted variable = 1 + 2*taper               
                              q_w : 1.2                  Substituted variable = 1 + taper                 
                      y_{\bar{c}} : 8.829      [m]       Spanwise location of mean aerodynamic chord      
                                                                                                          
WingBox, HorizontalTail, Aircraft |                                                                       
                          I_{cap} : 8.948e-06            Non-dim spar cap area moment of inertia          
                              M_r : 1.098e+05  [N]       Root moment per root chord                       
                          W_{cap} : 1177       [N]       Weight of spar caps                              
                          W_{web} : 200.9      [N]       Weight of shear web                              
                              \nu : 0.8612               Dummy variable = $(t^2 + t + 1)/(t+1)$           
                          t_{cap} : 0.001934             Non-dim. spar cap thickness                      
                          t_{web} : 0.001467             Non-dim. shear web thickness                     
                                                                                                          
  WingBox, VerticalTail, Aircraft |                                                                       
                                A : 1.07                 Aspect ratio                                     
                          I_{cap} : 5.524e-07            Non-dim spar cap area moment of inertia          
                              M_r : 9.118e+04  [N]       Root moment per root chord                       
                          W_{cap} : 1052       [N]       Weight of spar caps                              
                          W_{web} : 604.6      [N]       Weight of shear web                              
                              \nu : 0.8225               Dummy variable = $(t^2 + t + 1)/(t+1)$           
                          t_{cap} : 0.0001556            Non-dim. spar cap thickness                      
                          t_{web} : 0.0004594            Non-dim. shear web thickness                     
                                                                                                          
          WingBox, Wing, Aircraft |                                                                       
                          I_{cap} : 1.388e-05            Non-dim spar cap area moment of inertia          
                              M_r : 1.769e+06  [N]       Root moment per root chord                       
                          W_{cap} : 7.89e+04   [N]       Weight of spar caps                              
                          W_{web} : 1.043e+04  [N]       Weight of shear web                              
                              \nu : 0.8612               Dummy variable = $(t^2 + t + 1)/(t+1)$           
                          t_{cap} : 0.003051             Non-dim. spar cap thickness                      
                          t_{web} : 0.001792             Non-dim. shear web thickness                     

Constants
---------
                Aircraft |                                                                                
             C_{L_{max}} : 2.5                    Max lift coefficient                                    
                N_{lift} : 2                      Wing loading multiplier                                 
                   Range : 3000       [nmi]       Range                                                   
                    TSFC : 0.3        [lb/hr/lbf] Thrust specific fuel consumption                        
                     T_e : 1.29e+05   [N]         Engine thrust at takeoff                                
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
                     l_r : 5000       [ft]        Runway length                                           
                     r_h : 0.75                   Fractional wing thickness at spar web                   
                       w : 0.5                    Wingbox-width-to-chord ratio                            
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
                  \eta_w : 0.97                   Lift efficiency (diff b/w sectional, actual lift)       
           \tan(\Lambda) : 0.5774                 tangent of wing sweep                                   

Sensitivities
-------------
          Wing, Aircraft |                                                           
             C_{D_{0_w}} : 0.9663   Wing parasitic drag coefficient                  
            C_{L_{wmax}} : 0.3771   Lift coefficient (wing)                          
           \tan(\Lambda) : 0.2424   tangent of wing sweep                            
                  \eta_w : -0.9694  Lift efficiency (diff b/w sectional, actual lift)
          \alpha_{max,w} : -1.347   Max angle of attack                              
                                                                                     
  VerticalTail, Aircraft |                                                           
                     V_c : 0.02287  Cruise velocity                                  
            C_{L_{vmax}} : 0.02208  Max lift coefficient                             
                  \rho_c : 0.01141  Air density (35,000ft)                           
                       e : -0.01212 Span efficiency of vertical tail                 
               \rho_{TO} : -0.02156 Air density (SL))                                
              c_{l_{vt}} : -0.03258 Sectional lift force coefficient (engine out)    
                     V_1 : -0.08729 Minimum takeoff velocity                         
                                                                                     
   LandingGear, Aircraft |                                                           
              W_{0_{lg}} : 0.1132   Weight of aircraft excluding landing gear        
               f_{add,m} : 0.02062  Proportional added weight, main                  
                       K : 0.01631  Column effective length factor                   
       \tan(\theta_{TO}) : 0.01509  Takeoff pitch angle                              
                  n_{mg} : -0.05975 Number of main gear struts                       
                 n_{wps} : -0.06636 Number of wheels per strut                       
                                                                                     
HorizontalTail, Aircraft |                                                           
              \Delta x_w : 0.0182   Distance from aerodynamic centre to CG           
                                                                                     
      Fuselage, Aircraft |                                                           
                n_{seat} : 0.8685    Number of seats                                 
                      LF : 0.642    Load factor                                      
           W_{avg. pass} : 0.5665   Average passenger weight                         
                f_{padd} : 0.187    Other misc weight as fraction of payload weight  
                \Delta h : 0.172    Distance from floor to widest part of fuselage   
                     p_s : 0.1435   Seat pitch                                       
           \rho_{\infty} : 0.1149   Air density (35,000ft)                           
               W'_{seat} : 0.08308  Weight per seat                                  
             W_{checked} : 0.07553  Ave. checked bag weight                          
             \rho_{skin} : 0.05268  Skin density                                     
                \Delta p : 0.05268  Pressure difference across fuselage skin         
              f_{lugg,1} : 0.05035  Proportion of passengers with one suitcase       
               W_{cargo} : 0.04236  Cargo weight                                     
                 W_{fix} : 0.03859  Fixed weights (pilots, cockpit seats, navcom)    
             W'_{window} : 0.0307   Weight/length density of windows                 
              f_{lugg,2} : 0.02518  Proportion of passengers with two suitcases      
             \rho_{cone} : 0.01872  Cone material density                            
             W''_{floor} : 0.01689  Floor weight/area density                        
                N_{land} : 0.01636  Emergency landing load factor                    
            \rho_{floor} : 0.01636  Floor material density                           
                 f_{apu} : 0.01636  APU weight as fraction of payload weight         
              f_{string} : 0.01388  Fractional weight of stringers                   
             W''_{insul} : 0.0134   Weight/area density of insulation material       
          \sigma_{floor} : -0.01528 Max allowable cap stress                         
           \sigma_{skin} : -0.0714  Max allowable skin stress                        
                     SPR : -0.1435  Number of seats per row                          
                                                                                     
                Aircraft |                                                           
                       g : 1.775    Gravitational acceleration                       
                    TSFC : 1.303    Thrust specific fuel consumption                 
                   Range : 1.303    Range                                            
                  V_{ne} : 0.8098   Never exceed velocity                            
                N_{lift} : 0.3862   Wing loading multiplier                          
                  \rho_0 : 0.3828   Air density (0 ft)                               
              \rho_{cap} : 0.3263   Density of spar cap material                     
               f_{w,add} : 0.1059   Wing added weight fraction                       
                     r_h : 0.04425  Fractional wing thickness at spar web            
              \rho_{web} : 0.04425  Density of shear web material                    
                     T_e : 0.04364  Engine thrust at takeoff                         
                 y_{eng} : 0.0428   Engine moment arm                                
                 W_{eng} : 0.02892  Engine weight                                    
                     \mu : 0.02725  Dynamic viscosity (35,000ft)                     
                       w : -0.01564 Wingbox-width-to-chord ratio                     
      \sigma_{max,shear} : -0.04425 Allowable shear stress                           
            \sigma_{max} : -0.3419  Allowable tensile stress                         
                    \rho : -0.5476  Air density (35,000 ft)                          
              V_{\infty} : -2.133   Freestream velocity                              

