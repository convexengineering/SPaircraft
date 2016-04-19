Beginning signomial solve.
Solving took 7 GP solves and 222 seconds.
Warning: Constraint x_{CG_{vt}}_Aircraft >= 0.5*\Delta x_{lead}_VerticalTail, Aircraft + 0.5*\Delta x_{trail}_VerticalTail, Aircraft + x_{CG}_Aircraft is not tight because x_{CG_{vt}}_Aircraft [m] evaluated to 57.3981205496 meter but 0.5*\Delta x_{lead}_VerticalTail, Aircraft + 0.5*\Delta x_{trail}_VerticalTail, Aircraft + x_{CG}_Aircraft [m] evaluated to 56.5452797045 meter
Warning: Constraint x_{CG_{ht}}_Aircraft >= 0.5*\Delta x_{{lead}_h}_HorizontalTail, Aircraft + 0.5*\Delta x_{{trail}_h}_HorizontalTail, Aircraft + x_{CG}_Aircraft is not tight because x_{CG_{ht}}_Aircraft [m] evaluated to 56.9744099543 meter but 0.5*\Delta x_{{lead}_h}_HorizontalTail, Aircraft + 0.5*\Delta x_{{trail}_h}_HorizontalTail, Aircraft + x_{CG}_Aircraft [m] evaluated to 56.5452797028 meter

Cost
----
 1.601e+05 [N] 

Free Variables
--------------
                       Aircraft |                                                                       
                        C_{L_w} : 0.347                Lift coefficient (wing)                          
                     C_{L_{aw}} : 3.47                 Lift curve slope (wing)                          
                              D : 9.862e+04  [N]       Total aircraft drag (cruise)                     
                       D_{fuse} : 1.29e+04   [N]       Fuselage drag                                    
                         D_{ht} : 1719       [N]       Horizontal tail drag                             
                         D_{vt} : 1037       [N]       Vertical tail drag                               
                       D_{wing} : 8.297e+04  [N]       Wing drag                                        
                        I_{cap} : 6.06e-06             Non-dim spar cap area moment of inertia          
                    L_{v_{max}} : 6.421e+05  [N]       Maximum load for structural sizing               
                            M_r : 9.584e+05  [N]       Root moment per root chord                       
                            S_w : 133        [m**2]    Wing area                                        
                              W : 4.801e+05  [N]       Total aircraft weight                            
                        W_{cap} : 3.066e+04  [N]       Weight of spar caps                              
                      W_{empty} : 3.199e+05  [N]       Empty aircraft weight                            
                       W_{fuel} : 1.601e+05  [N]       Fuel weight                                      
                       W_{fuse} : 1.647e+05  [N]       Fuselage weight                                  
                         W_{ht} : 5.174e+04  [N]       Horizontal tail weight                           
                         W_{lg} : 1.585e+04  [N]       Weight of landing gear                           
                         W_{vt} : 2.587e+04  [N]       Vertical tail weight                             
                        W_{web} : 6297       [N]       Weight of shear web                              
                       W_{wing} : 5.174e+04  [N]       Wing weight                                      
                    \frac{L}{D} : 4.868                Lift/drag ratio                                  
                            \nu : 0.8612               Dummy variable = $(t^2 + t + 1)/(t+1)$           
                         b_{vt} : 3.071      [m]       Vertical tail span                               
                         c_{vt} : 9.742      [m]       Vertical tail root chord                         
                       h_{hold} : 0.9103     [m]       Hold height                                      
                       l_{fuse} : 61.42      [m]       Fuselage length                                  
                  p_{\lambda_v} : 1.6                  1 + 2*Tail taper ratio                           
                        t_{cap} : 0.001297             Non-dim. spar cap thickness                      
                        t_{web} : 0.001184             Non-dim. shear web thickness                     
                       w_{fuse} : 4.002      [m]       Fuselage width                                   
                    x_{CG_{fu}} : 18.69      [m]       x-location of fuselage CG                        
                    x_{CG_{ht}} : 56.97      [m]       Horizontal tail CG location                      
                    x_{CG_{lg}} : 20.44      [m]       x-location of landing gear CG                    
                    x_{CG_{vt}} : 57.4       [m]       x-location of tail CG                            
                         x_{CG} : 18.34      [m]       x-location of CG                                 
                         x_{up} : 29.61      [m]       Fuselage upsweep point                           
                        z_{bre} : 0.4066               Breguet parameter                                
                                                                                                        
             Fuselage, Aircraft |                                                                       
                      A_{floor} : 0.05656    [m**2]    Floor beam x-sectional area                      
                       A_{fuse} : 12.58      [m**2]    Fuselage x-sectional area                        
                       A_{hold} : 2.213      [m**2]    Cargo hold x-sectional area                      
                       A_{skin} : 0.01265    [m**2]    Skin cross sectional area                        
                   D_{friction} : 1.2e+04    [N]       Friction drag                                    
                    D_{upsweep} : 901.9      [N]       Drag due to fuse upsweep                         
                             FF : 1.055                Fuselage form factor                             
                      M_{floor} : 4.928e+05  [N*m]     Max bending moment in floor beams                
                      P_{floor} : 1.137e+06  [N]       Distributed floor load                           
                       R_{fuse} : 2.001      [m]       Fuselage radius                                  
                       S_{bulk} : 25.16      [m**2]    Bulkhead surface area                            
                      S_{floor} : 5.686e+05  [N]       Maximum shear in floor beams                     
                       S_{nose} : 54.12      [m**2]    Nose surface area                                
                       V_{bulk} : 0.02531    [m**3]    Bulkhead skin volume                             
                      V_{cabin} : 367.8      [m**3]    Cabin volume                                     
                      V_{cargo} : 6.796      [m**3]    Cargo volume                                     
                       V_{cone} : 0.136      [m**3]    Cone skin volume                                 
                        V_{cyl} : 0.3087     [m**3]    Cylinder skin volume                             
                      V_{floor} : 0.1961     [m**3]    Floor volume                                     
                       V_{hold} : 54.01      [m**3]    Hold volume                                      
                       V_{lugg} : 18.24      [m**3]    Luggage volume                                   
                       V_{nose} : 0.05444    [m**3]    Nose skin volume                                 
                        W_{apu} : 5657       [N]       APU weight                                       
                       W_{buoy} : 1772       [N]       Buoyancy weight                                  
                       W_{cone} : 6482       [N]       Cone weight                                      
                      W_{floor} : 1.11e+04   [N]       Floor weight                                     
                      W_{insul} : 4673       [N]       Insulation material weight                       
                       W_{lugg} : 1.79e+04   [N]       Passenger luggage weight                         
                       W_{padd} : 6.465e+04  [N]       Misc weights (galley, toilets, doors etc.)       
                       W_{pass} : 1.337e+05  [N]       Passenger weight                                 
                        W_{pay} : 1.616e+05  [N]       Payload weight                                   
                       W_{seat} : 2.79e+04   [N]       Seating weight                                   
                      W_{shell} : 1.852e+04  [N]       Shell weight                                     
                       W_{skin} : 1.029e+04  [N]       Skin weight                                      
                     W_{window} : 1.062e+04  [N]       Window weight                                    
                 \lambda_{cone} : 0.4                  Tailcone radius taper ratio (xshell2->xtail)     
                           \phi : 0.07981              Upsweep angle                                    
                   \rho_{cabin} : 0.8711     [kg/m**3] Air density in cabin                             
                       \sigma_x : 3.831e+07  [N/m**2]  Axial stress in skin                             
                \sigma_{\theta} : 1.034e+08  [N/m**2]  Skin hoop stress                                 
                    \tau_{cone} : 1.034e+08  [N/m**2]  Shear stress in cone                             
                              f : 15.35                Fineness ratio                                   
                      h_{floor} : 0.09084    [m]       Floor I-beam height                              
                       l_{cone} : 24.35      [m]       Cone length                                      
                      l_{floor} : 28.41      [m]       Floor length                                     
                       l_{nose} : 5.2        [m]       Nose length                                      
                      l_{shell} : 24.41      [m]       Shell length                                     
                       n_{pass} : 167                  Number of passengers                             
                       n_{rows} : 31                   Number of rows                                   
                      t_{shell} : 0.001358   [m]       Shell thickness                                  
                       t_{skin} : 0.001006   [m]       Skin thickness                                   
                      w_{floor} : 3.467      [m]       Floor width                                      
                         xVbulk : 0.7494     [m**4]    Volume moment of bulkhead                        
                          xVcyl : 5.656      [m**4]    Volume moment of cylinder                        
                         xVnose : 0.1415     [m**4]    Volume moment of nose                            
                          xWapu : 2.069e+05  [N*m]     Moment of APU                                    
                         xWcone : 3.092e+05  [N*m]     Moment of cone                                   
                          xWfix : 2.802e+04  [N*m]     Moment of fixed weights                          
                        xWfloor : 2.068e+05  [N*m]     Moment of floor weight                           
                         xWfuse : 3.079e+06  [N*m]     Fuselage moment                                  
                        xWinsul : 9.363e+04  [N*m]     Moment of insulation material                    
                         xWpadd : 1.14e+06   [N*m]     Moment of misc weights                           
                         xWseat : 4.998e+05  [N*m]     Moment of seats                                  
                        xWshell : 3.315e+05  [N*m]     Mass moment of shell                             
                         xWskin : 1.841e+05  [N*m]     Mass moment of skin                              
                       xWwindow : 1.983e+05  [N*m]     Mass moment of windows                           
                     x_{shell1} : 5.2        [m]       Start of cylinder section                        
                     x_{shell2} : 29.61      [m]       End of cylinder section                          
                                                                                                        
       HorizontalTail, Aircraft |                                                                       
                           AR_h : 0.9327               Horizontal tail aspect ratio                     
                        C_{D_h} : 0.005185             Horizontal tail drag coefficient                 
                    C_{D_{0_h}} : 0.00513              Horizontal tail parasitic drag coefficient       
                        C_{L_h} : 0.01273              Lift coefficient (htail)                         
                     C_{L_{ah}} : 1.367                Lift curve slope (htail)                         
                            K_f : 0.376                Empirical factor for fuselage-wing interference  
                    L_{{max}_h} : 1.052e+06  [N]       Maximum load                                     
                       Re_{c_h} : 4.262e+07            Cruise Reynolds number (Horizontal tail)         
                           S.M. : 0.05                 Stability margin                                 
                            S_h : 31.87      [m**2]    Horizontal tail area                             
            \Delta x_{{lead}_h} : 33.33      [m]       Distance from CG to horizontal tail leading edge 
           \Delta x_{{trail}_h} : 43.07      [m]       Distance from CG to horizontal tail trailing edge
                         \alpha : 0.009312             Horizontal tail angle of attack                  
                   \bar{c}_{ht} : 6.711      [m]       Mean aerodynamic chord (ht)                      
                 \bar{c}_{wing} : 11.75      [m]       Mean aerodynamic chord (wing)                    
                      \lambda_h : 0.2                  Horizontal tail taper ratio                      
                         \tau_h : 0.15                 Horizontal tail thickness/chord ratio            
                         b_{ht} : 5.452      [m]       Horizontal tail span                             
                     c_{root_h} : 9.742      [m]       Horizontal tail root chord                       
                      c_{tip_h} : 1.948      [m]       Horizontal tail tip chord                        
                            e_h : 0.9969               Oswald efficiency factor                         
                   f(\lambda_h) : 0.0033               Empirical efficiency function of taper           
                            l_h : 35.91      [m]       Horizontal tail moment arm                       
                         p_{ht} : 1.4                  Substituted variable = 1 + 2*taper               
                         q_{ht} : 1.2                  Substituted variable = 1 + taper                 
                            x_w : 20.34      [m]       Position of wing aerodynamic center              
               y_{\bar{c}_{ht}} : 1.558      [m]       Vertical location of mean aerodynamic chord      
                                                                                                        
          LandingGear, Aircraft |                                                                       
                              B : 11.95      [m]       Landing gear base                                
                       E_{land} : 3.809e+05  [J]       Max KE to be absorbed in landing                 
                        F_{w_m} : 7160                 Weight factor (main)                             
                        F_{w_n} : 643.8                Weight factor (nose)                             
                            I_m : 7.198e-06  [m**4]    Area moment of inertia (main strut)              
                            I_n : 1.156e-06  [m**4]    Area moment of inertia (nose strut)              
                            L_m : 6.437e+05  [N]       Max static load through main gear                
                            L_n : 1.609e+05  [N]       Min static load through nose gear                
                    L_{n_{dyn}} : 9.14e+04   [N]       Dyn. braking load, nose gear                     
                        L_{w_m} : 1.609e+05  [N]       Static load per wheel (main)                     
                        L_{w_n} : 8.044e+04  [N]       Static load per wheel (nose)                     
                           S_sa : 0.2959     [m]       Stroke of the shock absorber                     
                              T : 5.734      [m]       Main landing gear track                          
                         W_{mg} : 1.434e+04  [N]       Weight of main gear                              
                         W_{ms} : 1228       [N]       Weight of main struts                            
                         W_{mw} : 2377       [N]       Weight of main wheels (per strut)                
                         W_{ng} : 1505       [N]       Weight of nose gear                              
                         W_{ns} : 134.5      [N]       Weight of nose strut                             
                         W_{nw} : 548.2      [N]       Weight of nose wheels (total)                    
                       W_{wa,m} : 267.2      [lbf]     Wheel assembly weight for single main gear wheel 
                       W_{wa,n} : 61.62      [lbf]     Wheel assembly weight for single nose gear wheel 
                     \Delta x_m : 2.389      [m]       Distance b/w main gear and CG                    
                     \Delta x_n : 9.559      [m]       Distance b/w nose gear and CG                    
                     \tan(\phi) : 0.2679               Angle b/w main gear and CG                       
                     \tan(\psi) : 1.963                Tip over angles                                  
                    d_{nacelle} : 2.05       [m]       Nacelle diameter                                 
                       d_{oleo} : 0.3735     [m]       Diameter of oleo shock absorber                  
                        d_{t_m} : 44.5       [in]      Diameter of main gear tires                      
                        d_{t_n} : 35.6       [in]      Diameter of nose gear tires                      
                            l_m : 2.378      [m]       Length of main gear                              
                            l_n : 1.627      [m]       Length of nose gear                              
                       l_{oleo} : 0.7397     [m]       Length of oleo shock absorber                    
                            r_m : 0.04633    [m]       Radius of main gear struts                       
                            r_n : 0.04641    [m]       Radius of nose gear struts                       
                            t_m : 0.02303    [m]       Thickness of main gear strut wall                
                            t_n : 0.003681   [m]       Thickness of nose gear strut wall                
                        w_{t_m} : 0.4088     [m]       Width of main tires                              
                        w_{t_n} : 0.3271     [m]       Width of nose tires                              
                            x_m : 20.73      [m]       x-location of main gear                          
                            x_n : 8.787      [m]       x-location of nose gear                          
                            y_m : 2.867      [m]       y-location of main gear (symmetric)              
                                                                                                        
         VerticalTail, Aircraft |                                                                       
                        A_{fan} : 2.405      [m**2]    Engine reference area                            
                         A_{vt} : 0.485                Vertical tail aspect ratio                       
                    C_{D_{vis}} : 0.005127             Viscous drag coefficient                         
                     C_{L_{vt}} : 0.3546               Vertical tail lift coefficient                   
                         D_{wm} : 3112       [N]       Engine out windmill drag                         
                   L_{max_{vt}} : 1.284e+06  [N]       Maximum load for structural sizing               
                         L_{vt} : 1.784e+04  [N]       Vertical tail lift in engine out                 
                        Re_{vt} : 4.411e+07            Vertical tail reynolds number, cruise            
                              S : 38.89      [m**2]    Vertical tail reference area (full)              
                         S_{vt} : 19.45      [m**2]    Vertical tail ref. area (half)                   
                     W_{struct} : 5.174e+04  [N]       Full span weight                                 
                \Delta x_{lead} : 33.33      [m]       Distance from CG to vertical tail leading edge   
               \Delta x_{trail} : 43.07      [m]       Distance from CG to vertical tail trailing edge  
                   \bar{c}_{vt} : 6.944      [m]       Vertical tail mean aero chord                    
                   \lambda_{vt} : 0.3                  Vertical tail taper ratio                        
                      \tau_{vt} : 0.15                 Vertical tail thickness/chord ratio              
                              b : 6.142      [m]       Span                                             
                  c_{root_{vt}} : 9.742      [m]       Vertical tail root chord                         
                   c_{tip_{vt}} : 2.923      [m]       Vertical tail tip chord                          
                         l_{vt} : 35.76      [m]       Vertical tail moment arm                         
                         p_{vt} : 1.6                  Substituted variable = 1 + 2*taper               
                         q_{vt} : 1.3                  Substituted variable = 1 + taper                 
               z_{\bar{c}_{vt}} : 0.8317     [m]       Vertical location of mean aerodynamic chord      
                                                                                                        
                 Wing, Aircraft |                                                                       
                             AR : 3.892                Wing aspect ratio                                
                        C_{D_w} : 0.05998              Drag coefficient                                 
                    L_{max_{w}} : 4.222e+06  [N]       Maximum load                                     
                       \alpha_w : 0.1                  Wing angle of attack                             
                        \lambda : 0.2                  Wing taper ratio                                 
                         \tau_w : 0.15                 Wing thickness/chord ratio                       
                            b_w : 22.75      [m]       Wing span                                        
                       c_{root} : 9.742      [m]       Wing root chord                                  
                        c_{tip} : 1.948      [m]       Wing tip chord                                   
                            e_w : 0.9873               Oswald efficiency factor                         
                   f(\lambda_w) : 0.0033               Empirical efficiency function of taper           
                            p_w : 1.4                  Substituted variable = 1 + 2*taper               
                            q_w : 1.2                  Substituted variable = 1 + taper                 
                                                                                                        
WingBox, VerticalTail, Aircraft |                                                                       
                              A : 0.97                 Aspect ratio                                     

Constants
---------
                Aircraft |                                                                                
                N_{lift} : 2                      Wing loading multiplier                                 
                   Range : 3000       [nmi]       Range                                                   
                    TSFC : 0.3        [lb/hr/lbf] Thrust specific fuel consumption                        
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
             C_{D_{0_w}} : 1.129    Wing parasitic drag coefficient                         
            C_{L_{wmax}} : 0.8178   Lift coefficient (wing)                                 
           \tan(\Lambda) : 0.3223   tangent of wing sweep                                   
                  \eta_w : -1.289   Lift efficiency (diff between sectional and actual lift)
          \alpha_{max,w} : -2.147   Max angle of attack                                     
                                                                                            
  VerticalTail, Aircraft |                                                                  
                     T_e : 0.07272  Thrust per engine at takeoff                            
            C_{L_{vmax}} : 0.04033  Max lift coefficient                                    
                     V_c : 0.0336   Cruise velocity                                         
                  \rho_c : 0.01668  Air density (35,000ft)                                  
                       e : -0.02166 Span efficiency of vertical tail                        
               \rho_{TO} : -0.03239 Air density (SL))                                       
              c_{l_{vt}} : -0.05281 Sectional lift force coefficient (engine out)           
                     V_1 : -0.1454  Minimum takeoff velocity                                
                                                                                            
   LandingGear, Aircraft |                                                                  
              W_{0_{lg}} : 0.2116   Weight of aircraft excluding landing gear               
               f_{add,m} : 0.04437  Proportional added weight, main                         
                       K : 0.03056  Column effective length factor                          
       \tan(\theta_{TO}) : 0.02118  Takeoff pitch angle                                     
               \rho_{st} : 0.01612  Density of 4340 Steel                                   
                       E : -0.01528 Modulus of elasticity, 4340 steel                       
                  n_{mg} : -0.1011  Number of main gear struts                              
                 n_{wps} : -0.113   Number of wheels per strut                              
                                                                                            
HorizontalTail, Aircraft |                                                                  
              \Delta x_w : 0.01196  Distance from aerodynamic centre to CG                  
                                                                                            
      Fuselage, Aircraft |                                                                  
                n_{seat} : 0.8941    Number of seats                                        
                      LF : 0.436    Load factor                                             
                f_{padd} : 0.4022   Other misc weight as fraction of payload weight         
           W_{avg. pass} : 0.3847   Average passenger weight                                
                \Delta h : 0.2968   Distance from floor to widest part of fuselage          
                     p_s : 0.2798   Seat pitch                                              
               W'_{seat} : 0.1783   Weight per seat                                         
           \rho_{\infty} : 0.1629   Air density (35,000ft)                                  
             \rho_{skin} : 0.1152   Skin density                                            
                \Delta p : 0.1152   Pressure difference across fuselage skin                
                 W_{fix} : 0.08302  Fixed weights (pilots, cockpit seats, navcom)           
             W'_{window} : 0.06606  Weight/length density of windows                        
             W_{checked} : 0.0513   Ave. checked bag weight                                 
             \rho_{cone} : 0.04033  Cone material density                                   
             W''_{floor} : 0.03677  Floor weight/area density                               
                 f_{apu} : 0.0352   APU weight as fraction of payload weight                
              f_{lugg,1} : 0.0342   Proportion of passengers with one suitcase              
                N_{land} : 0.03231  Emergency landing load factor                           
            \rho_{floor} : 0.03231  Floor material density                                  
              f_{string} : 0.03025  Fractional weight of stringers                          
             W''_{insul} : 0.02907  Weight/area density of insulation material              
               W_{cargo} : 0.02877  Cargo weight                                            
               f_{frame} : 0.02161  Fractional frame weight                                 
               p_{cabin} : 0.01955  Cabin air pressure (8,000ft)                            
                f_{fadd} : 0.01728  Fractional added weight of local reinforcements         
              f_{lugg,2} : 0.0171   Proportion of passengers with two suitcases             
               T_{cabin} : -0.01955 Cabin temperature                                       
                       R : -0.01955 Universal gas constant                                  
          \sigma_{floor} : -0.02996 Max allowable cap stress                                
           \sigma_{skin} : -0.1556  Max allowable skin stress                               
                     SPR : -0.2798  Number of seats per row                                 
                                                                                            
                Aircraft |                                                                  
                       g : 2.629    Gravitational acceleration                              
                  V_{ne} : 1.716    Never exceed velocity                                   
                    TSFC : 1.61     Thrust specific fuel consumption                        
                   Range : 1.61     Range                                                   
                N_{lift} : 0.8178   Wing loading multiplier                                 
                  \rho_0 : 0.8178   Air density (0 ft)                                      
              \rho_{cap} : 0.6677   Density of spar cap material                            
               f_{w,add} : 0.2299   Wing added weight fraction                              
                     r_h : 0.1371   Fractional wing thickness at spar web                   
              \rho_{web} : 0.1371   Density of shear web material                           
                 y_{eng} : 0.06956  Engine moment arm                                       
                 W_{eng} : 0.06221  Engine weight                                           
                     \mu : 0.03985  Dynamic viscosity (35,000ft)                            
                 d_{fan} : 0.02386  Fan diameter                                            
                       w : -0.01304 Wingbox-width-to-chord ratio                            
      \sigma_{max,shear} : -0.1371  Allowable shear stress                                  
            \sigma_{max} : -0.6807  Allowable tensile stress                                
                    \rho : -1.218   Air density (35,000 ft)                                 
              V_{\infty} : -3.664   Freestream velocity                                     

