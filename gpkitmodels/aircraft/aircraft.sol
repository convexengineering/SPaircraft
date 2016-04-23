Beginning signomial solve.
Solving took 22 GP solves and 11.9 seconds.
Warning: Constraint [x_{CG_{vt}}_Aircraft >= 0.5*\D...] is not tight because the left hand side evaluated to 59.91297746 meter but the right hand side evaluated to 56.4954754037 meter (Allowable error: 0.0001%, Actual error: 6%)

Warning: Constraint [W_{lg}_Aircraft*x_{CG_{lg}}_Ai...] is not tight because the left hand side evaluated to 332548.338714 meter * newton but the right hand side evaluated to 325372.53469 meter * newton (Allowable error: 1.0%, Actual error: 2.2%)

Warning: Constraint [x_{CG_{ht}}_Aircraft >= 0.5*\D...] is not tight because the left hand side evaluated to 60.7781511303 meter but the right hand side evaluated to 59.7577805939 meter (Allowable error: 0.0001%, Actual error: 1.7%)


Cost
----
 1.913e+05 [N] 

Free Variables
--------------
                         Aircraft |                                                                       
                              C_D : 0.06748              Drag coefficient                                 
                              C_L : 0.3983               Lift coefficient                                 
                          C_{L_w} : 0.3983               Lift coefficient (wing)                          
                       C_{L_{aw}} : 3.983                Lift curve slope (wing)                          
                                D : 1.139e+05  [N]       Total aircraft drag (cruise)                     
                         D_{fuse} : 1.285e+04  [N]       Fuselage drag                                    
                           D_{ht} : 909.6      [N]       Horizontal tail drag                             
                           D_{vt} : 999        [N]       Vertical tail drag                               
                         D_{wing} : 9.915e+04  [N]       Wing drag                                        
                              L_w : 6.723e+05  [N]       Wing lift                                        
                      L_{v_{max}} : 6.386e+05  [N]       Maximum load for structural sizing               
                              S_w : 162.2      [m**2]    Wing reference area                              
                                W : 6.723e+05  [N]       Total aircraft weight                            
                         W_{fuel} : 1.913e+05  [N]       Fuel weight                                      
                         W_{fuse} : 1.647e+05  [N]       Fuselage weight                                  
                           W_{ht} : 1947       [N]       Horizontal tail weight                           
                           W_{lg} : 1.627e+04  [N]       Weight of landing gear                           
                          W_{pay} : 1.616e+05  [N]       Payload weight                                   
                           W_{vt} : 1157       [N]       Vertical tail weight                             
                         W_{wing} : 1.252e+05  [N]       Wing weight                                      
                           W_{zf} : 4.809e+05  [N]       Zero fuel weight                                 
                   \bar{c}_{wing} : 6.087      [m]       Mean aerodynamic chord (wing)                    
                      \frac{L}{D} : 5.902                Lift/drag ratio                                  
                           b_{vt} : 3.214      [m]       Vertical tail span                               
                           c_{vt} : 9.258      [m]       Vertical tail root chord                         
                         h_{hold} : 0.9046     [m]       Hold height                                      
                         l_{fuse} : 61.12      [m]       Fuselage length                                  
                    p_{\lambda_v} : 1.6                  1 + 2*Tail taper ratio                           
                         w_{fuse} : 3.971      [m]       Fuselage width                                   
                              x_w : 21.25      [m]       Position of wing aerodynamic center              
                     x_{CG_{eng}} : 21.25      [m]       x-location of engine CG                          
                      x_{CG_{fu}} : 57.28      [m]       x-location of fuselage CG                        
                      x_{CG_{ht}} : 60.78      [m]       Horizontal tail CG location                      
                      x_{CG_{lg}} : 20.44      [m]       x-location of landing gear CG                    
                      x_{CG_{vt}} : 59.91      [m]       x-location of tail CG                            
                    x_{CG_{wing}} : 21.25      [m]       x-location of wing CG                            
                           x_{CG} : 19.25      [m]       x-location of CG                                 
                           x_{up} : 29.61      [m]       Fuselage upsweep point                           
                          z_{bre} : 0.3354               Breguet parameter                                
                                                                                                          
               Fuselage, Aircraft |                                                                       
                        A_{floor} : 0.06239    [m**2]    Floor beam x-sectional area                      
                         A_{fuse} : 12.38      [m**2]    Fuselage x-sectional area                        
                         A_{hold} : 2.177      [m**2]    Cargo hold x-sectional area                      
                         A_{skin} : 0.01245    [m**2]    Skin cross sectional area                        
                     D_{friction} : 1.186e+04  [N]       Friction drag                                    
                      D_{upsweep} : 984.9      [N]       Drag due to fuse upsweep                         
                               FF : 1.055                Fuselage form factor                             
                        M_{floor} : 4.877e+05  [N*m]     Max bending moment in floor beams                
                        P_{floor} : 1.137e+06  [N]       Distributed floor load                           
                         R_{fuse} : 1.985      [m]       Fuselage radius                                  
                         S_{bulk} : 24.77      [m**2]    Bulkhead surface area                            
                        S_{floor} : 5.686e+05  [N]       Maximum shear in floor beams                     
                         S_{nose} : 53.66      [m**2]    Nose surface area                                
                         V_{bulk} : 0.02472    [m**3]    Bulkhead skin volume                             
                        V_{cabin} : 361.9      [m**3]    Cabin volume                                     
                        V_{cargo} : 6.796      [m**3]    Cargo volume                                     
                         V_{cone} : 0.1355     [m**3]    Cone skin volume                                 
                          V_{cyl} : 0.3039     [m**3]    Cylinder skin volume                             
                        V_{floor} : 0.214      [m**3]    Floor volume                                     
                         V_{hold} : 53.13      [m**3]    Hold volume                                      
                         V_{lugg} : 18.24      [m**3]    Luggage volume                                   
                         V_{nose} : 0.05355    [m**3]    Nose skin volume                                 
                          W_{apu} : 5657       [N]       APU weight                                       
                         W_{buoy} : 1744       [N]       Buoyancy weight                                  
                         W_{cone} : 6462       [N]       Cone weight                                      
                        W_{floor} : 1.151e+04  [N]       Floor weight                                     
                        W_{insul} : 4634       [N]       Insulation material weight                       
                         W_{lugg} : 1.79e+04   [N]       Passenger luggage weight                         
                         W_{padd} : 6.465e+04  [N]       Misc weights (galley, toilets, doors etc.)       
                         W_{pass} : 1.337e+05  [N]       Passenger weight                                 
                         W_{seat} : 2.79e+04   [N]       Seating weight                                   
                        W_{shell} : 1.822e+04  [N]       Shell weight                                     
                         W_{skin} : 1.012e+04  [N]       Skin weight                                      
                       W_{window} : 1.062e+04  [N]       Window weight                                    
                   \lambda_{cone} : 0.4                  Tailcone radius taper ratio (xshell2->xtail)     
                             \phi : 0.08319              Upsweep angle                                    
                     \rho_{cabin} : 0.8711     [kg/m**3] Air density in cabin                             
                         \sigma_x : 3.831e+07  [N/m**2]  Axial stress in skin                             
                  \sigma_{\theta} : 1.034e+08  [N/m**2]  Skin hoop stress                                 
                      \tau_{cone} : 1.034e+08  [N/m**2]  Shear stress in cone                             
                                f : 15.39                Fineness ratio                                   
                        h_{floor} : 0.0809     [m]       Floor I-beam height                              
                         l_{cone} : 23.15      [m]       Cone length                                      
                        l_{floor} : 28.38      [m]       Floor length                                     
                         l_{nose} : 5.2        [m]       Nose length                                      
                        l_{shell} : 24.41      [m]       Shell length                                     
                         n_{pass} : 167                  Number of passengers                             
                         n_{rows} : 31                   Number of rows                                   
                        t_{shell} : 0.001347   [m]       Shell thickness                                  
                         t_{skin} : 0.000998   [m]       Skin thickness                                   
                        w_{floor} : 3.431      [m]       Floor width                                      
                           xVbulk : 0.732      [m**4]    Volume moment of bulkhead                        
                            xVcyl : 11.37      [m**4]    Volume moment of cylinder                        
                           xVnose : 0.1392     [m**4]    Volume moment of nose                            
                            xWapu : 2.069e+05  [N*m]     Moment of APU                                    
                           xWcone : 7.002e+05  [N*m]     Moment of cone                                   
                            xWfix : 2.802e+04  [N*m]     Moment of fixed weights                          
                          xWfloor : 5.316e+05  [N*m]     Moment of floor weight                           
                           xWfuse : 9.437e+06  [N*m]     Fuselage moment                                  
                          xWinsul : 3.094e+05  [N*m]     Moment of insulation material                    
                           xWpadd : 1.719e+06  [N*m]     Moment of misc weights                           
                           xWseat : 9.434e+05  [N*m]     Moment of seats                                  
                          xWshell : 1.151e+06  [N*m]     Mass moment of shell                             
                           xWskin : 6.393e+05  [N*m]     Mass moment of skin                              
                         xWwindow : 5.058e+05  [N*m]     Mass moment of windows                           
                       x_{shell1} : 5.2        [m]       Start of cylinder section                        
                       x_{shell2} : 29.61      [m]       End of cylinder section                          
                                                                                                          
         HorizontalTail, Aircraft |                                                                       
                             AR_h : 4.617                Horizontal tail aspect ratio                     
                          C_{D_h} : 0.00704              Horizontal tail drag coefficient                 
                      C_{D_{0_h}} : 0.005307             Horizontal tail parasitic drag coefficient       
                          C_{L_h} : 0.1573               Lift coefficient (htail)                         
                       C_{L_{ah}} : 3.696                Lift curve slope (htail)                         
                              K_f : 0.4015               Empirical factor for fuselage-wing interference  
                              L_h : 2.033e+04  [N]       Horizontal tail downforce                        
                      L_{{max}_h} : 4.101e+05  [N]       Maximum load                                     
                         Re_{c_h} : 1.196e+07            Cruise Reynolds number (Horizontal tail)         
                             S.M. : 0.05                 Stability margin                                 
                              S_h : 12.42      [m**2]    Horizontal tail area                             
              \Delta x_{{lead}_h} : 39.15      [m]       Distance from CG to horizontal tail leading edge 
             \Delta x_{{trail}_h} : 41.88      [m]       Distance from CG to horizontal tail trailing edge
                           \alpha : 0.04257              Horizontal tail angle of attack                  
                     \bar{c}_{ht} : 1.883      [m]       Mean aerodynamic chord (ht)                      
                        \lambda_h : 0.2                  Horizontal tail taper ratio                      
                           \tau_h : 0.15                 Horizontal tail thickness/chord ratio            
                           b_{ht} : 7.572      [m]       Horizontal tail span                             
                       c_{root_h} : 2.734      [m]       Horizontal tail root chord                       
                        c_{tip_h} : 0.5467     [m]       Horizontal tail tip chord                        
                              e_h : 0.985                Oswald efficiency factor                         
                     f(\lambda_h) : 0.0033               Empirical efficiency function of taper           
                              l_h : 40.87      [m]       Horizontal tail moment arm                       
                           p_{ht} : 1.4                  Substituted variable = 1 + 2*taper               
                           q_{ht} : 1.2                  Substituted variable = 1 + taper                 
                 y_{\bar{c}_{ht}} : 2.164      [m]       Vertical location of mean aerodynamic chord      
                                                                                                          
            LandingGear, Aircraft |                                                                       
                                B : 7.16       [m]       Landing gear base                                
                         E_{land} : 3.809e+05  [J]       Max KE to be absorbed in landing                 
                          F_{w_m} : 7161                 Weight factor (main)                             
                          F_{w_n} : 643.8                Weight factor (nose)                             
                              I_m : 7.289e-06  [m**4]    Area moment of inertia (main strut)              
                              I_n : 1.271e-06  [m**4]    Area moment of inertia (nose strut)              
                              L_m : 6.437e+05  [N]       Max static load through main gear                
                              L_n : 1.609e+05  [N]       Min static load through nose gear                
                      L_{n_{dyn}} : 1.53e+05   [N]       Dyn. braking load, nose gear                     
                          L_{w_m} : 1.609e+05  [N]       Static load per wheel (main)                     
                          L_{w_n} : 8.044e+04  [N]       Static load per wheel (nose)                     
                             S_sa : 0.2959     [m]       Stroke of the shock absorber                     
                                T : 6.078      [m]       Main landing gear track                          
                           W_{mg} : 1.473e+04  [N]       Weight of main gear                              
                           W_{ms} : 1422       [N]       Weight of main struts                            
                           W_{mw} : 2377       [N]       Weight of main wheels (per strut)                
                           W_{ng} : 1538       [N]       Weight of nose gear                              
                           W_{ns} : 167.4      [N]       Weight of nose strut                             
                           W_{nw} : 548.2      [N]       Weight of nose wheels (total)                    
                         W_{wa,m} : 267.2      [lbf]     Wheel assembly weight for single main gear wheel 
                         W_{wa,n} : 61.62      [lbf]     Wheel assembly weight for single nose gear wheel 
                       \Delta x_m : 1.432      [m]       Distance b/w main gear and CG                    
                       \Delta x_n : 5.729      [m]       Distance b/w nose gear and CG                    
                       \tan(\phi) : 0.2679               Angle b/w main gear and CG                       
                       \tan(\psi) : 1.963                Tip over angles                                  
                      d_{nacelle} : 2.05       [m]       Nacelle diameter                                 
                         d_{oleo} : 0.3735     [m]       Diameter of oleo shock absorber                  
                          d_{t_m} : 44.5       [in]      Diameter of main gear tires                      
                          d_{t_n} : 35.6       [in]      Diameter of nose gear tires                      
                              l_m : 2.393      [m]       Length of main gear                              
                              l_n : 1.627      [m]       Length of nose gear                              
                         l_{oleo} : 0.7397     [m]       Length of oleo shock absorber                    
                              r_m : 0.04346    [m]       Radius of main gear struts                       
                              r_n : 0.04362    [m]       Radius of nose gear struts                       
                              t_m : 0.02826    [m]       Thickness of main gear strut wall                
                              t_n : 0.004874   [m]       Thickness of nose gear strut wall                
                          w_{t_m} : 0.4088     [m]       Width of main tires                              
                          w_{t_n} : 0.3271     [m]       Width of nose tires                              
                              x_m : 20.68      [m]       x-location of main gear                          
                              x_n : 13.52      [m]       x-location of nose gear                          
                              y_m : 3.039      [m]       y-location of main gear (symmetric)              
                                                                                                          
           VerticalTail, Aircraft |                                                                       
                          A_{fan} : 2.405      [m**2]    Engine reference area                            
                           A_{vt} : 0.534                Vertical tail aspect ratio                       
                      C_{D_{vis}} : 0.004965             Viscous drag coefficient                         
                       C_{L_{vt}} : 0.3643               Vertical tail lift coefficient                   
                           D_{wm} : 3112       [N]       Engine out windmill drag                         
                     L_{max_{vt}} : 1.277e+06  [N]       Maximum load for structural sizing               
                           L_{vt} : 1.823e+04  [N]       Vertical tail lift in engine out                 
                          Re_{vt} : 4.192e+07            Vertical tail reynolds number, cruise            
                                S : 38.68      [m**2]    Vertical tail reference area (full)              
                           S_{vt} : 19.34      [m**2]    Vertical tail ref. area (half)                   
                       W_{struct} : 2313       [N]       Full span weight                                 
                  \Delta x_{lead} : 32.62      [m]       Distance from CG to vertical tail leading edge   
                 \Delta x_{trail} : 41.88      [m]       Distance from CG to vertical tail trailing edge  
                     \bar{c}_{vt} : 6.6        [m]       Vertical tail mean aero chord                    
                     \lambda_{vt} : 0.3                  Vertical tail taper ratio                        
                        \tau_{vt} : 0.1296               Vertical tail thickness/chord ratio              
                                b : 6.427      [m]       Span                                             
                    c_{root_{vt}} : 9.258      [m]       Vertical tail root chord                         
                     c_{tip_{vt}} : 2.777      [m]       Vertical tail tip chord                          
                           l_{vt} : 35         [m]       Vertical tail moment arm                         
                           p_{vt} : 1.6                  Substituted variable = 1 + 2*taper               
                           q_{vt} : 1.3                  Substituted variable = 1 + taper                 
                 z_{\bar{c}_{vt}} : 0.8703     [m]       Vertical location of mean aerodynamic chord      
                                                                                                          
                   Wing, Aircraft |                                                                       
                               AR : 5.889                Wing aspect ratio                                
                          C_{D_w} : 0.05874              Drag coefficient                                 
                      L_{max_{w}} : 5.151e+06  [N]       Maximum load                                     
                         \alpha_w : 0.1                  Wing angle of attack                             
                          \lambda : 0.2                  Wing taper ratio                                 
                           \tau_w : 0.15                 Wing thickness/chord ratio                       
                              b_w : 30.91      [m]       Wing span                                        
                         c_{root} : 8.836      [m]       Wing root chord                                  
                          c_{tip} : 1.767      [m]       Wing tip chord                                   
                              e_w : 0.9809               Oswald efficiency factor                         
                     f(\lambda_w) : 0.0033               Empirical efficiency function of taper           
                              p_w : 1.4                  Substituted variable = 1 + 2*taper               
                              q_w : 1.2                  Substituted variable = 1 + taper                 
                      y_{\bar{c}} : 8.831      [m]       Spanwise location of mean aerodynamic chord      
                                                                                                          
WingBox, HorizontalTail, Aircraft |                                                                       
                          I_{cap} : 8.868e-06            Non-dim spar cap area moment of inertia          
                              M_r : 1.105e+05  [N]       Root moment per root chord                       
                          W_{cap} : 1187       [N]       Weight of spar caps                              
                          W_{web} : 203.6      [N]       Weight of shear web                              
                              \nu : 0.8612               Dummy variable = $(t^2 + t + 1)/(t+1)$           
                          t_{cap} : 0.001916             Non-dim. spar cap thickness                      
                          t_{web} : 0.001461             Non-dim. shear web thickness                     
                                                                                                          
  WingBox, VerticalTail, Aircraft |                                                                       
                                A : 1.068                Aspect ratio                                     
                          I_{cap} : 5.501e-07            Non-dim spar cap area moment of inertia          
                              M_r : 9.094e+04  [N]       Root moment per root chord                       
                          W_{cap} : 1049       [N]       Weight of spar caps                              
                          W_{web} : 603.3      [N]       Weight of shear web                              
                              \nu : 0.8225               Dummy variable = $(t^2 + t + 1)/(t+1)$           
                          t_{cap} : 0.0001551            Non-dim. spar cap thickness                      
                          t_{web} : 0.0004589            Non-dim. shear web thickness                     
                                                                                                          
          WingBox, Wing, Aircraft |                                                                       
                          I_{cap} : 1.387e-05            Non-dim spar cap area moment of inertia          
                              M_r : 1.769e+06  [N]       Root moment per root chord                       
                          W_{cap} : 7.896e+04  [N]       Weight of spar caps                              
                          W_{web} : 1.044e+04  [N]       Weight of shear web                              
                              \nu : 0.8612               Dummy variable = $(t^2 + t + 1)/(t+1)$           
                          t_{cap} : 0.003048             Non-dim. spar cap thickness                      
                          t_{web} : 0.001791             Non-dim. shear web thickness                     

Constants
---------
                         Aircraft |                                                                                
                            Range : 3000       [nmi]       Range                                                   
                             TSFC : 0.3        [lb/hr/lbf] Thrust specific fuel consumption                        
                       V_{\infty} : 234        [m/s]       Freestream velocity                                     
                           V_{ne} : 144        [m/s]       Never exceed velocity                                   
                          W_{eng} : 1e+04      [N]         Engine weight                                           
                              \mu : 1.4e-05    [N*s/m**2]  Dynamic viscosity (35,000ft)                            
                             \rho : 0.38       [kg/m**3]   Air density                                             
                           \rho_0 : 1.225      [kg/m**3]   Air density (0 ft)                                      
                          d_{fan} : 1.75       [m]         Fan diameter                                            
                                g : 9.81       [m/s**2]    Gravitational acceleration                              
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
                                g : 9.81       [m/s**2]    Gravitational acceleration                              
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
                                g : 9.81       [m/s**2]    Gravitational acceleration                              
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
                           \eta_w : 0.97                   Lift efficiency (diff b/w sectional, actual lift)       
                    \tan(\Lambda) : 0.5774                 tangent of wing sweep                                   
                                                                                                                   
WingBox, HorizontalTail, Aircraft |                                                                                
                         N_{lift} : 2                      Wing loading multiplier                                 
                       \rho_{cap} : 2700       [kg/m**3]   Density of spar cap material                            
                       \rho_{web} : 2700       [kg/m**3]   Density of shear web material                           
               \sigma_{max,shear} : 1.67e+08   [Pa]        Allowable shear stress                                  
                     \sigma_{max} : 2.5e+08    [Pa]        Allowable tensile stress                                
                        f_{w,add} : 0.4                    Wing added weight fraction                              
                                g : 9.81       [m/s**2]    Gravitational acceleration                              
                              r_h : 0.75                   Fractional wing thickness at spar web                   
                                w : 0.5                    Wingbox-width-to-chord ratio                            
                                                                                                                   
  WingBox, VerticalTail, Aircraft |                                                                                
                         N_{lift} : 2                      Wing loading multiplier                                 
                       \rho_{cap} : 2700       [kg/m**3]   Density of spar cap material                            
                       \rho_{web} : 2700       [kg/m**3]   Density of shear web material                           
               \sigma_{max,shear} : 1.67e+08   [Pa]        Allowable shear stress                                  
                     \sigma_{max} : 2.5e+08    [Pa]        Allowable tensile stress                                
                        f_{w,add} : 0.4                    Wing added weight fraction                              
                                g : 9.81       [m/s**2]    Gravitational acceleration                              
                              r_h : 0.75                   Fractional wing thickness at spar web                   
                                w : 0.5                    Wingbox-width-to-chord ratio                            
                                                                                                                   
          WingBox, Wing, Aircraft |                                                                                
                         N_{lift} : 2                      Wing loading multiplier                                 
                       \rho_{cap} : 2700       [kg/m**3]   Density of spar cap material                            
                       \rho_{web} : 2700       [kg/m**3]   Density of shear web material                           
               \sigma_{max,shear} : 1.67e+08   [Pa]        Allowable shear stress                                  
                     \sigma_{max} : 2.5e+08    [Pa]        Allowable tensile stress                                
                        f_{w,add} : 0.4                    Wing added weight fraction                              
                                g : 9.81       [m/s**2]    Gravitational acceleration                              
                              r_h : 0.75                   Fractional wing thickness at spar web                   
                                w : 0.5                    Wingbox-width-to-chord ratio                            

Sensitivities
-------------
 WingBox, Wing, Aircraft |                                                           
                N_{lift} : 0.3768   Wing loading multiplier                          
                       g : 0.3614   Gravitational acceleration                       
              \rho_{cap} : 0.3192   Density of spar cap material                     
               f_{w,add} : 0.1033   Wing added weight fraction                       
                     r_h : 0.0422   Fractional wing thickness at spar web            
              \rho_{web} : 0.0422   Density of shear web material                    
                       w : -0.01547 Wingbox-width-to-chord ratio                     
      \sigma_{max,shear} : -0.0422  Allowable shear stress                           
            \sigma_{max} : -0.3346  Allowable tensile stress                         
                                                                                     
          Wing, Aircraft |                                                           
             C_{D_{0_w}} : 0.9652   Wing parasitic drag coefficient                  
            C_{L_{wmax}} : 0.3768   Lift coefficient (wing)                          
           \tan(\Lambda) : 0.2421   tangent of wing sweep                            
                  \eta_w : -0.9683  Lift efficiency (diff b/w sectional, actual lift)
          \alpha_{max,w} : -1.345   Max angle of attack                              
                                                                                     
  VerticalTail, Aircraft |                                                           
                     T_e : 0.04349  Thrust per engine at takeoff                     
                     V_c : 0.0228   Cruise velocity                                  
            C_{L_{vmax}} : 0.022    Max lift coefficient                             
                  \rho_c : 0.01138  Air density (35,000ft)                           
                       e : -0.01209 Span efficiency of vertical tail                 
               \rho_{TO} : -0.02148 Air density (SL))                                
              c_{l_{vt}} : -0.03245 Sectional lift force coefficient (engine out)    
                     V_1 : -0.08697 Minimum takeoff velocity                         
                                                                                     
   LandingGear, Aircraft |                                                           
              W_{0_{lg}} : 0.1135   Weight of aircraft excluding landing gear        
               f_{add,m} : 0.02059  Proportional added weight, main                  
                       K : 0.01643  Column effective length factor                   
       \tan(\theta_{TO}) : 0.01577  Takeoff pitch angle                              
                  n_{mg} : -0.05987 Number of main gear struts                       
                 n_{wps} : -0.06654 Number of wheels per strut                       
                                                                                     
HorizontalTail, Aircraft |                                                           
              \Delta x_w : 0.01815  Distance from aerodynamic centre to CG           
                                                                                     
      Fuselage, Aircraft |                                                           
                n_{seat} : 0.869     Number of seats                                 
                      LF : 0.6411   Load factor                                      
           W_{avg. pass} : 0.5657   Average passenger weight                         
                f_{padd} : 0.1867   Other misc weight as fraction of payload weight  
                \Delta h : 0.1725   Distance from floor to widest part of fuselage   
                     p_s : 0.1449   Seat pitch                                       
           \rho_{\infty} : 0.1159   Air density (35,000ft)                           
                       g : 0.09267  Gravitational acceleration                       
               W'_{seat} : 0.08297  Weight per seat                                  
             W_{checked} : 0.07542  Ave. checked bag weight                          
             \rho_{skin} : 0.05261  Skin density                                     
                \Delta p : 0.05261  Pressure difference across fuselage skin         
              f_{lugg,1} : 0.05028  Proportion of passengers with one suitcase       
               W_{cargo} : 0.0423   Cargo weight                                     
                 W_{fix} : 0.03853  Fixed weights (pilots, cockpit seats, navcom)    
             W'_{window} : 0.03066  Weight/length density of windows                 
              f_{lugg,2} : 0.02514  Proportion of passengers with two suitcases      
             \rho_{cone} : 0.01866  Cone material density                            
             W''_{floor} : 0.01687  Floor weight/area density                        
                N_{land} : 0.01637  Emergency landing load factor                    
            \rho_{floor} : 0.01637  Floor material density                           
                 f_{apu} : 0.01633  APU weight as fraction of payload weight         
              f_{string} : 0.01386  Fractional weight of stringers                   
             W''_{insul} : 0.01338  Weight/area density of insulation material       
          \sigma_{floor} : -0.01529 Max allowable cap stress                         
           \sigma_{skin} : -0.07127 Max allowable skin stress                        
                     SPR : -0.1449  Number of seats per row                          
                                                                                     
                Aircraft |                                                           
                       g : 1.303    Gravitational acceleration                       
                    TSFC : 1.303    Thrust specific fuel consumption                 
                   Range : 1.303    Range                                            
                  V_{ne} : 0.8092   Never exceed velocity                            
                  \rho_0 : 0.3826   Air density (0 ft)                               
                 y_{eng} : 0.0427   Engine moment arm                                
                 W_{eng} : 0.02887  Engine weight                                    
                     \mu : 0.02749  Dynamic viscosity (35,000ft)                     
                    \rho : -0.5472  Air density                                      
              V_{\infty} : -2.13    Freestream velocity                              

