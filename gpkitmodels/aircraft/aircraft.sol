Beginning signomial solve.
Solving took 22 GP solves and 11 seconds.
Warning: Constraint x_{CG_{vt}}_Aircraft >= 0.5*\Delta x_{lead}_VerticalTail, Aircraft + 0.5*\Delta x_{trail}_VerticalTail, Aircraft + x_{CG}_Aircraft is not tight because x_{CG_{vt}}_Aircraft [m] evaluated to 57.7220358915 meter but 0.5*\Delta x_{lead}_VerticalTail, Aircraft + 0.5*\Delta x_{trail}_VerticalTail, Aircraft + x_{CG}_Aircraft [m] evaluated to 54.3500377989 meter
Warning: Constraint W_{lg}_Aircraft*x_{CG_{lg}}_Aircraft >= W_{mg}_LandingGear, Aircraft*x_m_LandingGear, Aircraft + W_{ng}_LandingGear, Aircraft*x_n_LandingGear, Aircraft is not tight because W_{lg}_Aircraft*x_{CG_{lg}}_Aircraft [N*m] evaluated to 333610.481634 meter * newton but W_{mg}_LandingGear, Aircraft*x_m_LandingGear, Aircraft + W_{ng}_LandingGear, Aircraft*x_n_LandingGear, Aircraft [N*m] evaluated to 327047.969394 meter * newton
Warning: Constraint x_{CG_{ht}}_Aircraft >= 0.5*\Delta x_{{lead}_h}_HorizontalTail, Aircraft + 0.5*\Delta x_{{trail}_h}_HorizontalTail, Aircraft + x_{CG}_Aircraft is not tight because x_{CG_{ht}}_Aircraft [m] evaluated to 58.7085177234 meter but 0.5*\Delta x_{{lead}_h}_HorizontalTail, Aircraft + 0.5*\Delta x_{{trail}_h}_HorizontalTail, Aircraft + x_{CG}_Aircraft [m] evaluated to 58.0437013651 meter
Warning: Constraint 0.667*\lambda_Wing, Aircraft**2*c_{root}_Wing, Aircraft*q_w_Wing, Aircraft**-1 + 0.667*\lambda_Wing, Aircraft*c_{root}_Wing, Aircraft*q_w_Wing, Aircraft**-1 + 0.667*c_{root}_Wing, Aircraft*q_w_Wing, Aircraft**-1 >= \bar{c}_{wing}_Aircraft is not tight because 0.667*\lambda_Wing, Aircraft**2*c_{root}_Wing, Aircraft*q_w_Wing, Aircraft**-1 + 0.667*\lambda_Wing, Aircraft*c_{root}_Wing, Aircraft*q_w_Wing, Aircraft**-1 + 0.667*c_{root}_Wing, Aircraft*q_w_Wing, Aircraft**-1 [m] evaluated to 4.03051541707 meter but \bar{c}_{wing}_Aircraft [m] evaluated to 0.000273698205522 meter

Cost
----
 1.077e+05 [N] 

Free Variables
--------------
                         Aircraft |                                                                       
                              C_D : 0.07453              Drag coefficient                                 
                              C_L : 0.413                Lift coefficient                                 
                          C_{L_w} : 0.413                Lift coefficient (wing)                          
                       C_{L_{aw}} : 4.13                 Lift curve slope (wing)                          
                                D : 6.48e+04   [N]       Total aircraft drag (cruise)                     
                         D_{fuse} : 1.252e+04  [N]       Fuselage drag                                    
                           D_{ht} : 639.6      [N]       Horizontal tail drag                             
                           D_{vt} : 1050       [N]       Vertical tail drag                               
                         D_{wing} : 5.059e+04  [N]       Wing drag                                        
                              L_w : 3.591e+05  [N]       Wing lift                                        
                      L_{v_{max}} : 6.691e+05  [N]       Maximum load for structural sizing               
                              S_w : 83.57      [m**2]    Wing reference area                              
                                W : 3.591e+05  [N]       Total aircraft weight                            
                        W_{empty} : 2.514e+05  [N]       Empty aircraft weight                            
                         W_{fuel} : 1.077e+05  [N]       Fuel weight                                      
                         W_{fuse} : 1.654e+05  [N]       Fuselage weight                                  
                           W_{ht} : 1192       [N]       Horizontal tail weight                           
                           W_{lg} : 1.632e+04  [N]       Weight of landing gear                           
                           W_{vt} : 1324       [N]       Vertical tail weight                             
                         W_{wing} : 5.714e+04  [N]       Wing weight                                      
                   \bar{c}_{wing} : 0.0002737  [m]       Mean aerodynamic chord (wing)                    
                      \frac{L}{D} : 5.541                Lift/drag ratio                                  
                           b_{vt} : 3.401      [m]       Vertical tail span                               
                           c_{vt} : 9.166      [m]       Vertical tail root chord                         
                         h_{hold} : 0.9042     [m]       Hold height                                      
                         l_{fuse} : 58.93      [m]       Fuselage length                                  
                    p_{\lambda_v} : 1.6                  1 + 2*Tail taper ratio                           
                         w_{fuse} : 3.968      [m]       Fuselage width                                   
                              x_w : 21.35      [m]       Position of wing aerodynamic center              
                     x_{CG_{eng}} : 21.35      [m]       x-location of engine CG                          
                      x_{CG_{fu}} : 30.09      [m]       x-location of fuselage CG                        
                      x_{CG_{ht}} : 58.71      [m]       Horizontal tail CG location                      
                      x_{CG_{lg}} : 20.44      [m]       x-location of landing gear CG                    
                      x_{CG_{vt}} : 57.72      [m]       x-location of tail CG                            
                    x_{CG_{wing}} : 21.35      [m]       x-location of wing CG                            
                           x_{CG} : 19.35      [m]       x-location of CG                                 
                           x_{up} : 29.61      [m]       Fuselage upsweep point                           
                          z_{bre} : 0.3572               Breguet parameter                                
                                                                                                          
               Fuselage, Aircraft |                                                                       
                        A_{floor} : 0.06299    [m**2]    Floor beam x-sectional area                      
                         A_{fuse} : 12.37      [m**2]    Fuselage x-sectional area                        
                         A_{hold} : 2.174      [m**2]    Cargo hold x-sectional area                      
                         A_{skin} : 0.01243    [m**2]    Skin cross sectional area                        
                     D_{friction} : 1.152e+04  [N]       Friction drag                                    
                      D_{upsweep} : 1006       [N]       Drag due to fuse upsweep                         
                               FF : 1.055                Fuselage form factor                             
                        M_{floor} : 4.872e+05  [N*m]     Max bending moment in floor beams                
                        P_{floor} : 1.137e+06  [N]       Distributed floor load                           
                         R_{fuse} : 1.984      [m]       Fuselage radius                                  
                         S_{bulk} : 24.74      [m**2]    Bulkhead surface area                            
                        S_{floor} : 5.686e+05  [N]       Maximum shear in floor beams                     
                         S_{nose} : 53.62      [m**2]    Nose surface area                                
                         V_{bulk} : 0.02467    [m**3]    Bulkhead skin volume                             
                        V_{cabin} : 361.4      [m**3]    Cabin volume                                     
                        V_{cargo} : 6.796      [m**3]    Cargo volume                                     
                         V_{cone} : 0.1489     [m**3]    Cone skin volume                                 
                          V_{cyl} : 0.3035     [m**3]    Cylinder skin volume                             
                        V_{floor} : 0.2159     [m**3]    Floor volume                                     
                         V_{hold} : 53.06      [m**3]    Hold volume                                      
                         V_{lugg} : 18.24      [m**3]    Luggage volume                                   
                         V_{nose} : 0.05348    [m**3]    Nose skin volume                                 
                          W_{apu} : 5657       [N]       APU weight                                       
                         W_{buoy} : 1741       [N]       Buoyancy weight                                  
                         W_{cone} : 7097       [N]       Cone weight                                      
                        W_{floor} : 1.155e+04  [N]       Floor weight                                     
                        W_{insul} : 4630       [N]       Insulation material weight                       
                         W_{lugg} : 1.79e+04   [N]       Passenger luggage weight                         
                         W_{padd} : 6.465e+04  [N]       Misc weights (galley, toilets, doors etc.)       
                         W_{pass} : 1.337e+05  [N]       Passenger weight                                 
                          W_{pay} : 1.616e+05  [N]       Payload weight                                   
                         W_{seat} : 2.79e+04   [N]       Seating weight                                   
                        W_{shell} : 1.82e+04   [N]       Shell weight                                     
                         W_{skin} : 1.011e+04  [N]       Skin weight                                      
                       W_{window} : 1.062e+04  [N]       Window weight                                    
                   \lambda_{cone} : 0.4                  Tailcone radius taper ratio (xshell2->xtail)     
                             \phi : 0.08394              Upsweep angle                                    
                     \rho_{cabin} : 0.8711     [kg/m**3] Air density in cabin                             
                         \sigma_x : 3.831e+07  [N/m**2]  Axial stress in skin                             
                  \sigma_{\theta} : 1.034e+08  [N/m**2]  Skin hoop stress                                 
                      \tau_{cone} : 1.034e+08  [N/m**2]  Shear stress in cone                             
                                f : 14.85                Fineness ratio                                   
                        h_{floor} : 0.08       [m]       Floor I-beam height                              
                         l_{cone} : 22.92      [m]       Cone length                                      
                        l_{floor} : 28.38      [m]       Floor length                                     
                         l_{nose} : 5.2        [m]       Nose length                                      
                        l_{shell} : 24.41      [m]       Shell length                                     
                         n_{pass} : 167                  Number of passengers                             
                         n_{rows} : 31                   Number of rows                                   
                        t_{shell} : 0.001346   [m]       Shell thickness                                  
                         t_{skin} : 0.0009974  [m]       Skin thickness                                   
                        w_{floor} : 3.428      [m]       Floor width                                      
                           xVbulk : 0.7305     [m**4]    Volume moment of bulkhead                        
                            xVcyl : 7.78       [m**4]    Volume moment of cylinder                        
                           xVnose : 0.139      [m**4]    Volume moment of nose                            
                            xWapu : 2.069e+05  [N*m]     Moment of APU                                    
                           xWcone : 4.673e+05  [N*m]     Moment of cone                                   
                            xWfix : 2.802e+04  [N*m]     Moment of fixed weights                          
                          xWfloor : 3.29e+05   [N*m]     Moment of floor weight                           
                           xWfuse : 4.976e+06  [N*m]     Fuselage moment                                  
                          xWinsul : 1.748e+05  [N*m]     Moment of insulation material                    
                           xWpadd : 1.314e+06  [N*m]     Moment of misc weights                           
                           xWseat : 6.47e+05   [N*m]     Moment of seats                                  
                          xWshell : 6.024e+05  [N*m]     Mass moment of shell                             
                           xWskin : 3.346e+05  [N*m]     Mass moment of skin                              
                         xWwindow : 3.094e+05  [N*m]     Mass moment of windows                           
                       x_{shell1} : 5.2        [m]       Start of cylinder section                        
                       x_{shell2} : 29.61      [m]       End of cylinder section                          
                                                                                                          
         HorizontalTail, Aircraft |                                                                       
                             AR_h : 6.056                Horizontal tail aspect ratio                     
                          C_{D_h} : 0.008911             Horizontal tail drag coefficient                 
                      C_{D_{0_h}} : 0.005414             Horizontal tail parasitic drag coefficient       
                          C_{L_h} : 0.2554               Lift coefficient (htail)                         
                       C_{L_{ah}} : 4.013                Lift curve slope (htail)                         
                              K_f : 0.425                Empirical factor for fuselage-wing interference  
                              L_h : 1.833e+04  [N]       Horizontal tail downforce                        
                      L_{{max}_h} : 2.278e+05  [N]       Maximum load                                     
                         Re_{c_h} : 7.784e+06            Cruise Reynolds number (Horizontal tail)         
                             S.M. : 0.05                 Stability margin                                 
                              S_h : 6.899      [m**2]    Horizontal tail area                             
              \Delta x_{{lead}_h} : 37.81      [m]       Distance from CG to horizontal tail leading edge 
             \Delta x_{{trail}_h} : 39.58      [m]       Distance from CG to horizontal tail trailing edge
                           \alpha : 0.06364              Horizontal tail angle of attack                  
                     \bar{c}_{ht} : 1.225      [m]       Mean aerodynamic chord (ht)                      
                        \lambda_h : 0.2                  Horizontal tail taper ratio                      
                           \tau_h : 0.15                 Horizontal tail thickness/chord ratio            
                           b_{ht} : 6.464      [m]       Horizontal tail span                             
                       c_{root_h} : 1.779      [m]       Horizontal tail root chord                       
                        c_{tip_h} : 0.3558     [m]       Horizontal tail tip chord                        
                              e_h : 0.9804               Oswald efficiency factor                         
                     f(\lambda_h) : 0.0033               Empirical efficiency function of taper           
                              l_h : 39.18      [m]       Horizontal tail moment arm                       
                           p_{ht} : 1.4                  Substituted variable = 1 + 2*taper               
                           q_{ht} : 1.2                  Substituted variable = 1 + taper                 
                 y_{\bar{c}_{ht}} : 1.847      [m]       Vertical location of mean aerodynamic chord      
                                                                                                          
            LandingGear, Aircraft |                                                                       
                                B : 6.549      [m]       Landing gear base                                
                         E_{land} : 3.809e+05  [J]       Max KE to be absorbed in landing                 
                          F_{w_m} : 7161                 Weight factor (main)                             
                          F_{w_n} : 643.8                Weight factor (nose)                             
                              I_m : 7.321e-06  [m**4]    Area moment of inertia (main strut)              
                              I_n : 1.293e-06  [m**4]    Area moment of inertia (nose strut)              
                              L_m : 6.437e+05  [N]       Max static load through main gear                
                              L_n : 1.609e+05  [N]       Min static load through nose gear                
                      L_{n_{dyn}} : 1.675e+05  [N]       Dyn. braking load, nose gear                     
                          L_{w_m} : 1.609e+05  [N]       Static load per wheel (main)                     
                          L_{w_n} : 8.044e+04  [N]       Static load per wheel (nose)                     
                             S_sa : 0.2959     [m]       Stroke of the shock absorber                     
                                T : 6.197      [m]       Main landing gear track                          
                           W_{mg} : 1.478e+04  [N]       Weight of main gear                              
                           W_{ms} : 1445       [N]       Weight of main struts                            
                           W_{mw} : 2377       [N]       Weight of main wheels (per strut)                
                           W_{ng} : 1546       [N]       Weight of nose gear                              
                           W_{ns} : 175.1      [N]       Weight of nose strut                             
                           W_{nw} : 548.2      [N]       Weight of nose wheels (total)                    
                         W_{wa,m} : 267.2      [lbf]     Wheel assembly weight for single main gear wheel 
                         W_{wa,n} : 61.62      [lbf]     Wheel assembly weight for single nose gear wheel 
                       \Delta x_m : 1.31       [m]       Distance b/w main gear and CG                    
                       \Delta x_n : 5.24       [m]       Distance b/w nose gear and CG                    
                       \tan(\phi) : 0.2679               Angle b/w main gear and CG                       
                       \tan(\psi) : 1.963                Tip over angles                                  
                      d_{nacelle} : 2.05       [m]       Nacelle diameter                                 
                         d_{oleo} : 0.3735     [m]       Diameter of oleo shock absorber                  
                          d_{t_m} : 44.5       [in]      Diameter of main gear tires                      
                          d_{t_n} : 35.6       [in]      Diameter of nose gear tires                      
                              l_m : 2.399      [m]       Length of main gear                              
                              l_n : 1.627      [m]       Length of nose gear                              
                         l_{oleo} : 0.7397     [m]       Length of oleo shock absorber                    
                              r_m : 0.04326    [m]       Radius of main gear struts                       
                              r_n : 0.04301    [m]       Radius of nose gear struts                       
                              t_m : 0.02878    [m]       Thickness of main gear strut wall                
                              t_n : 0.005171   [m]       Thickness of nose gear strut wall                
                          w_{t_m} : 0.4088     [m]       Width of main tires                              
                          w_{t_n} : 0.3271     [m]       Width of nose tires                              
                              x_m : 20.66      [m]       x-location of main gear                          
                              x_n : 14.11      [m]       x-location of nose gear                          
                              y_m : 3.099      [m]       y-location of main gear (symmetric)              
                                                                                                          
           VerticalTail, Aircraft |                                                                       
                          A_{fan} : 2.405      [m**2]    Engine reference area                            
                           A_{vt} : 0.5707               Vertical tail aspect ratio                       
                      C_{D_{vis}} : 0.004981             Viscous drag coefficient                         
                       C_{L_{vt}} : 0.3708               Vertical tail lift coefficient                   
                           D_{wm} : 3112       [N]       Engine out windmill drag                         
                     L_{max_{vt}} : 1.338e+06  [N]       Maximum load for structural sizing               
                           L_{vt} : 1.944e+04  [N]       Vertical tail lift in engine out                 
                          Re_{vt} : 4.15e+07             Vertical tail reynolds number, cruise            
                                S : 40.52      [m**2]    Vertical tail reference area (full)              
                           S_{vt} : 20.26      [m**2]    Vertical tail ref. area (half)                   
                       W_{struct} : 2649       [N]       Full span weight                                 
                  \Delta x_{lead} : 30.42      [m]       Distance from CG to vertical tail leading edge   
                 \Delta x_{trail} : 39.58      [m]       Distance from CG to vertical tail trailing edge  
                     \bar{c}_{vt} : 6.534      [m]       Vertical tail mean aero chord                    
                     \lambda_{vt} : 0.3                  Vertical tail taper ratio                        
                        \tau_{vt} : 0.1317               Vertical tail thickness/chord ratio              
                                b : 6.801      [m]       Span                                             
                    c_{root_{vt}} : 9.166      [m]       Vertical tail root chord                         
                     c_{tip_{vt}} : 2.75       [m]       Vertical tail tip chord                          
                           l_{vt} : 32.82      [m]       Vertical tail moment arm                         
                           p_{vt} : 1.6                  Substituted variable = 1 + 2*taper               
                           q_{vt} : 1.3                  Substituted variable = 1 + taper                 
                 z_{\bar{c}_{vt}} : 0.921      [m]       Vertical location of mean aerodynamic chord      
                                                                                                          
                   Wing, Aircraft |                                                                       
                               AR : 6.782                Wing aspect ratio                                
                          C_{D_w} : 0.05819              Drag coefficient                                 
                      L_{max_{w}} : 2.654e+06  [N]       Maximum load                                     
                         \alpha_w : 0.1                  Wing angle of attack                             
                          \lambda : 0.2                  Wing taper ratio                                 
                           \tau_w : 0.15                 Wing thickness/chord ratio                       
                              b_w : 23.81      [m]       Wing span                                        
                         c_{root} : 5.851      [m]       Wing root chord                                  
                          c_{tip} : 1.17       [m]       Wing tip chord                                   
                              e_w : 0.9781               Oswald efficiency factor                         
                     f(\lambda_w) : 0.0033               Empirical efficiency function of taper           
                              p_w : 1.4                  Substituted variable = 1 + 2*taper               
                              q_w : 1.2                  Substituted variable = 1 + taper                 
                      y_{\bar{c}} : 6.802      [m]       Spanwise location of mean aerodynamic chord      
                                                                                                          
WingBox, HorizontalTail, Aircraft |                                                                       
                          I_{cap} : 1.526e-05            Non-dim spar cap area moment of inertia          
                              M_r : 8.048e+04  [N]       Root moment per root chord                       
                          W_{cap} : 754.7      [N]       Weight of spar caps                              
                          W_{web} : 96.55      [N]       Weight of shear web                              
                              \nu : 0.8612               Dummy variable = $(t^2 + t + 1)/(t+1)$           
                          t_{cap} : 0.003369             Non-dim. spar cap thickness                      
                          t_{web} : 0.001916             Non-dim. shear web thickness                     
                                                                                                          
  WingBox, VerticalTail, Aircraft |                                                                       
                                A : 1.141                Aspect ratio                                     
                          I_{cap} : 6.386e-07            Non-dim spar cap area moment of inertia          
                              M_r : 1.018e+05  [N]       Root moment per root chord                       
                          W_{cap} : 1223       [N]       Weight of spar caps                              
                          W_{web} : 668.8      [N]       Weight of shear web                              
                              \nu : 0.8225               Dummy variable = $(t^2 + t + 1)/(t+1)$           
                          t_{cap} : 0.0001744            Non-dim. spar cap thickness                      
                          t_{web} : 0.0004826            Non-dim. shear web thickness                     
                                                                                                          
          WingBox, Wing, Aircraft |                                                                       
                          I_{cap} : 1.84e-05             Non-dim spar cap area moment of inertia          
                              M_r : 1.05e+06   [N]       Root moment per root chord                       
                          W_{cap} : 3.667e+04  [N]       Weight of spar caps                              
                          W_{web} : 4142       [N]       Weight of shear web                              
                              \nu : 0.8612               Dummy variable = $(t^2 + t + 1)/(t+1)$           
                          t_{cap} : 0.00411              Non-dim. spar cap thickness                      
                          t_{web} : 0.002063             Non-dim. shear web thickness                     

Constants
---------
                         Aircraft |                                                                                
                            Range : 3000       [nmi#exact] Range                                                   
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
                N_{lift} : 0.273    Wing loading multiplier                          
                       g : 0.2574   Gravitational acceleration                       
              \rho_{cap} : 0.2313   Density of spar cap material                     
               f_{w,add} : 0.07354  Wing added weight fraction                       
                     r_h : 0.02612  Fractional wing thickness at spar web            
              \rho_{web} : 0.02612  Density of shear web material                    
                       w : -0.01564 Wingbox-width-to-chord ratio                     
      \sigma_{max,shear} : -0.02612 Allowable shear stress                           
            \sigma_{max} : -0.2469  Allowable tensile stress                         
                                                                                     
          Wing, Aircraft |                                                           
             C_{D_{0_w}} : 0.826    Wing parasitic drag coefficient                  
            C_{L_{wmax}} : 0.273    Lift coefficient (wing)                          
           \tan(\Lambda) : 0.2045   tangent of wing sweep                            
                  \eta_w : -0.8179  Lift efficiency (diff b/w sectional, actual lift)
          \alpha_{max,w} : -1.084   Max angle of attack                              
                                                                                     
  VerticalTail, Aircraft |                                                           
                     T_e : 0.07594  Thrust per engine at takeoff                     
                     V_c : 0.0398   Cruise velocity                                  
            C_{L_{vmax}} : 0.03795  Max lift coefficient                             
                  \rho_c : 0.01985  Air density (35,000ft)                           
                       e : -0.0201  Span efficiency of vertical tail                 
               \rho_{TO} : -0.03799 Air density (SL))                                
              c_{l_{vt}} : -0.05767 Sectional lift force coefficient (engine out)    
                     V_1 : -0.1519  Minimum takeoff velocity                         
                                                                                     
   LandingGear, Aircraft |                                                           
              W_{0_{lg}} : 0.1797   Weight of aircraft excluding landing gear        
               f_{add,m} : 0.03212  Proportional added weight, main                  
                       K : 0.02604  Column effective length factor                   
       \tan(\theta_{TO}) : 0.02594  Takeoff pitch angle                              
               \rho_{st} : 0.01381  Density of 4340 Steel                            
                       g : 0.01381  Gravitational acceleration                       
                       E : -0.01302 Modulus of elasticity, 4340 steel                
                  n_{mg} : -0.09561 Number of main gear struts                       
                 n_{wps} : -0.1062  Number of wheels per strut                       
                                                                                     
HorizontalTail, Aircraft |                                                           
              \Delta x_w : 0.01677  Distance from aerodynamic centre to CG           
                                                                                     
      Fuselage, Aircraft |                                                           
                n_{seat} : 0.6765    Number of seats                                 
                      LF : 0.3176   Load factor                                      
                f_{padd} : 0.2912   Other misc weight as fraction of payload weight  
           W_{avg. pass} : 0.2802   Average passenger weight                         
                \Delta h : 0.275    Distance from floor to widest part of fuselage   
                     p_s : 0.2295   Seat pitch                                       
           \rho_{\infty} : 0.1881   Air density (35,000ft)                           
                       g : 0.1475   Gravitational acceleration                       
               W'_{seat} : 0.1295   Weight per seat                                  
             \rho_{skin} : 0.08197  Skin density                                     
                \Delta p : 0.08197  Pressure difference across fuselage skin         
                 W_{fix} : 0.06011  Fixed weights (pilots, cockpit seats, navcom)    
             W'_{window} : 0.04783  Weight/length density of windows                 
             W_{checked} : 0.03736  Ave. checked bag weight                          
             \rho_{cone} : 0.03197  Cone material density                            
             W''_{floor} : 0.02629  Floor weight/area density                        
                N_{land} : 0.02576  Emergency landing load factor                    
            \rho_{floor} : 0.02576  Floor material density                           
                 f_{apu} : 0.02548  APU weight as fraction of payload weight         
              f_{lugg,1} : 0.02491  Proportion of passengers with one suitcase       
              f_{string} : 0.02215  Fractional weight of stringers                   
               W_{cargo} : 0.02095  Cargo weight                                     
             W''_{insul} : 0.02086  Weight/area density of insulation material       
               f_{frame} : 0.01582  Fractional frame weight                          
               p_{cabin} : 0.01391  Cabin air pressure (8,000ft)                     
                f_{fadd} : 0.01266  Fractional added weight of local reinforcements  
              f_{lugg,2} : 0.01245  Proportion of passengers with two suitcases      
               T_{cabin} : -0.01391 Cabin temperature                                
                       R : -0.01391 Universal gas constant                           
          \sigma_{floor} : -0.02408 Max allowable cap stress                         
           \sigma_{skin} : -0.1139  Max allowable skin stress                        
                     SPR : -0.2295  Number of seats per row                          
                                                                                     
                Aircraft |                                                           
                       g : 1.231    Gravitational acceleration                       
                    TSFC : 1.231    Thrust specific fuel consumption                 
                   Range : 1.231    Range                                            
                  V_{ne} : 0.6332   Never exceed velocity                            
                  \rho_0 : 0.2786   Air density (0 ft)                               
                 y_{eng} : 0.07503  Engine moment arm                                
                 W_{eng} : 0.04505  Engine weight                                    
                     \mu : 0.04426  Dynamic viscosity (35,000ft)                     
                 d_{fan} : 0.01503  Fan diameter                                     
                    \rho : -0.3906  Air density                                      
              V_{\infty} : -1.58    Freestream velocity                              

