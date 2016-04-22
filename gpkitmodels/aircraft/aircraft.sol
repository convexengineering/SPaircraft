Beginning signomial solve.
W_{fuel}_Aircraft has no upper bound
W_{fuel}_Aircraft has no upper bound
W_{fuel}_Aircraft has no upper bound
W_{fuel}_Aircraft has no upper bound
W_{fuel}_Aircraft has no upper bound
W_{fuel}_Aircraft has no upper bound
Solving took 24 GP solves and 11.8 seconds.
Warning: Constraint x_{CG_{vt}}_Aircraft >= 0.5*\Delta x_{lead}_VerticalTail, Aircraft + 0.5*\Delta x_{trail}_VerticalTail, Aircraft + x_{CG}_Aircraft is not tight because x_{CG_{vt}}_Aircraft [m] evaluated to 56.9693877678 meter but 0.5*\Delta x_{lead}_VerticalTail, Aircraft + 0.5*\Delta x_{trail}_VerticalTail, Aircraft + x_{CG}_Aircraft [m] evaluated to 53.5771489161 meter
Warning: Constraint W_{lg}_Aircraft*x_{CG_{lg}}_Aircraft >= W_{mg}_LandingGear, Aircraft*x_m_LandingGear, Aircraft + W_{ng}_LandingGear, Aircraft*x_n_LandingGear, Aircraft is not tight because W_{lg}_Aircraft*x_{CG_{lg}}_Aircraft [N*m] evaluated to 335262.568087 meter * newton but W_{mg}_LandingGear, Aircraft*x_m_LandingGear, Aircraft + W_{ng}_LandingGear, Aircraft*x_n_LandingGear, Aircraft [N*m] evaluated to 328498.716949 meter * newton
Warning: Constraint x_{CG_{ht}}_Aircraft >= 0.5*\Delta x_{{lead}_h}_HorizontalTail, Aircraft + 0.5*\Delta x_{{trail}_h}_HorizontalTail, Aircraft + x_{CG}_Aircraft is not tight because x_{CG_{ht}}_Aircraft [m] evaluated to 57.8789988471 meter but 0.5*\Delta x_{{lead}_h}_HorizontalTail, Aircraft + 0.5*\Delta x_{{trail}_h}_HorizontalTail, Aircraft + x_{CG}_Aircraft [m] evaluated to 56.977938243 meter
Warning: Constraint 0.5*b_w_Wing, Aircraft*c_{root}_Wing, Aircraft + 0.5*b_w_Wing, Aircraft*c_{tip}_Wing, Aircraft >= S_w_Aircraft is not tight because 0.5*b_w_Wing, Aircraft*c_{root}_Wing, Aircraft + 0.5*b_w_Wing, Aircraft*c_{tip}_Wing, Aircraft [m**2] evaluated to 208.18879665 meter ** 2 but S_w_Aircraft [m**2] evaluated to 92.2408030543 meter ** 2

Cost
----
 1.352e+05 [N] 

Free Variables
--------------
                         Aircraft |                                                                                
                              C_D : 0.07276                       Drag coefficient                                 
                          C_{L_w} : 0.4141                        Lift coefficient (wing)                          
                       C_{L_{aw}} : 4.141                         Lift curve slope (wing)                          
                                D : 6.982e+04  [N]                Total aircraft drag (cruise)                     
                         D_{fuse} : 1.237e+04  [N]                Fuselage drag                                    
                           D_{ht} : 589.7      [N]                Horizontal tail drag                             
                           D_{vt} : 1068       [N]                Vertical tail drag                               
                         D_{wing} : 5.579e+04  [N]                Wing drag                                        
                              L_w : 3.974e+05  [N]                Wing lift                                        
                      L_{v_{max}} : 6.829e+05  [N]                Maximum load for structural sizing               
                              S_w : 92.24      [m**2]             Wing reference area                              
                                W : 3.974e+05  [N]                Total aircraft weight                            
                        W_{empty} : 2.622e+05  [N]                Empty aircraft weight                            
                         W_{fuel} : 1.352e+05  [N]                Fuel weight                                      
                         W_{fuse} : 1.657e+05  [N]                Fuselage weight                                  
                           W_{ht} : 1300       [N]                Horizontal tail weight                           
                           W_{lg} : 1.64e+04   [N]                Weight of landing gear                           
                           W_{vt} : 1395       [N]                Vertical tail weight                             
                         W_{wing} : 6.74e+04   [N]                Wing weight                                      
                          Z_{bre} : 4.695e+04  [m**0.5*s/kg**0.5] Breguet parameter                                
                   \bar{c}_{wing} : 9.503      [m]                Mean aerodynamic chord (wing)                    
                           b_{vt} : 3.452      [m]                Vertical tail span                               
                           c_{vt} : 9.215      [m]                Vertical tail root chord                         
                         h_{hold} : 0.9029     [m]                Hold height                                      
                         l_{fuse} : 58.18      [m]                Fuselage length                                  
                    p_{\lambda_v} : 1.6                           1 + 2*Tail taper ratio                           
                         w_{fuse} : 3.962      [m]                Fuselage width                                   
                              x_w : 21.32      [m]                Position of wing aerodynamic center              
                     x_{CG_{eng}} : 21.32      [m]                x-location of engine CG                          
                      x_{CG_{fu}} : 33.02      [m]                x-location of fuselage CG                        
                      x_{CG_{ht}} : 57.88      [m]                Horizontal tail CG location                      
                      x_{CG_{lg}} : 20.44      [m]                x-location of landing gear CG                    
                      x_{CG_{vt}} : 56.97      [m]                x-location of tail CG                            
                    x_{CG_{wing}} : 21.32      [m]                x-location of wing CG                            
                           x_{CG} : 19.32      [m]                x-location of CG                                 
                           x_{up} : 29.61      [m]                Fuselage upsweep point                           
                                                                                                                   
               Fuselage, Aircraft |                                                                                
                        A_{floor} : 0.06427    [m**2]             Floor beam x-sectional area                      
                         A_{fuse} : 12.33      [m**2]             Fuselage x-sectional area                        
                         A_{hold} : 2.166      [m**2]             Cargo hold x-sectional area                      
                         A_{skin} : 0.0124     [m**2]             Skin cross sectional area                        
                     D_{friction} : 1.138e+04  [N]                Friction drag                                    
                      D_{upsweep} : 986.2      [N]                Drag due to fuse upsweep                         
                               FF : 1.056                         Fuselage form factor                             
                        M_{floor} : 4.862e+05  [N*m]              Max bending moment in floor beams                
                        P_{floor} : 1.137e+06  [N]                Distributed floor load                           
                         R_{fuse} : 1.981      [m]                Fuselage radius                                  
                         S_{bulk} : 24.66      [m**2]             Bulkhead surface area                            
                        S_{floor} : 5.686e+05  [N]                Maximum shear in floor beams                     
                         S_{nose} : 53.52      [m**2]             Nose surface area                                
                         V_{bulk} : 0.02456    [m**3]             Bulkhead skin volume                             
                        V_{cabin} : 360.3      [m**3]             Cabin volume                                     
                        V_{cargo} : 6.796      [m**3]             Cargo volume                                     
                         V_{cone} : 0.1553     [m**3]             Cone skin volume                                 
                          V_{cyl} : 0.3026     [m**3]             Cylinder skin volume                             
                        V_{floor} : 0.2198     [m**3]             Floor volume                                     
                         V_{hold} : 52.88      [m**3]             Hold volume                                      
                         V_{lugg} : 18.24      [m**3]             Luggage volume                                   
                         V_{nose} : 0.0533     [m**3]             Nose skin volume                                 
                          W_{apu} : 5657       [N]                APU weight                                       
                         W_{buoy} : 1736       [N]                Buoyancy weight                                  
                         W_{cone} : 7405       [N]                Cone weight                                      
                        W_{floor} : 1.164e+04  [N]                Floor weight                                     
                        W_{insul} : 4622       [N]                Insulation material weight                       
                         W_{lugg} : 1.79e+04   [N]                Passenger luggage weight                         
                         W_{padd} : 6.465e+04  [N]                Misc weights (galley, toilets, doors etc.)       
                         W_{pass} : 1.337e+05  [N]                Passenger weight                                 
                          W_{pay} : 1.616e+05  [N]                Payload weight                                   
                         W_{seat} : 2.79e+04   [N]                Seating weight                                   
                        W_{shell} : 1.814e+04  [N]                Shell weight                                     
                         W_{skin} : 1.008e+04  [N]                Skin weight                                      
                       W_{window} : 1.062e+04  [N]                Window weight                                    
                   \lambda_{cone} : 0.4                           Tailcone radius taper ratio (xshell2->xtail)     
                             \phi : 0.08338                       Upsweep angle                                    
                     \rho_{cabin} : 0.8711     [kg/m**3]          Air density in cabin                             
                         \sigma_x : 3.831e+07  [N/m**2]           Axial stress in skin                             
                  \sigma_{\theta} : 1.034e+08  [N/m**2]           Skin hoop stress                                 
                      \tau_{cone} : 1.034e+08  [N/m**2]           Shear stress in cone                             
                                f : 14.69                         Fineness ratio                                   
                        h_{floor} : 0.07814    [m]                Floor I-beam height                              
                         l_{cone} : 23.04      [m]                Cone length                                      
                        l_{floor} : 28.37      [m]                Floor length                                     
                         l_{nose} : 5.2        [m]                Nose length                                      
                        l_{shell} : 24.41      [m]                Shell length                                     
                         n_{pass} : 167                           Number of passengers                             
                         n_{rows} : 31                            Number of rows                                   
                        t_{shell} : 0.001344   [m]                Shell thickness                                  
                         t_{skin} : 0.0009958  [m]                Skin thickness                                   
                        w_{floor} : 3.42       [m]                Floor width                                      
                           xVbulk : 0.7271     [m**4]             Volume moment of bulkhead                        
                            xVcyl : 8.271      [m**4]             Volume moment of cylinder                        
                           xVnose : 0.1386     [m**4]             Volume moment of nose                            
                            xWapu : 2.069e+05  [N*m]              Moment of APU                                    
                           xWcone : 5.126e+05  [N*m]              Moment of cone                                   
                            xWfix : 2.802e+04  [N*m]              Moment of fixed weights                          
                          xWfloor : 3.57e+05   [N*m]              Moment of floor weight                           
                           xWfuse : 5.471e+06  [N*m]              Fuselage moment                                  
                          xWinsul : 1.93e+05   [N*m]              Moment of insulation material                    
                           xWpadd : 1.359e+06  [N*m]              Moment of misc weights                           
                           xWseat : 6.824e+05  [N*m]              Moment of seats                                  
                          xWshell : 6.691e+05  [N*m]              Mass moment of shell                             
                           xWskin : 3.717e+05  [N*m]              Mass moment of skin                              
                         xWwindow : 3.347e+05  [N*m]              Mass moment of windows                           
                       x_{shell1} : 5.2        [m]                Start of cylinder section                        
                       x_{shell2} : 29.61      [m]                End of cylinder section                          
                                                                                                                   
         HorizontalTail, Aircraft |                                                                                
                             AR_h : 4.569                         Horizontal tail aspect ratio                     
                          C_{D_h} : 0.005915                      Horizontal tail drag coefficient                 
                      C_{D_{0_h}} : 0.005335                      Horizontal tail parasitic drag coefficient       
                          C_{L_h} : 0.09059                       Lift coefficient (htail)                         
                       C_{L_{ah}} : 3.683                         Lift curve slope (htail)                         
                              K_f : 0.4317                        Empirical factor for fuselage-wing interference  
                              L_h : 9032       [N]                Horizontal tail downforce                        
                      L_{{max}_h} : 3.164e+05  [N]                Maximum load                                     
                         Re_{c_h} : 1.056e+07                     Cruise Reynolds number (Horizontal tail)         
                             S.M. : 0.05                          Stability margin                                 
                              S_h : 9.583      [m**2]             Horizontal tail area                             
              \Delta x_{{lead}_h} : 36.46      [m]                Distance from CG to horizontal tail leading edge 
             \Delta x_{{trail}_h} : 38.87      [m]                Distance from CG to horizontal tail trailing edge
                           \alpha : 0.0246                        Horizontal tail angle of attack                  
                     \bar{c}_{ht} : 1.663      [m]                Mean aerodynamic chord (ht)                      
                        \lambda_h : 0.2                           Horizontal tail taper ratio                      
                           \tau_h : 0.15                          Horizontal tail thickness/chord ratio            
                           b_{ht} : 6.617      [m]                Horizontal tail span                             
                       c_{root_h} : 2.414      [m]                Horizontal tail root chord                       
                        c_{tip_h} : 0.4828     [m]                Horizontal tail tip chord                        
                              e_h : 0.9851                        Oswald efficiency factor                         
                     f(\lambda_h) : 0.0033                        Empirical efficiency function of taper           
                              l_h : 37.96      [m]                Horizontal tail moment arm                       
                           p_{ht} : 1.4                           Substituted variable = 1 + 2*taper               
                           q_{ht} : 1.2                           Substituted variable = 1 + taper                 
                 y_{\bar{c}_{ht}} : 1.891      [m]                Vertical location of mean aerodynamic chord      
                                                                                                                   
            LandingGear, Aircraft |                                                                                
                                B : 6.751      [m]                Landing gear base                                
                         E_{land} : 3.809e+05  [J]                Max KE to be absorbed in landing                 
                          F_{w_m} : 7161                          Weight factor (main)                             
                          F_{w_n} : 643.8                         Weight factor (nose)                             
                              I_m : 7.309e-06  [m**4]             Area moment of inertia (main strut)              
                              I_n : 1.287e-06  [m**4]             Area moment of inertia (nose strut)              
                              L_m : 6.437e+05  [N]                Max static load through main gear                
                              L_n : 1.609e+05  [N]                Min static load through nose gear                
                      L_{n_{dyn}} : 1.624e+05  [N]                Dyn. braking load, nose gear                     
                          L_{w_m} : 1.609e+05  [N]                Static load per wheel (main)                     
                          L_{w_n} : 8.044e+04  [N]                Static load per wheel (nose)                     
                             S_sa : 0.2959     [m]                Stroke of the shock absorber                     
                                T : 6.153      [m]                Main landing gear track                          
                           W_{mg} : 1.486e+04  [N]                Weight of main gear                              
                           W_{ms} : 1486       [N]                Weight of main struts                            
                           W_{mw} : 2377       [N]                Weight of main wheels (per strut)                
                           W_{ng} : 1543       [N]                Weight of nose gear                              
                           W_{ns} : 172.4      [N]                Weight of nose strut                             
                           W_{nw} : 548.2      [N]                Weight of nose wheels (total)                    
                         W_{wa,m} : 267.2      [lbf]              Wheel assembly weight for single main gear wheel 
                         W_{wa,n} : 61.62      [lbf]              Wheel assembly weight for single nose gear wheel 
                       \Delta x_m : 1.35       [m]                Distance b/w main gear and CG                    
                       \Delta x_n : 5.402      [m]                Distance b/w nose gear and CG                    
                       \tan(\phi) : 0.2679                        Angle b/w main gear and CG                       
                       \tan(\psi) : 1.963                         Tip over angles                                  
                      d_{nacelle} : 2.05       [m]                Nacelle diameter                                 
                         d_{oleo} : 0.3735     [m]                Diameter of oleo shock absorber                  
                          d_{t_m} : 44.5       [in]               Diameter of main gear tires                      
                          d_{t_n} : 35.6       [in]               Diameter of nose gear tires                      
                              l_m : 2.397      [m]                Length of main gear                              
                              l_n : 1.627      [m]                Length of nose gear                              
                         l_{oleo} : 0.7397     [m]                Length of oleo shock absorber                    
                              r_m : 0.04261    [m]                Radius of main gear struts                       
                              r_n : 0.04325    [m]                Radius of nose gear struts                       
                              t_m : 0.03006    [m]                Thickness of main gear strut wall                
                              t_n : 0.005063   [m]                Thickness of nose gear strut wall                
                          w_{t_m} : 0.4088     [m]                Width of main tires                              
                          w_{t_n} : 0.3271     [m]                Width of nose tires                              
                              x_m : 20.67      [m]                x-location of main gear                          
                              x_n : 13.91      [m]                x-location of nose gear                          
                              y_m : 3.077      [m]                y-location of main gear (symmetric)              
                                                                                                                   
           VerticalTail, Aircraft |                                                                                
                          A_{fan} : 2.405      [m**2]             Engine reference area                            
                           A_{vt} : 0.5763                        Vertical tail aspect ratio                       
                      C_{D_{vis}} : 0.004966                      Viscous drag coefficient                         
                       C_{L_{vt}} : 0.3717                        Vertical tail lift coefficient                   
                           D_{wm} : 3112       [N]                Engine out windmill drag                         
                     L_{max_{vt}} : 1.366e+06  [N]                Maximum load for structural sizing               
                           L_{vt} : 1.989e+04  [N]                Vertical tail lift in engine out                 
                          Re_{vt} : 4.172e+07                     Vertical tail reynolds number, cruise            
                                S : 41.36      [m**2]             Vertical tail reference area (full)              
                           S_{vt} : 20.68      [m**2]             Vertical tail ref. area (half)                   
                       W_{struct} : 2791       [N]                Full span weight                                 
                  \Delta x_{lead} : 29.65      [m]                Distance from CG to vertical tail leading edge   
                 \Delta x_{trail} : 38.87      [m]                Distance from CG to vertical tail trailing edge  
                     \bar{c}_{vt} : 6.569      [m]                Vertical tail mean aero chord                    
                     \lambda_{vt} : 0.3                           Vertical tail taper ratio                        
                        \tau_{vt} : 0.1297                        Vertical tail thickness/chord ratio              
                                b : 6.904      [m]                Span                                             
                    c_{root_{vt}} : 9.215      [m]                Vertical tail root chord                         
                     c_{tip_{vt}} : 2.765      [m]                Vertical tail tip chord                          
                           l_{vt} : 32.08      [m]                Vertical tail moment arm                         
                           p_{vt} : 1.6                           Substituted variable = 1 + 2*taper               
                           q_{vt} : 1.3                           Substituted variable = 1 + taper                 
                 z_{\bar{c}_{vt}} : 0.935      [m]                Vertical location of mean aerodynamic chord      
                                                                                                                   
                   Wing, Aircraft |                                                                                
                               AR : 6.859                         Wing aspect ratio                                
                          C_{D_w} : 0.05814                       Drag coefficient                                 
                      L_{max_{w}} : 2.929e+06  [N]                Maximum load                                     
                         \alpha_w : 0.1                           Wing angle of attack                             
                          \lambda : 0.2                           Wing taper ratio                                 
                           \tau_w : 0.15                          Wing thickness/chord ratio                       
                              b_w : 25.15      [m]                Wing span                                        
                         c_{root} : 13.79      [m]                Wing root chord                                  
                          c_{tip} : 2.759      [m]                Wing tip chord                                   
                              e_w : 0.9779                        Oswald efficiency factor                         
                     f(\lambda_w) : 0.0033                        Empirical efficiency function of taper           
                              p_w : 1.4                           Substituted variable = 1 + 2*taper               
                              q_w : 1.2                           Substituted variable = 1 + taper                 
                      y_{\bar{c}} : 7.186      [m]                Spanwise location of mean aerodynamic chord      
                                                                                                                   
WingBox, HorizontalTail, Aircraft |                                                                                
                          I_{cap} : 8.685e-06                     Non-dim spar cap area moment of inertia          
                              M_r : 8.434e+04  [N]                Root moment per root chord                       
                          W_{cap} : 791.6      [N]                Weight of spar caps                              
                          W_{web} : 137.3      [N]                Weight of shear web                              
                              \nu : 0.8612                        Dummy variable = $(t^2 + t + 1)/(t+1)$           
                          t_{cap} : 0.001875                      Non-dim. spar cap thickness                      
                          t_{web} : 0.001445                      Non-dim. shear web thickness                     
                                                                                                                   
  WingBox, VerticalTail, Aircraft |                                                                                
                                A : 1.153                         Aspect ratio                                     
                          I_{cap} : 6.41e-07                      Non-dim spar cap area moment of inertia          
                              M_r : 1.049e+05  [N]                Root moment per root chord                       
                          W_{cap} : 1300       [N]                Weight of spar caps                              
                          W_{web} : 693        [N]                Weight of shear web                              
                              \nu : 0.8225                        Dummy variable = $(t^2 + t + 1)/(t+1)$           
                          t_{cap} : 0.0001807                     Non-dim. spar cap thickness                      
                          t_{web} : 0.0004951                     Non-dim. shear web thickness                     
                                                                                                                   
          WingBox, Wing, Aircraft |                                                                                
                          I_{cap} : 1.882e-05                     Non-dim spar cap area moment of inertia          
                              M_r : 1.172e+06  [N]                Root moment per root chord                       
                          W_{cap} : 4.332e+04  [N]                Weight of spar caps                              
                          W_{web} : 4830       [N]                Weight of shear web                              
                              \nu : 0.8612                        Dummy variable = $(t^2 + t + 1)/(t+1)$           
                          t_{cap} : 0.00421                       Non-dim. spar cap thickness                      
                          t_{web} : 0.002087                      Non-dim. shear web thickness                     

Constants
---------
                         Aircraft |                                                                                
                            Range : 3000       [nmi#exact] Range                                                   
                       V_{\infty} : 234        [m/s]       Freestream velocity                                     
                           V_{ne} : 144        [m/s]       Never exceed velocity                                   
                          W_{eng} : 1e+04      [N]         Engine weight                                           
                              \mu : 1.4e-05    [N*s/m**2]  Dynamic viscosity (35,000ft)                            
                             \rho : 0.38       [kg/m**3]   Air density                                             
                           \rho_0 : 1.225      [kg/m**3]   Air density (0 ft)                                      
                              c_t : 9e-05      [1/s]       Thrust specific fuel consumption                        
                          d_{fan} : 1.75       [m]         Fan diameter                                            
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
                N_{lift} : 0.3226   Wing loading multiplier                                 
                       g : 0.3036   Gravitational acceleration                              
              \rho_{cap} : 0.2732   Density of spar cap material                            
               f_{w,add} : 0.08675  Wing added weight fraction                              
                     r_h : 0.03046  Fractional wing thickness at spar web                   
              \rho_{web} : 0.03046  Density of shear web material                           
                       w : -0.01898 Wingbox-width-to-chord ratio                            
      \sigma_{max,shear} : -0.03046 Allowable shear stress                                  
            \sigma_{max} : -0.2921  Allowable tensile stress                                
                                                                                            
          Wing, Aircraft |                                                                  
             C_{D_{0_w}} : 0.9911   Wing parasitic drag coefficient                         
            C_{L_{wmax}} : 0.3226   Lift coefficient (wing)                                 
           \tan(\Lambda) : 0.2444   tangent of wing sweep                                   
                  \eta_w : -0.9775  Lift efficiency (diff b/w sectional, actual lift)       
          \alpha_{max,w} : -1.296   Max angle of attack                                     
                                                                                            
  VerticalTail, Aircraft |                                                                  
                     T_e : 0.0804   Thrust per engine at takeoff                            
                     V_c : 0.04404  Cruise velocity                                         
            C_{L_{vmax}} : 0.03965  Max lift coefficient                                    
                  \rho_c : 0.02198  Air density (35,000ft)                                  
                       e : -0.02113 Span efficiency of vertical tail                        
               \rho_{TO} : -0.04074 Air density (SL))                                       
              c_{l_{vt}} : -0.06121 Sectional lift force coefficient (engine out)           
                     V_1 : -0.1608  Minimum takeoff velocity                                
                                                                                            
   LandingGear, Aircraft |                                                                  
              W_{0_{lg}} : 0.1853   Weight of aircraft excluding landing gear               
               f_{add,m} : 0.03213  Proportional added weight, main                         
       \tan(\theta_{TO}) : 0.02878  Takeoff pitch angle                                     
                       K : 0.02677  Column effective length factor                          
               \rho_{st} : 0.01416  Density of 4340 Steel                                   
                       g : 0.01416  Gravitational acceleration                              
                       E : -0.01338 Modulus of elasticity, 4340 steel                       
                  n_{mg} : -0.1004  Number of main gear struts                              
                 n_{wps} : -0.1114  Number of wheels per strut                              
                                                                                            
HorizontalTail, Aircraft |                                                                  
              \Delta x_w : 0.01622  Distance from aerodynamic centre to CG                  
                  \eta_h : -0.01217 Lift efficiency (diff between sectional and actual lift)
                                                                                            
      Fuselage, Aircraft |                                                                  
                n_{seat} : 0.6845    Number of seats                                        
                      LF : 0.318    Load factor                                             
                f_{padd} : 0.2912   Other misc weight as fraction of payload weight         
                \Delta h : 0.2878   Distance from floor to widest part of fuselage          
           W_{avg. pass} : 0.2806   Average passenger weight                                
                     p_s : 0.237    Seat pitch                                              
           \rho_{\infty} : 0.2024   Air density (35,000ft)                                  
                       g : 0.1491   Gravitational acceleration                              
               W'_{seat} : 0.1295   Weight per seat                                         
             \rho_{skin} : 0.0817   Skin density                                            
                \Delta p : 0.0817   Pressure difference across fuselage skin                
                 W_{fix} : 0.06011  Fixed weights (pilots, cockpit seats, navcom)           
             W'_{window} : 0.04783  Weight/length density of windows                        
             W_{checked} : 0.03741  Ave. checked bag weight                                 
             \rho_{cone} : 0.03336  Cone material density                                   
             W''_{floor} : 0.02623  Floor weight/area density                               
                N_{land} : 0.02623  Emergency landing load factor                           
            \rho_{floor} : 0.02623  Floor material density                                  
                 f_{apu} : 0.02548  APU weight as fraction of payload weight                
              f_{lugg,1} : 0.02494  Proportion of passengers with one suitcase              
              f_{string} : 0.02237  Fractional weight of stringers                          
               W_{cargo} : 0.02098  Cargo weight                                            
             W''_{insul} : 0.02082  Weight/area density of insulation material              
               f_{frame} : 0.01598  Fractional frame weight                                 
               p_{cabin} : 0.01387  Cabin air pressure (8,000ft)                            
                f_{fadd} : 0.01278  Fractional added weight of local reinforcements         
              f_{lugg,2} : 0.01247  Proportion of passengers with two suitcases             
               T_{cabin} : -0.01387 Cabin temperature                                       
                       R : -0.01387 Universal gas constant                                  
          \sigma_{floor} : -0.02454 Max allowable cap stress                                
           \sigma_{skin} : -0.1151  Max allowable skin stress                               
                     SPR : -0.237   Number of seats per row                                 
                                                                                            
                Aircraft |                                                                  
                     c_t : 1.442    Thrust specific fuel consumption                        
                   Range : 1.442    Range                                                   
                  V_{ne} : 0.7365   Never exceed velocity                                   
                  \rho_0 : 0.3286   Air density (0 ft)                                      
                 y_{eng} : 0.07993  Engine moment arm                                       
                     \mu : 0.0476   Dynamic viscosity (35,000ft)                            
                 W_{eng} : 0.04505  Engine weight                                           
                 d_{fan} : 0.01385  Fan diameter                                            
                    \rho : -0.4593  Air density                                             
              V_{\infty} : -1.896   Freestream velocity                                     

