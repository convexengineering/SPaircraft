Beginning signomial solve.
Solving took 4 GP solves and 0.634 seconds.

Cost
----
 1.338e+04 [N] 

Free Variables
--------------
            AR_h : 3.421             Horizontal tail aspect ratio                     
         C_{D_h} : 0.01456           Horizontal tail drag coefficient                 
     C_{D_{0_h}} : 0.005179          Horizontal tail parasitic drag coefficient       
         C_{L_h} : 0.284             Lift coefficient (htail)                         
      C_{L_{ah}} : 3.288             Lift curve slope (htail)                         
             D_h : 7946       [N]    Horizontal tail drag                             
             K_f : 0.7831            Empirical factor for fuselage-wing interference  
            Re_c : 2.607e+07         Cruise Reynolds number (Horizontal tail          
            S.M. : 0.05              Stability margin                                 
             S_h : 49.88      [m**2] Horizontal tail area                             
               W : 1.087e+04  [N]    Horizontal tail weight                           
 \Delta x_{lead} : 14.73      [m]    Distance from CG to horizontal tail leading edge 
\Delta x_{trail} : 20         [m]    Distance from CG to horizontal tail trailing edge
          \alpha : 0.08636           Horizontal tail angle of attack                  
    \bar{c}_{ht} : 4.001      [m]    Mean aerodynamic chord (ht)                      
            \tau : 0.15              Horizontal tail thickness/chord ratio            
          b_{ht} : 13.06      [m]    Horizontal tail span                             
        c_{root} : 5.267      [m]    Horizontal tail root chord                       
         c_{tip} : 2.37       [m]    Horizontal tail tip chord                        
             l_h : 17.65      [m]    Horizontal tail moment arm                       
               p : 1.9               Substituted variable = 1 + 2*taper               
               q : 1.45              Substituted variable = 1 + taper                 
             x_w : 22         [m]    Position of wing aerodynamic center              
     y_{\bar{c}} : 3.323      [m]    Vertical location of mean aerodynamic chord      
                                                                                      
         WingBox |                                                                    
         I_{cap} : 5.858e-06         Non-dim spar cap area moment of inertia          
             M_r : 2.708e+05  [N]    Root moment per root chord                       
         W_{cap} : 6625       [N]    Weight of spar caps                              
         W_{web} : 1141       [N]    Weight of shear web                              
             \nu : 0.786             Dummy variable = $(t^2 + t + 1)/(t+1)$           
         t_{cap} : 0.001253          Non-dim. spar cap thickness                      
         t_{web} : 0.0009595         Non-dim. shear web thickness                     

Constants
---------
           C_{L_w} : 0.5                  Lift coefficient (wing)                                 
        C_{L_{aw}} : 6.283                Lift curve slope (wing)                                 
        C_{m_{ac}} : 0.1                  Moment coefficient about aerodynamic centre (wing)      
      C_{m_{fuse}} : 0.1                  Moment coefficient (fuselage)                           
        S.M._{min} : 0.05                 Minimum stability margin                                
               S_w : 125       [m**2]     Wing area                                               
        V_{\infty} : 240       [m/s]      Freestream velocity                                     
        \Delta x_w : 2         [m]        Distance from aerodynamic centre to CG                  
      \alpha_{max} : 0.1                  Max angle of attack (htail)                             
    \bar{c}_{wing} : 5         [m]        Mean aerodynamic chord (wing)                           
              \eta : 0.97                 Lift efficiency (diff between sectional and actual lift)
           \lambda : 0.45                 Horizontal tail taper ratio                             
               \mu : 1.4e-05   [N*s/m**2] Dynamic viscosity (35,000ft)                            
              \rho : 0.38      [kg/m**3]  Air density (35,000 ft)                                 
   \tan(\Lambda_h) : 0.5774               tangent of horizontal tail sweep                        
               e_h : 0.8                  Oswald efficiency factor                                
          l_{fuse} : 40        [m]        Fuselage length                                         
               w_f : 6         [m]        Fuselage width                                          
            x_{CG} : 20        [m]        CG location                                             
                                                                                                  
           WingBox |                                                                              
           L_{max} : 1e+06     [N]        Maximum wing load                                       
          N_{lift} : 2                    Wing loading multiplier                                 
        \rho_{cap} : 2700      [kg/m**3]  Density of spar cap material                            
        \rho_{web} : 2700      [kg/m**3]  Density of shear web material                           
\sigma_{max,shear} : 1.67e+08  [Pa]       Allowable shear stress                                  
      \sigma_{max} : 2.5e+08   [Pa]       Allowable tensile stress                                
         f_{w,add} : 0.4                  Wing added weight fraction                              
                 g : 9.81      [m/s**2]   Gravitational acceleration                              
               r_h : 0.75                 Fractional wing thickness at spar web                   
                 w : 0.5                  Wingbox-width-to-chord ratio                            

Sensitivities
-------------
           WingBox |                                                                  
          N_{lift} : 0.4128   Wing loading multiplier                                 
           L_{max} : 0.4128   Maximum wing load                                       
                 g : 0.4062   Gravitational acceleration                              
        \rho_{cap} : 0.3465   Density of spar cap material                            
         f_{w,add} : 0.1161   Wing added weight fraction                              
               r_h : 0.05969  Fractional wing thickness at spar web                   
        \rho_{web} : 0.05969  Density of shear web material                           
\sigma_{max,shear} : -0.05969 Allowable shear stress                                  
      \sigma_{max} : -0.3531  Allowable tensile stress                                
                                                                                      
        V_{\infty} : 1.183    Freestream velocity                                     
            x_{CG} : 0.9992   CG location                                             
               S_w : 0.8133   Wing area                                               
              \rho : 0.5888   Air density (35,000 ft)                                 
        \Delta x_w : 0.4297   Distance from aerodynamic centre to CG                  
    \bar{c}_{wing} : 0.3879   Mean aerodynamic chord (wing)                           
           C_{L_w} : 0.3826   Lift coefficient (wing)                                 
           \lambda : 0.3239   Horizontal tail taper ratio                             
      C_{m_{fuse}} : 0.1913   Moment coefficient (fuselage)                           
        C_{m_{ac}} : 0.1913   Moment coefficient about aerodynamic centre (wing)      
               w_f : 0.06147  Fuselage width                                          
        C_{L_{aw}} : 0.04816  Lift curve slope (wing)                                 
              \eta : -0.04412 Lift efficiency (diff between sectional and actual lift)
   \tan(\Lambda_h) : -0.0807  tangent of horizontal tail sweep                        
               e_h : -0.3826  Oswald efficiency factor                                
          l_{fuse} : -0.9728  Fuselage length                                         

