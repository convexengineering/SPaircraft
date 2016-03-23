Beginning signomial solve.
Solving took 4 GP solves and 0.607 seconds.

Cost
----
 2.546e+04 [N] 

Free Variables
--------------
            AR_h : 4.748             Horizontal tail aspect ratio                     
         C_{D_h} : 0.0387            Horizontal tail drag coefficient                 
         C_{L_h} : 0.3222            Lift coefficient (htail)                         
      C_{L_{ah}} : 3.731             Lift curve slope (htail)                         
             D_h : 1.752e+04  [N]    Horizontal tail drag                             
             K_f : 0.7831            Empirical factor for fuselage-wing interference  
            S.M. : 0.05              Stability margin                                 
             S_h : 41.36      [m**2] Horizontal tail area                             
               W : 1.589e+04  [N]    Horizontal tail weight                           
 \Delta x_{lead} : 15.93      [m]    Distance from CG to horizontal tail leading edge 
\Delta x_{trail} : 20         [m]    Distance from CG to horizontal tail trailing edge
          \alpha : 0.08636           Horizontal tail angle of attack                  
    \bar{c}_{ht} : 3.093      [m]    Mean aerodynamic chord (ht)                      
          b_{ht} : 14.01      [m]    Horizontal tail span                             
        c_{root} : 4.071      [m]    Horizontal tail root chord                       
         c_{tip} : 1.832      [m]    Horizontal tail tip chord                        
             l_h : 18.76      [m]    Horizontal tail moment arm                       
               p : 1.9               Substituted variable = 1 + 2*taper               
               q : 1.45              Substituted variable = 1 + taper                 
             x_w : 22         [m]    Position of wing aerodynamic center              
     y_{\bar{c}} : 3.565      [m]    Vertical location of mean aerodynamic chord      
                                                                                      
         WingBox |                                                                    
         I_{cap} : 1.361e-05         Non-dim spar cap area moment of inertia          
             M_r : 3.759e+05  [N]    Root moment per root chord                       
         W_{cap} : 1.012e+04  [N]    Weight of spar caps                              
         W_{web} : 1224       [N]    Weight of shear web                              
             \nu : 0.786             Dummy variable = $(t^2 + t + 1)/(t+1)$           
            \tau : 0.15              Thickness to chord ratio                         
         t_{cap} : 0.002988          Non-dim. spar cap thickness                      
         t_{web} : 0.001606          Non-dim. shear web thickness                     

Constants
---------
       C_{D_{0_h}} : 0.03                Horizontal tail parasitic drag coefficient              
           C_{L_w} : 0.5                 Lift coefficient (wing)                                 
        C_{L_{aw}} : 6.283               Lift curve slope (wing)                                 
        C_{m_{ac}} : 0.1                 Moment coefficient about aerodynamic centre (wing)      
      C_{m_{fuse}} : 0.1                 Moment coefficient (fuselage)                           
        S.M._{min} : 0.05                Minimum stability margin                                
               S_w : 125       [m**2]    Wing area                                               
        V_{\infty} : 240       [m/s]     Freestream velocity                                     
        \Delta x_w : 2         [m]       Distance from aerodynamic centre to CG                  
      \alpha_{max} : 0.1                 Max angle of attack (htail)                             
    \bar{c}_{wing} : 5         [m]       Mean aerodynamic chord (wing)                           
              \eta : 0.97                Lift efficiency (diff between sectional and actual lift)
           \lambda : 0.45                Horizontal tail taper ratio                             
              \rho : 0.38      [kg/m**3] Air density (35,000 ft)                                 
   \tan(\Lambda_h) : 0.5774              tangent of horizontal tail sweep                        
               e_h : 0.8                 Oswald efficiency factor                                
          l_{fuse} : 40        [m]       Fuselage length                                         
               w_f : 6         [m]       Fuselage width                                          
            x_{CG} : 20        [m]       CG location                                             
                                                                                                 
           WingBox |                                                                             
           L_{max} : 1e+06     [N]       Maximum wing load                                       
          N_{lift} : 2                   Wing loading multiplier                                 
        \rho_{cap} : 2700      [kg/m**3] Density of spar cap material                            
        \rho_{web} : 2700      [kg/m**3] Density of shear web material                           
\sigma_{max,shear} : 1.67e+08  [Pa]      Allowable shear stress                                  
      \sigma_{max} : 2.5e+08   [Pa]      Allowable tensile stress                                
         f_{w,add} : 0.4                 Wing added weight fraction                              
                 g : 9.81      [m/s**2]  Gravitational acceleration                              
               r_h : 0.75                Fractional wing thickness at spar web                   
                 w : 0.5                 Wingbox-width-to-chord ratio                            

Sensitivities
-------------
           WingBox |                                                                  
          N_{lift} : 0.3252   Wing loading multiplier                                 
           L_{max} : 0.3252   Maximum wing load                                       
                 g : 0.312    Gravitational acceleration                              
        \rho_{cap} : 0.2783   Density of spar cap material                            
         f_{w,add} : 0.08914  Wing added weight fraction                              
               r_h : 0.03366  Fractional wing thickness at spar web                   
        \rho_{web} : 0.03366  Density of shear web material                           
                 w : -0.0132  Wingbox-width-to-chord ratio                            
\sigma_{max,shear} : -0.03366 Allowable shear stress                                  
      \sigma_{max} : -0.2915  Allowable tensile stress                                
                                                                                      
        V_{\infty} : 1.376    Freestream velocity                                     
            x_{CG} : 1.215    CG location                                             
              \rho : 0.688    Air density (35,000 ft)                                 
               S_w : 0.645    Wing area                                               
       C_{D_{0_h}} : 0.5333   Horizontal tail parasitic drag coefficient              
        \Delta x_w : 0.483    Distance from aerodynamic centre to CG                  
               w_f : 0.4284   Fuselage width                                          
        C_{L_{aw}} : 0.3356   Lift curve slope (wing)                                 
           \lambda : 0.269    Horizontal tail taper ratio                             
    \bar{c}_{wing} : 0.192    Mean aerodynamic chord (wing)                           
           C_{L_w} : 0.1547   Lift coefficient (wing)                                 
      C_{m_{fuse}} : 0.07734  Moment coefficient (fuselage)                           
        C_{m_{ac}} : 0.07734  Moment coefficient about aerodynamic centre (wing)      
        S.M._{min} : 0.03729  Minimum stability margin                                
               e_h : -0.1547  Oswald efficiency factor                                
              \eta : -0.3664  Lift efficiency (diff between sectional and actual lift)
          l_{fuse} : -1.031   Fuselage length                                         

