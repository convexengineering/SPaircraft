Beginning signomial solve.
Solving took 4 GP solves and 0.567 seconds.

Cost
----
 2.489e+04 [N] 

Free Variables
--------------
            AR_h : 4.735             Horizontal tail aspect ratio                     
         C_{D_h} : 0.03948           Horizontal tail drag coefficient                 
         C_{L_h} : 0.3358            Lift coefficient (htail)                         
      C_{L_{ah}} : 3.728             Lift curve slope (htail)                         
             D_h : 1.713e+04  [N]    Horizontal tail drag                             
            S.M. : 0.05              Stability margin                                 
             S_h : 39.64      [m**2] Horizontal tail area                             
               W : 1.552e+04  [N]    Horizontal tail weight                           
 \Delta x_{lead} : 16.01      [m]    Distance from CG to horizontal tail leading edge 
\Delta x_{trail} : 20         [m]    Distance from CG to horizontal tail trailing edge
          \alpha : 0.09009           Horizontal tail angle of attack                  
    \bar{c}_{ht} : 3.032      [m]    Mean aerodynamic chord (ht)                      
          b_{ht} : 13.7       [m]    Horizontal tail span                             
        c_{root} : 3.991      [m]    Horizontal tail root chord                       
         c_{tip} : 1.796      [m]    Horizontal tail tip chord                        
             l_h : 18.78      [m]    Horizontal tail moment arm                       
               p : 1.9               Substituted variable = 1 + 2*taper               
               q : 1.45              Substituted variable = 1 + taper                 
     y_{\bar{c}} : 3.485      [m]    Vertical location of mean aerodynamic chord      
                                                                                      
         WingBox |                                                                    
         I_{cap} : 1.412e-05         Non-dim spar cap area moment of inertia          
             M_r : 3.749e+05  [N]    Root moment per root chord                       
         W_{cap} : 9888       [N]    Weight of spar caps                              
         W_{web} : 1197       [N]    Weight of shear web                              
             \nu : 0.786             Dummy variable = $(t^2 + t + 1)/(t+1)$           
            \tau : 0.15              Thickness to chord ratio                         
         t_{cap} : 0.003106          Non-dim. spar cap thickness                      
         t_{web} : 0.001671          Non-dim. shear web thickness                     

Constants
---------
       C_{D_{0_h}} : 0.03                Horizontal tail parasitic drag coefficient              
           C_{L_w} : 0.5                 Lift coefficient (wing)                                 
        C_{L_{aw}} : 6.283               Lift curve slope (wing)                                 
        C_{m_{ac}} : 0.1                 Moment coefficient about aerodynamic centre (wing)      
      C_{m_{fuse}} : 0.1                 Moment coefficient (fuselage)                           
               K_f : 0.7                 Empirical factor for fuselage-wing interference         
        S.M._{min} : 0.05                Minimum stability margin                                
               S_w : 125       [m**2]    Wing area                                               
        V_{\infty} : 240       [m/s]     Freestream velocity                                     
      \alpha_{max} : 0.1                 Max angle of attack (htail)                             
    \bar{c}_{wing} : 5         [m]       Mean aerodynamic chord (wing)                           
              \eta : 0.97                Lift efficiency (diff between sectional and actual lift)
           \lambda : 0.45                Horizontal tail taper ratio                             
              \rho : 0.38      [kg/m**3] Air density (35,000 ft)                                 
   \tan(\Lambda_h) : 0.5774              tangent of horizontal tail sweep                        
               e_h : 0.8                 Oswald efficiency factor                                
          l_{fuse} : 40        [m]       Fuselage length                                         
               w_f : 6         [m]       Fuselage width                                          
               x_w : 2         [m]       Distance from aerodynamic centre to CG                  
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
          N_{lift} : 0.3256   Wing loading multiplier                                 
           L_{max} : 0.3256   Maximum wing load                                       
                 g : 0.3118   Gravitational acceleration                              
        \rho_{cap} : 0.2781   Density of spar cap material                            
         f_{w,add} : 0.08909  Wing added weight fraction                              
               r_h : 0.03367  Fractional wing thickness at spar web                   
        \rho_{web} : 0.03367  Density of shear web material                           
                 w : -0.01376 Wingbox-width-to-chord ratio                            
\sigma_{max,shear} : -0.03367 Allowable shear stress                                  
      \sigma_{max} : -0.2919  Allowable tensile stress                                
                                                                                      
        V_{\infty} : 1.376    Freestream velocity                                     
            x_{CG} : 0.914    CG location                                             
              \rho : 0.6882   Air density (35,000 ft)                                 
               S_w : 0.6665   Wing area                                               
       C_{D_{0_h}} : 0.523    Horizontal tail parasitic drag coefficient              
               x_w : 0.464    Distance from aerodynamic centre to CG                  
               w_f : 0.3834   Fuselage width                                          
        C_{L_{aw}} : 0.3361   Lift curve slope (wing)                                 
           \lambda : 0.2702   Horizontal tail taper ratio                             
    \bar{c}_{wing} : 0.2025   Mean aerodynamic chord (wing)                           
               K_f : 0.1917   Empirical factor for fuselage-wing interference         
           C_{L_w} : 0.1652   Lift coefficient (wing)                                 
      C_{m_{fuse}} : 0.0826   Moment coefficient (fuselage)                           
        C_{m_{ac}} : 0.0826   Moment coefficient about aerodynamic centre (wing)      
        S.M._{min} : 0.03735  Minimum stability margin                                
               e_h : -0.1652  Oswald efficiency factor                                
              \eta : -0.3513  Lift efficiency (diff between sectional and actual lift)
          l_{fuse} : -0.7223  Fuselage length                                         

