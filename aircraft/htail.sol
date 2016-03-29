Beginning signomial solve.
Solving took 4 GP solves and 0.783 seconds.

Cost
----
 6778 [N] 

Free Variables
--------------
            AR_h : 3.259             Horizontal tail aspect ratio                     
         C_{D_h} : 0.006232          Horizontal tail drag coefficient                 
     C_{D_{0_h}} : 0.005161          Horizontal tail parasitic drag coefficient       
         C_{L_h} : 0.1042            Lift coefficient (htail)                         
      C_{L_{ah}} : 3.217             Lift curve slope (htail)                         
             D_h : 3705       [N]    Horizontal tail drag                             
             K_f : 0.7831            Empirical factor for fuselage-wing interference  
            Re_c : 3.053e+07         Cruise Reynolds number (Horizontal tail)         
            S.M. : 0.05              Stability margin                                 
             S_h : 54.32      [m**2] Horizontal tail area                             
               W : 6146       [N]    Horizontal tail weight                           
 \Delta x_{lead} : 13.2       [m]    Distance from CG to horizontal tail leading edge 
\Delta x_{trail} : 20         [m]    Distance from CG to horizontal tail trailing edge
          \alpha : 0.03239           Horizontal tail angle of attack                  
    \bar{c}_{ht} : 4.687      [m]    Mean aerodynamic chord (ht)                      
         \lambda : 0.2               Horizontal tail taper ratio                      
            \tau : 0.15              Horizontal tail thickness/chord ratio            
          b_{ht} : 13.31      [m]    Horizontal tail span                             
        c_{root} : 6.804      [m]    Horizontal tail root chord                       
         c_{tip} : 1.361      [m]    Horizontal tail tip chord                        
             e_h : 0.9893            Oswald efficiency factor                         
    f\(\lambda\) : 0.003309          Empirical efficiency function of taper           
             l_h : 16.56      [m]    Horizontal tail moment arm                       
               p : 1.4               Substituted variable = 1 + 2*taper               
               q : 1.2               Substituted variable = 1 + taper                 
             x_w : 22         [m]    Position of wing aerodynamic center              
     y_{\bar{c}} : 3.802      [m]    Vertical location of mean aerodynamic chord      
                                                                                      
         WingBox |                                                                    
         I_{cap} : 2.464e-06         Non-dim spar cap area moment of inertia          
             M_r : 1.901e+05  [N]    Root moment per root chord                       
         W_{cap} : 3518       [N]    Weight of spar caps                              
         W_{web} : 872.5      [N]    Weight of shear web                              
             \nu : 0.8612            Dummy variable = $(t^2 + t + 1)/(t+1)$           
         t_{cap} : 0.0005215         Non-dim. spar cap thickness                      
         t_{web} : 0.0005749         Non-dim. shear web thickness                     

Constants
---------
           C_{L_w} : 0.5                  Lift coefficient (wing)                                 
        C_{L_{aw}} : 6.283                Lift curve slope (wing)                                 
      C_{m_{fuse}} : 0.05                 Moment coefficient (fuselage)                           
        S.M._{min} : 0.05                 Minimum stability margin                                
               S_w : 125       [m**2]     Wing area                                               
        V_{\infty} : 240       [m/s]      Freestream velocity                                     
        \Delta x_w : 2         [m]        Distance from aerodynamic centre to CG                  
      \alpha_{max} : 0.1                  Max angle of attack (htail)                             
    \bar{c}_{wing} : 5         [m]        Mean aerodynamic chord (wing)                           
              \eta : 0.97                 Lift efficiency (diff between sectional and actual lift)
               \mu : 1.4e-05   [N*s/m**2] Dynamic viscosity (35,000ft)                            
              \rho : 0.38      [kg/m**3]  Air density (35,000 ft)                                 
   \tan(\Lambda_h) : 0.5774               tangent of horizontal tail sweep                        
          l_{fuse} : 40        [m]        Fuselage length                                         
               w_f : 6         [m]        Fuselage width                                          
            x_{CG} : 20        [m]        CG location                                             
      |C_{m_{ac}}| : 0.1                  Moment coefficient about aerodynamic centre (wing)      
                                                                                                  
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
          N_{lift} : 0.4562   Wing loading multiplier                                 
           L_{max} : 0.4562   Maximum wing load                                       
                 g : 0.4534   Gravitational acceleration                              
        \rho_{cap} : 0.3633   Density of spar cap material                            
         f_{w,add} : 0.1295   Wing added weight fraction                              
               r_h : 0.0901   Fractional wing thickness at spar web                   
        \rho_{web} : 0.0901   Density of shear web material                           
\sigma_{max,shear} : -0.0901  Allowable shear stress                                  
      \sigma_{max} : -0.3661  Allowable tensile stress                                
                                                                                      
            x_{CG} : 1.395    CG location                                             
        V_{\infty} : 1.084    Freestream velocity                                     
              \rho : 0.5373   Air density (35,000 ft)                                 
               w_f : 0.5193   Fuselage width                                          
               S_w : 0.4068   Wing area                                               
        C_{L_{aw}} : 0.4068   Lift curve slope (wing)                                 
        \Delta x_w : 0.3979   Distance from aerodynamic centre to CG                  
    \bar{c}_{wing} : 0.04521  Mean aerodynamic chord (wing)                           
        S.M._{min} : 0.04521  Minimum stability margin                                
   \tan(\Lambda_h) : -0.02297 tangent of horizontal tail sweep                        
              \eta : -0.3611  Lift efficiency (diff between sectional and actual lift)
          l_{fuse} : -1.171   Fuselage length                                         

