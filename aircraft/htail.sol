Beginning signomial solve.
Solving took 5 GP solves and 0.822 seconds.

Cost
----
 1.021e+04 [N] 

Free Variables
--------------
            AR_h : 4.333             Horizontal tail aspect ratio                     
         C_{D_h} : 0.01245           Horizontal tail drag coefficient                 
     C_{D_{0_h}} : 0.005189          Horizontal tail parasitic drag coefficient       
         C_{L_h} : 0.3122            Lift coefficient (htail)                         
      C_{L_{ah}} : 3.614             Lift curve slope (htail)                         
             D_h : 6109       [N]    Horizontal tail drag                             
             K_f : 0.7831            Empirical factor for fuselage-wing interference  
            Re_c : 2.406e+07         Cruise Reynolds number (Horizontal tail)         
            S.M. : 0.05              Stability margin                                 
             S_h : 44.84      [m**2] Horizontal tail area                             
               W : 8199       [N]    Horizontal tail weight                           
 \Delta x_{lead} : 14.64      [m]    Distance from CG to horizontal tail leading edge 
\Delta x_{trail} : 20         [m]    Distance from CG to horizontal tail trailing edge
          \alpha : 0.08636           Horizontal tail angle of attack                  
    \bar{c}_{ht} : 3.693      [m]    Mean aerodynamic chord (ht)                      
         \lambda : 0.2               Horizontal tail taper ratio                      
            \tau : 0.15              Horizontal tail thickness/chord ratio            
          b_{ht} : 13.94      [m]    Horizontal tail span                             
        c_{root} : 5.361      [m]    Horizontal tail root chord                       
         c_{tip} : 1.072      [m]    Horizontal tail tip chord                        
             e_h : 0.9859            Oswald efficiency factor                         
    f\(\lambda\) : 0.0033            Empirical efficiency function of taper           
             l_h : 17.86      [m]    Horizontal tail moment arm                       
               p : 1.4               Substituted variable = 1 + 2*taper               
               q : 1.2               Substituted variable = 1 + taper                 
             x_w : 22         [m]    Position of wing aerodynamic center              
     y_{\bar{c}} : 3.983      [m]    Vertical location of mean aerodynamic chord      
                                                                                      
         WingBox |                                                                    
         I_{cap} : 5.276e-06         Non-dim spar cap area moment of inertia          
             M_r : 2.528e+05  [N]    Root moment per root chord                       
         W_{cap} : 4942       [N]    Weight of spar caps                              
         W_{web} : 913.9      [N]    Weight of shear web                              
             \nu : 0.8612            Dummy variable = $(t^2 + t + 1)/(t+1)$           
         t_{cap} : 0.001127          Non-dim. spar cap thickness                      
         t_{web} : 0.0009258         Non-dim. shear web thickness                     

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
               \mu : 1.4e-05   [N*s/m**2] Dynamic viscosity (35,000ft)                            
              \rho : 0.38      [kg/m**3]  Air density (35,000 ft)                                 
   \tan(\Lambda_h) : 0.5774               tangent of horizontal tail sweep                        
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
          N_{lift} : 0.4073   Wing loading multiplier                                 
           L_{max} : 0.4073   Maximum wing load                                       
                 g : 0.4016   Gravitational acceleration                              
        \rho_{cap} : 0.3389   Density of spar cap material                            
         f_{w,add} : 0.1147   Wing added weight fraction                              
               r_h : 0.06267  Fractional wing thickness at spar web                   
        \rho_{web} : 0.06267  Density of shear web material                           
\sigma_{max,shear} : -0.06267 Allowable shear stress                                  
      \sigma_{max} : -0.3446  Allowable tensile stress                                
                                                                                      
        V_{\infty} : 1.191    Freestream velocity                                     
            x_{CG} : 1.019    CG location                                             
               S_w : 0.7851   Wing area                                               
              \rho : 0.5922   Air density (35,000 ft)                                 
        \Delta x_w : 0.4342   Distance from aerodynamic centre to CG                  
    \bar{c}_{wing} : 0.3587   Mean aerodynamic chord (wing)                           
           C_{L_w} : 0.349    Lift coefficient (wing)                                 
      C_{m_{fuse}} : 0.1745   Moment coefficient (fuselage)                           
        C_{m_{ac}} : 0.1745   Moment coefficient about aerodynamic centre (wing)      
               w_f : 0.1112   Fuselage width                                          
        C_{L_{aw}} : 0.08712  Lift curve slope (wing)                                 
   \tan(\Lambda_h) : -0.08545 tangent of horizontal tail sweep                        
              \eta : -0.09112 Lift efficiency (diff between sectional and actual lift)
          l_{fuse} : -0.9713  Fuselage length                                         

