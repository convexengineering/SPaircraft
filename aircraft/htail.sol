Using solver 'mosek'
Solving for 21 variables.
Solving took 0.0158 seconds.

Cost
----
 2.335e+04 [N] 

Free Variables
--------------
      AR_h : 4.018             Horizontal tail aspect ratio           
   C_{D_h} : 0.03992           Horizontal tail drag coefficient       
   C_{L_h} : 0.3165            Lift coefficient (htail)               
C_{L_{ah}} : 3.514             Lift curve slope (htail)               
       D_h : 1.725e+04  [N]    Horizontal tail drag                   
      S.M. : 0.05              Stability margin                       
       S_h : 39.49      [m**2] Horizontal tail area                   
         W : 1.22e+04   [N]    Horizontal tail weight                 
    \alpha : 0.09009           Horizontal tail angle of attack        
       l_h : 20         [m]    Horizontal tail moment arm             
                                                                      
   WingBox |                                                          
   I_{cap} : 1.021e-05         Non-dim spar cap area moment of inertia
       M_r : 3.181e+05  [N]    Root moment per root chord             
   W_{cap} : 7611       [N]    Weight of spar caps                    
   W_{web} : 1100       [N]    Weight of shear web                    
       \nu : 0.786             Dummy variable = $(t^2 + t + 1)/(t+1)$ 
      \tau : 0.15              Thickness to chord ratio               
         b : 12.6       [m]    Span                                   
         p : 1.9               Substituted variable = 1 + 2*taper     
         q : 1.45              Substituted variable = 1 + taper       
   t_{cap} : 0.002215          Non-dim. spar cap thickness            
   t_{web} : 0.001423          Non-dim. shear web thickness           

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
           \bar{c} : 5         [m]       Mean aerodynamic chord (wing)                           
              \eta : 0.97                Lift efficiency (diff between sectional and actual lift)
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
             taper : 0.45                Taper ratio                                             
                 w : 0.5                 Wingbox-width-to-chord ratio                            

Sensitivities
-------------
           WingBox |                                                                  
          N_{lift} : 0.269    Wing loading multiplier                                 
           L_{max} : 0.269    Maximum wing load                                       
                 g : 0.2611   Gravitational acceleration                              
             taper : 0.2518   Taper ratio                                             
        \rho_{cap} : 0.2281   Density of spar cap material                            
         f_{w,add} : 0.07461  Wing added weight fraction                              
               r_h : 0.03299  Fractional wing thickness at spar web                   
        \rho_{web} : 0.03299  Density of shear web material                           
\sigma_{max,shear} : -0.03299 Allowable shear stress                                  
      \sigma_{max} : -0.236   Allowable tensile stress                                
                                                                                      
        V_{\infty} : 1.478    Freestream velocity                                     
              \rho : 0.7389   Air density (35,000 ft)                                 
               S_w : 0.6821   Wing area                                               
       C_{D_{0_h}} : 0.5552   Horizontal tail parasitic drag coefficient              
               x_w : 0.4634   Distance from aerodynamic centre to CG                  
               w_f : 0.3591   Fuselage width                                          
        C_{L_{aw}} : 0.3148   Lift curve slope (wing)                                 
           \bar{c} : 0.2186   Mean aerodynamic chord (wing)                           
           C_{L_w} : 0.1837   Lift coefficient (wing)                                 
               K_f : 0.1795   Empirical factor for fuselage-wing interference         
      C_{m_{fuse}} : 0.09183  Moment coefficient (fuselage)                           
        C_{m_{ac}} : 0.09183  Moment coefficient about aerodynamic centre (wing)      
   \tan(\Lambda_h) : 0.0759   tangent of horizontal tail sweep                        
        S.M._{min} : 0.03497  Minimum stability margin                                
               e_h : -0.1837  Oswald efficiency factor                                
              \eta : -0.3036  Lift efficiency (diff between sectional and actual lift)
          l_{fuse} : -0.6821  Fuselage length                                         

