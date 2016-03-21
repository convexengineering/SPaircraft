Using solver 'mosek'
Solving for 9 variables.
Solving took 0.0163 seconds.

Cost
----
 1.926e+04 [N] 

Free Variables
--------------
   C_{D_h} : 0.03993           Horizontal tail drag coefficient
   C_{L_h} : 0.316             Lift coefficient (htail)        
C_{L_{ah}} : 3.508             Lift curve slope (htail)        
       D_h : 1.729e+04  [N]    Horizontal tail drag            
      S.M. : 0.05              Stability margin                
       S_h : 39.56      [m**2] Horizontal tail area            
         W : 3956       [N]    Horizontal tail weight          
    \alpha : 0.09009           Horizontal tail angle of attack 
       l_h : 20         [m]    Horizontal tail moment arm      

Constants
---------
           AR_h : 4                 Horizontal tail aspect ratio                            
    C_{D_{0_h}} : 0.03              Horizontal tail parasitic drag coefficient              
        C_{L_w} : 0.5               Lift coefficient (wing)                                 
     C_{L_{aw}} : 6.283             Lift curve slope (wing)                                 
     C_{m_{ac}} : 0.1               Moment coefficient about aerodynamic centre (wing)      
   C_{m_{fuse}} : 0.1               Moment coefficient (fuselage)                           
            K_f : 0.7               Empirical factor for fuselage-wing interference         
     S.M._{min} : 0.05              Minimum stability margin                                
            S_w : 125     [m**2]    Wing area                                               
     V_{\infty} : 240     [m/s]     Freestream velocity                                     
   \alpha_{max} : 0.1               Max angle of attack (htail)                             
        \bar{c} : 5       [m]       Mean aerodynamic chord (wing)                           
           \eta : 0.97              Lift efficiency (diff between sectional and actual lift)
           \rho : 0.38    [kg/m**3] Air density (35,000 ft)                                 
\tan(\Lambda_h) : 0.5774            tangent of horizontal tail sweep                        
            e_h : 0.8               Oswald efficiency factor                                
       l_{fuse} : 40      [m]       Fuselage length                                         
            w_f : 6       [m]       Fuselage width                                          
            x_w : 2       [m]       Distance from aerodynamic centre to CG                  
         x_{CG} : 20      [m]       CG location                                             

Sensitivities
-------------
     V_{\infty} : 1.795   Freestream velocity                                     
           \rho : 0.8973  Air density (35,000 ft)                                 
            S_w : 0.7989  Wing area                                               
    C_{D_{0_h}} : 0.6741  Horizontal tail parasitic drag coefficient              
            x_w : 0.5365  Distance from aerodynamic centre to CG                  
            w_f : 0.4021  Fuselage width                                          
     C_{L_{aw}} : 0.3525  Lift curve slope (wing)                                 
        \bar{c} : 0.2624  Mean aerodynamic chord (wing)                           
        C_{L_w} : 0.2232  Lift coefficient (wing)                                 
            K_f : 0.2011  Empirical factor for fuselage-wing interference         
   C_{m_{fuse}} : 0.1116  Moment coefficient (fuselage)                           
     C_{m_{ac}} : 0.1116  Moment coefficient about aerodynamic centre (wing)      
\tan(\Lambda_h) : 0.0848  tangent of horizontal tail sweep                        
     S.M._{min} : 0.03917 Minimum stability margin                                
            e_h : -0.2232 Oswald efficiency factor                                
           \eta : -0.3392 Lift efficiency (diff between sectional and actual lift)
           AR_h : -0.4376 Horizontal tail aspect ratio                            
       l_{fuse} : -0.7989 Fuselage length                                         

