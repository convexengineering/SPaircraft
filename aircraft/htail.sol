Using solver 'mosek'
Solving for 7 variables.
Solving took 0.0171 seconds.

Cost
----
 3956 [N] 

Free Variables
--------------
   C_{L_h} : 0.3508         Lift coefficient (htail)       
C_{L_{ah}} : 3.508          Lift curve slope (htail)       
      S.M. : 0.05           Stability margin               
       S_h : 39.56   [m**2] Horizontal tail area           
         W : 3956    [N]    Horizontal tail weight         
    \alpha : 0.1            Horizontal tail angle of attack
       l_h : 20      [m]    Horizontal tail moment arm     

Constants
---------
           AR_h : 4              Horizontal tail aspect ratio                            
        C_{L_w} : 0.5            Lift coefficient (wing)                                 
     C_{L_{aw}} : 6.283          Lift curve slope (wing)                                 
     C_{m_{ac}} : 0.1            Moment coefficient about aerodynamic centre (wing)      
   C_{m_{fuse}} : 0.1            Moment coefficient (fuselage)                           
            K_f : 0.7            Empirical factor for fuselage-wing interference         
     S.M._{min} : 0.05           Minimum stability margin                                
            S_w : 125     [m**2] Wing area                                               
   \alpha_{max} : 0.1            Max angle of attack (htail)                             
        \bar{c} : 5       [m]    Mean aerodynamic chord (wing)                           
           \eta : 0.97           Lift efficiency (diff between sectional and actual lift)
\tan(\Lambda_h) : 0.5774         tangent of horizontal tail sweep                        
       l_{fuse} : 40      [m]    Fuselage length                                         
            w_f : 6       [m]    Fuselage width                                          
            x_w : 2       [m]    Distance from aerodynamic centre to CG                  
         x_{CG} : 20      [m]    CG location                                             

Sensitivities
-------------
            w_f : 0.7264  Fuselage width                                          
            S_w : 0.6368  Wing area                                               
     C_{L_{aw}} : 0.6368  Lift curve slope (wing)                                 
            x_w : 0.566   Distance from aerodynamic centre to CG                  
            K_f : 0.3632  Empirical factor for fuselage-wing interference         
\tan(\Lambda_h) : 0.1532  tangent of horizontal tail sweep                        
        \bar{c} : 0.07075 Mean aerodynamic chord (wing)                           
     S.M._{min} : 0.07075 Minimum stability margin                                
           AR_h : -0.3872 Horizontal tail aspect ratio                            
           \eta : -0.6128 Lift efficiency (diff between sectional and actual lift)
       l_{fuse} : -0.6368 Fuselage length                                         

