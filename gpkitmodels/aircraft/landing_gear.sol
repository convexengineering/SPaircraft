Beginning signomial solve.
Solving took 7 GP solves and 2.67 seconds.

Cost
----
 1.463e+04 [N] 

Free Variables
--------------
          B : 14.05      [m]    Landing gear base                               
   E_{land} : 3.809e+05  [J]    Max KE to be absorbed in landing                
    F_{w_m} : 7159              Weight factor (main)                            
    F_{w_n} : 643.7             Weight factor (nose)                            
        I_m : 7.934e-06  [m**4] Area moment of inertia (main strut)             
        I_n : 9.614e-07  [m**4] Area moment of inertia (nose strut)             
        L_m : 6.435e+05  [N]    Max static load through main gear               
        L_n : 1.609e+05  [N]    Min static load through nose gear               
L_{n_{dyn}} : 7.77e+04   [N]    Dyn. braking load, nose gear                    
    L_{w_m} : 1.609e+05  [N]    Static load per wheel (main)                    
    L_{w_n} : 8.044e+04  [N]    Static load per wheel (nose)                    
       S_sa : 0.2959     [m]    Stroke of the shock absorber                    
          T : 5.688      [m]    Main landing gear track                         
     W_{lg} : 1.463e+04  [N]    Weight of landing gear                          
     W_{mg} : 1.314e+04  [N]    Weight of main gear                             
     W_{ms} : 626.4      [N]    Weight of main struts                           
     W_{mw} : 2377       [N]    Weight of main wheels (per strut)               
     W_{ng} : 1498       [N]    Weight of nose gear                             
     W_{ns} : 127.2      [N]    Weight of nose strut                            
     W_{nw} : 548.1      [N]    Weight of nose wheels (total)                   
   W_{wa,m} : 267.2      [lbf]  Wheel assembly weight for single main gear wheel
   W_{wa,n} : 61.61      [lbf]  Wheel assembly weight for single nose gear wheel
 \Delta x_m : 2.809      [m]    Distance b/w main gear and CG                   
 \Delta x_n : 11.24      [m]    Distance b/w nose gear and CG                   
 \tan(\phi) : 0.2679            Angle b/w main gear and CG                      
 \tan(\psi) : 1.963             Tip over angles                                 
d_{nacelle} : 2.05       [m]    Nacelle diameter                                
   d_{oleo} : 0.3735     [m]    Diameter of oleo shock absorber                 
    d_{t_m} : 44.5       [in]   Diameter of main gear tires                     
    d_{t_n} : 35.6       [in]   Diameter of nose gear tires                     
        l_m : 2.376      [m]    Length of main gear                             
        l_n : 1.627      [m]    Length of nose gear                             
   l_{oleo} : 0.7399     [m]    Length of oleo shock absorber                   
        r_m : 0.06809    [m]    Radius of main gear struts                      
        r_n : 0.04352    [m]    Radius of nose gear struts                      
        t_m : 0.008002   [m]    Thickness of main gear strut wall               
        t_n : 0.003713   [m]    Thickness of nose gear strut wall               
    w_{t_m} : 0.4088     [m]    Width of main tires                             
    w_{t_n} : 0.327      [m]    Width of nose tires                             
        x_m : 20.83      [m]    x-location of main gear                         
        x_n : 6.788      [m]    x-location of nose gear                         
     x_{CG} : 18.02      [m]    x-location of CG incl. LG                       
        y_m : 2.844      [m]    y-location of main gear (symmetric)             

Constants
---------
                E : 205        [GPa]       Modulus of elasticity, 4340 steel          
                K : 2                      Column effective length factor             
              N_s : 2                      Factor of safety                           
              W_0 : 8.044e+05  [N]         Weight of aircraft excluding landing gear  
           \eta_s : 0.8                    Shock absorber efficiency                  
     \lambda_{LG} : 2.5                    Ratio of max to static load                
        \rho_{st} : 7850       [kg/m**3]   Density of 4340 Steel                      
     \sigma_{y_c} : 4.7e+08    [Pa]        Compressive yield strength 4340 steel      
     \tan(\gamma) : 0.08749                Tangent, dihedral angle                    
 \tan(\phi_{min}) : 0.2679                 Lower bound on phi                         
 \tan(\psi_{max}) : 1.963                  Upper bound on psi                         
\tan(\theta_{TO}) : 0.2679                 Takeoff pitch angle                        
          d_{fan} : 1.75       [m]         Fan diameter                               
        f_{add,m} : 1.5                    Proportional added weight, main            
        f_{add,n} : 1.5                    Proportional added weight, nose            
                g : 9.81       [m/s**2]    Gravitational acceleration                 
         h_{hold} : 1          [m]         Hold height                                
      h_{nacelle} : 0.5        [m]         Min. nacelle clearance                     
           n_{mg} : 2                      Number of main gear struts                 
          n_{wps} : 2                      Number of wheels per strut                 
         p_{oleo} : 1800       [lbf/in**2] Oleo pressure                              
      t_{nacelle} : 0.15       [m]         Nacelle thickness                          
          w_{ult} : 10         [ft/s]      Ultimate velocity of descent               
         x_{CG_0} : 18         [m]         x-location of CG excl. LG                  
           x_{up} : 28         [m]         Fuselage upsweep point                     
          y_{eng} : 4.83       [m]         Spanwise loc. of engines                   
           z_{CG} : 2          [m]         CG height relative to bottom of fuselage   
         z_{wing} : 0.5        [m]         Height of wing relative to base of fuselage

Sensitivities
-------------
             W_0 : 0.8196   Weight of aircraft excluding landing gear
       f_{add,m} : 0.4872   Proportional added weight, main          
         n_{wps} : 0.1804   Number of wheels per strut               
          n_{mg} : 0.1438   Number of main gear struts               
       \rho_{st} : 0.0943   Density of 4340 Steel                    
               g : 0.0943   Gravitational acceleration               
             N_s : 0.0943   Factor of safety                         
    \lambda_{LG} : 0.0856   Ratio of max to static load              
         d_{fan} : 0.07756  Fan diameter                             
       f_{add,n} : 0.05619  Proportional added weight, nose          
     h_{nacelle} : 0.02216  Min. nacelle clearance                   
     t_{nacelle} : 0.0133   Nacelle thickness                        
\tan(\psi_{max}) : -0.0101  Upper bound on psi                       
         y_{eng} : -0.01873 Spanwise loc. of engines                 
    \sigma_{y_c} : -0.0943  Compressive yield strength 4340 steel    

