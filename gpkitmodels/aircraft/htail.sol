Beginning signomial solve.
Solving took 5 GP solves and 0.659 seconds.

Cost
----
 1.84e+04 [N] 

Free Variables
--------------
                AR_h : 4.417      Horizontal tail aspect ratio                     
             C_{D_h} : 0.005527      Horizontal tail drag coefficient                 
         C_{D_{0_h}} : 0.005162      Horizontal tail parasitic drag coefficient       
             C_{L_h} : 0.07063      Lift coefficient (htail)                         
        C_{L_{ah_0}} : 4.362      Isolated lift curve slope (htail)                
          C_{L_{ah}} : 2.181      Lift curve slope (htail)                         
              D_{ht} : 4327  [N] Horizontal tail drag                             
                 K_f : 0.7831      Empirical factor for fuselage-wing interference  
                 L_h : 5.53e+04  [N] Horizontal tail downforce                        
         L_{{max}_h} : 2.536e+06  [N] Maximum load                                     
            Re_{c_h} : 3.01e+07      Cruise Reynolds number (Horizontal tail)         
                S.M. : 0.05      Stability margin                                 
                 S_h : 76.79  [m**2] Horizontal tail area                             
          V_{\infty} : 231.7  [m/s] Freestream velocity                              
              W_{ht} : 2.814e+04  [N] Horizontal tail weight                           
 \Delta x_{{lead}_h} : 13.05  [m] Distance from CG to horizontal tail leading edge 
\Delta x_{{trail}_h} : 20  [m] Distance from CG to horizontal tail trailing edge
              \alpha : 0.03239      Horizontal tail angle of attack                  
        \bar{c}_{ht} : 4.787  [m] Mean aerodynamic chord (ht)                      
           \lambda_h : 0.2      Horizontal tail taper ratio                      
              \tau_h : 0.15      Horizontal tail thickness/chord ratio            
              b_{ht} : 18.42  [m] Horizontal tail span                             
          c_{root_h} : 6.949  [m] Horizontal tail root chord                       
           c_{tip_h} : 1.39  [m] Horizontal tail tip chord                        
                 e_h : 0.9856      Oswald efficiency factor                         
        f(\lambda_h) : 0.0033      Empirical efficiency function of taper           
                 l_h : 17.29  [m] Horizontal tail moment arm                       
              p_{ht} : 1.4      Substituted variable = 1 + 2*taper               
              q_{ht} : 1.2      Substituted variable = 1 + taper                 
                 x_w : 22  [m] Position of wing aerodynamic center              
    y_{\bar{c}_{ht}} : 5.262  [m] Vertical location of mean aerodynamic chord      
                               
             WingBox |                                                          
             I_{cap} : 8.117e-06      Non-dim spar cap area moment of inertia          
                 M_r : 6.533e+05  [N] Root moment per root chord                       
             W_{cap} : 1.704e+04  [N] Weight of spar caps                              
             W_{web} : 3062  [N] Weight of shear web                              
                 \nu : 0.8612      Dummy variable = $(t^2 + t + 1)/(t+1)$           
             t_{cap} : 0.001749      Non-dim. spar cap thickness                      
             t_{web} : 0.001397      Non-dim. shear web thickness                     

Constants
---------
                AR : 9              Wing aspect ratio     
           C_{L_w} : 0.5            Lift coefficient (wing)
        C_{L_{aw}} : 6.283          Lift curve slope (wing)
      C_{L_{hmax}} : 2.6            Max lift coefficient  
      C_{m_{fuse}} : 0.05           Moment coefficient (fuselage)
                 M : 0.78           Cruise Mach number    
        S.M._{min} : 0.05           Minimum stability margin
               S_w : 125       [m**2] Wing area             
            V_{ne} : 144       [m/s] Never exceed velocity 
        \Delta x_w : 2         [m]  Distance from aerodynamic centre to CG
    \alpha_{max,h} : 0.1            Max angle of attack, htail
    \bar{c}_{wing} : 5         [m]  Mean aerodynamic chord (wing)
            \eta_h : 0.97           Lift efficiency (diff between sectional and actual lift)
         \eta_{ht} : 0.9            Tail efficiency       
               \mu : 1.4e-05   [N*s/m**2] Dynamic viscosity (35,000ft)
              \rho : 0.38      [kg/m**3] Air density (35,000 ft)
            \rho_0 : 1.225     [kg/m**3] Air density (0 ft)    
   \tan(\Lambda_h) : 0.5774         tangent of horizontal tail sweep
                 a : 297       [m/s] Speed of sound (35,000 ft)
          l_{fuse} : 40        [m]  Fuselage length       
          w_{fuse} : 6         [m]  Fuselage width        
            x_{CG} : 20        [m]  CG location           
      |C_{m_{ac}}| : 0.1            Moment coefficient about aerodynamic centre (wing)
                                    
           WingBox |                                      
          N_{lift} : 2              Wing loading multiplier
        \rho_{cap} : 2700      [kg/m**3] Density of spar cap material
        \rho_{web} : 2700      [kg/m**3] Density of shear web material
\sigma_{max,shear} : 1.67e+08  [Pa] Allowable shear stress
      \sigma_{max} : 2.5e+08   [Pa] Allowable tensile stress
         f_{w,add} : 0.4            Wing added weight fraction
                 g : 9.81      [m/s**2] Gravitational acceleration
               r_h : 0.75           Fractional wing thickness at spar web
                 w : 0.5            Wingbox-width-to-chord ratio

Sensitivities
-------------
           WingBox |                               
          N_{lift} : 0.7821  Wing loading multiplier
                 g : 0.7648  Gravitational acceleration
        \rho_{cap} : 0.6482  Density of spar cap material
         f_{w,add} : 0.2185  Wing added weight fraction
               r_h : 0.1165  Fractional wing thickness at spar web
        \rho_{web} : 0.1165  Density of shear web material
                 w : -0.01731 Wingbox-width-to-chord ratio
\sigma_{max,shear} : -0.1165 Allowable shear stress
      \sigma_{max} : -0.6656 Allowable tensile stress
                             
            V_{ne} : 1.564   Never exceed velocity 
          w_{fuse} : 1.142   Fuselage width        
               S_w : 0.8951  Wing area             
        C_{L_{aw}} : 0.8951  Lift curve slope (wing)
        \Delta x_w : 0.8755  Distance from aerodynamic centre to CG
            x_{CG} : 0.7984  CG location           
            \rho_0 : 0.7821  Air density (0 ft)    
      C_{L_{hmax}} : 0.7821  Max lift coefficient  
                 a : 0.4659  Speed of sound (35,000 ft)
              \rho : 0.2307  Air density (35,000 ft)
   \tan(\Lambda_h) : 0.1015  tangent of horizontal tail sweep
    \bar{c}_{wing} : 0.09945 Mean aerodynamic chord (wing)
        S.M._{min} : 0.09945 Minimum stability margin
                 M : -0.1989 Cruise Mach number    
            \eta_h : -0.7946 Lift efficiency (diff between sectional and actual lift)
         \eta_{ht} : -1.466  Tail efficiency       
          l_{fuse} : -2.039  Fuselage length       

