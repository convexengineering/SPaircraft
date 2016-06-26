Beginning signomial solve.
Solving took 5 GP solves and 0.633 seconds.

Cost
----
 3589 [N] 

Free Variables
--------------
                AR_h : 6.244      Horizontal tail aspect ratio                     
             C_{D_h} : 0.005964      Horizontal tail drag coefficient                 
         C_{D_{0_h}} : 0.005237      Horizontal tail parasitic drag coefficient       
             C_{L_h} : 0.1182      Lift coefficient (htail)                         
        C_{L_{ah_0}} : 5.008      Isolated lift curve slope (htail)                
          C_{L_{ah}} : 2.913      Lift curve slope (htail)                         
              D_{ht} : 2176  [N] Horizontal tail drag                             
                 K_f : 0.6274      Empirical factor for fuselage-wing interference  
                 L_h : 4.313e+04  [N] Horizontal tail downforce                        
         L_{{max}_h} : 1.136e+06  [N] Maximum load                                     
            Re_{c_h} : 1.728e+07      Cruise Reynolds number (Horizontal tail)         
                S.M. : 0.05      Stability margin                                 
                 S_h : 35.78  [m**2] Horizontal tail area                             
          V_{\infty} : 231.7  [m/s] Freestream velocity                              
              W_{ht} : 1.413e+04  [N] Horizontal tail weight                           
 \Delta x_{{lead}_h} : 19.01  [m] Distance from CG to horizontal tail leading edge 
\Delta x_{{trail}_h} : 23  [m] Distance from CG to horizontal tail trailing edge
              \alpha : 0.04059      Horizontal tail angle of attack                  
        \bar{c}_{ht} : 2.748  [m] Mean aerodynamic chord (ht)                      
           \lambda_h : 0.2      Horizontal tail taper ratio                      
              \tau_h : 0.15      Horizontal tail thickness/chord ratio            
              b_{ht} : 14.95  [m] Horizontal tail span                             
          c_{root_h} : 3.99  [m] Horizontal tail root chord                       
           c_{tip_h} : 0.7979  [m] Horizontal tail tip chord                        
                 e_h : 0.9798      Oswald efficiency factor                         
        f(\lambda_h) : 0.0033      Empirical efficiency function of taper           
              l_{ht} : 22.16  [m] Horizontal tail moment arm                       
              p_{ht} : 1.4      Substituted variable = 1 + 2*taper               
              q_{ht} : 1.2      Substituted variable = 1 + taper                 
                 x_w : 19  [m] Position of wing aerodynamic center              
    y_{\bar{c}_{ht}} : 4.27  [m] Spanwise location of mean aerodynamic chord      
                               
             WingBox |                                                          
             I_{cap} : 1.56e-05      Non-dim spar cap area moment of inertia          
                 M_r : 4.137e+05  [N] Root moment per root chord                       
             W_{cap} : 8981  [N] Weight of spar caps                              
             W_{web} : 1113  [N] Weight of shear web                              
                 \nu : 0.8612      Dummy variable = $(t^2 + t + 1)/(t+1)$           
             t_{cap} : 0.003448      Non-dim. spar cap thickness                      
             t_{web} : 0.001899      Non-dim. shear web thickness                     

Constants
---------
              AR_w : 9              Wing aspect ratio     
           C_{L_w} : 0.5            Lift coefficient (wing)
        C_{L_{aw}} : 5              Lift curve slope (wing)
      C_{L_{hmax}} : 2.5            Max lift coefficient  
      C_{m_{fuse}} : 0.05           Moment coefficient (fuselage)
                 M : 0.78           Cruise Mach number    
        S.M._{min} : 0.05           Minimum stability margin
               S_w : 125       [m**2] Wing area             
            V_{ne} : 144       [m/s] Never exceed velocity 
        \Delta x_w : 2         [m]  Distance from aerodynamic centre to CG
    \alpha_{max,h} : 0.1            Max angle of attack, htail
         \bar{c}_w : 5         [m]  Mean aerodynamic chord (wing)
            \eta_h : 0.97           Lift efficiency (diff between sectional and actual lift)
         \eta_{ht} : 0.9            Tail efficiency       
               \mu : 1.4e-05   [N*s/m**2] Dynamic viscosity (35,000 ft)
              \rho : 0.38      [kg/m**3] Air density (35,000 ft)
            \rho_0 : 1.225     [kg/m**3] Air density (0 ft)    
\tan(\Lambda_{ht}) : 0.5774         tangent of horizontal tail sweep
                 a : 297       [m/s] Speed of sound (35,000 ft)
          l_{fuse} : 40        [m]  Fuselage length       
          w_{fuse} : 6         [m]  Fuselage width        
            x_{CG} : 17        [m]  CG location           
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
          N_{lift} : 0.4132   Wing loading multiplier
                 g : 0.3938   Gravitational acceleration
        \rho_{cap} : 0.3503   Density of spar cap material
         f_{w,add} : 0.1125   Wing added weight fraction
               r_h : 0.04342  Fractional wing thickness at spar web
        \rho_{web} : 0.04342  Density of shear web material
                 w : -0.01945 Wingbox-width-to-chord ratio
\sigma_{max,shear} : -0.04342 Allowable shear stress
      \sigma_{max} : -0.3698  Allowable tensile stress
                              
                 a : 1.196    Speed of sound (35,000 ft)
          w_{fuse} : 0.8319   Fuselage width        
            V_{ne} : 0.8264   Never exceed velocity 
               S_w : 0.6475   Wing area             
        C_{L_{aw}} : 0.6475   Lift curve slope (wing)
        \Delta x_w : 0.6407   Distance from aerodynamic centre to CG
                 M : 0.6095   Cruise Mach number    
              \rho : 0.5894   Air density (35,000 ft)
            x_{CG} : 0.5535   CG location           
            \rho_0 : 0.4132   Air density (0 ft)    
      C_{L_{hmax}} : 0.4132   Max lift coefficient  
\tan(\Lambda_{ht}) : 0.1865   tangent of horizontal tail sweep
         \bar{c}_w : 0.07194  Mean aerodynamic chord (wing)
        S.M._{min} : 0.07194  Minimum stability margin
               \mu : 0.01683  Dynamic viscosity (35,000 ft)
            \eta_h : -0.699   Lift efficiency (diff between sectional and actual lift)
         \eta_{ht} : -1.063   Tail efficiency       
          l_{fuse} : -1.46    Fuselage length       

