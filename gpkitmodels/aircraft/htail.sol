Beginning signomial solve.
Solving took 4 GP solves and 4.95 seconds.

Cost
----
 9086 [N] 

Free Variables
--------------
                AR_h : 3.853             Horizontal tail aspect ratio                     
             C_{D_h} : 0.006225          Horizontal tail drag coefficient                 
         C_{D_{0_h}} : 0.005177          Horizontal tail parasitic drag coefficient       
             C_{L_h} : 0.1119            Lift coefficient (htail)                         
          C_{L_{ah}} : 3.457             Lift curve slope (htail)                         
              D_{ht} : 3285       [N]    Horizontal tail drag                             
                 K_f : 0.7831            Empirical factor for fuselage-wing interference  
         L_{{max}_h} : 1.592e+06  [N]    Maximum load                                     
            Re_{c_h} : 2.646e+07         Cruise Reynolds number (Horizontal tail)         
                S.M. : 0.05              Stability margin                                 
                 S_h : 48.22      [m**2] Horizontal tail area                             
              W_{ht} : 1.16e+04   [N]    Horizontal tail weight                           
 \Delta x_{{lead}_h} : 14.1       [m]    Distance from CG to horizontal tail leading edge 
\Delta x_{{trail}_h} : 20         [m]    Distance from CG to horizontal tail trailing edge
              \alpha : 0.03239           Horizontal tail angle of attack                  
        \bar{c}_{ht} : 4.062      [m]    Mean aerodynamic chord (ht)                      
           \lambda_h : 0.2               Horizontal tail taper ratio                      
              \tau_h : 0.15              Horizontal tail thickness/chord ratio            
              b_{ht} : 13.63      [m]    Horizontal tail span                             
        c_{{root}_h} : 5.896      [m]    Horizontal tail root chord                       
         c_{{tip}_h} : 1.179      [m]    Horizontal tail tip chord                        
                 e_h : 0.9874            Oswald efficiency factor                         
        f(\lambda_h) : 0.003309          Empirical efficiency function of taper           
                 l_h : 17.37      [m]    Horizontal tail moment arm                       
              p_{ht} : 1.4               Substituted variable = 1 + 2*taper               
              q_{ht} : 1.2               Substituted variable = 1 + taper                 
                 x_w : 22         [m]    Position of wing aerodynamic center              
         y_{\bar{c}} : 3.894      [m]    Vertical location of mean aerodynamic chord      
                                                                                          
             WingBox |                                                                    
             I_{cap} : 6.177e-06         Non-dim spar cap area moment of inertia          
                 M_r : 3.579e+05  [N]    Root moment per root chord                       
             W_{cap} : 6863       [N]    Weight of spar caps                              
             W_{web} : 1423       [N]    Weight of shear web                              
                 \nu : 0.8612            Dummy variable = $(t^2 + t + 1)/(t+1)$           
             t_{cap} : 0.001323          Non-dim. spar cap thickness                      
             t_{web} : 0.001219          Non-dim. shear web thickness                     

Constants
---------
        C_{L_{aw}} : 6.283                Lift curve slope (wing)                                 
      C_{L_{hmax}} : 2.6                  Max lift coefficient                                    
        S.M._{min} : 0.05                 Minimum stability margin                                
               S_w : 125       [m**2]     Wing area                                               
        V_{\infty} : 240       [m/s]      Freestream velocity                                     
            V_{ne} : 144       [m/s]      Never exceed velocity                                   
        \Delta x_w : 2         [m]        Distance from aerodynamic centre to CG                  
      \alpha_{max} : 0.1                  Max angle of attack (htail)                             
    \bar{c}_{wing} : 5         [m]        Mean aerodynamic chord (wing)                           
            \eta_h : 0.97                 Lift efficiency (diff between sectional and actual lift)
               \mu : 1.4e-05   [N*s/m**2] Dynamic viscosity (35,000ft)                            
              \rho : 0.38      [kg/m**3]  Air density (35,000 ft)                                 
            \rho_0 : 1.225     [kg/m**3]  Air density (0 ft)                                      
   \tan(\Lambda_h) : 0.5774               tangent of horizontal tail sweep                        
          l_{fuse} : 40        [m]        Fuselage length                                         
          w_{fuse} : 6         [m]        Fuselage width                                          
            x_{CG} : 20        [m]        CG location                                             
                                                                                                  
           WingBox |                                                                              
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
          N_{lift} : 0.649    Wing loading multiplier                                 
                 g : 0.6384   Gravitational acceleration                              
        \rho_{cap} : 0.5288   Density of spar cap material                            
         f_{w,add} : 0.1824   Wing added weight fraction                              
               r_h : 0.1096   Fractional wing thickness at spar web                   
        \rho_{web} : 0.1096   Density of shear web material                           
                 w : -0.01054 Wingbox-width-to-chord ratio                            
\sigma_{max,shear} : -0.1096  Allowable shear stress                                  
      \sigma_{max} : -0.5393  Allowable tensile stress                                
                                                                                      
            V_{ne} : 1.298    Never exceed velocity                                   
          w_{fuse} : 1.014    Fuselage width                                          
               S_w : 0.7947   Wing area                                               
        C_{L_{aw}} : 0.7947   Lift curve slope (wing)                                 
        \Delta x_w : 0.7772   Distance from aerodynamic centre to CG                  
        V_{\infty} : 0.7162   Freestream velocity                                     
            x_{CG} : 0.7088   CG location                                             
            \rho_0 : 0.649    Air density (0 ft)                                      
      C_{L_{hmax}} : 0.649    Max lift coefficient                                    
              \rho : 0.3546   Air density (35,000 ft)                                 
    \bar{c}_{wing} : 0.0883   Mean aerodynamic chord (wing)                           
        S.M._{min} : 0.0883   Minimum stability margin                                
   \tan(\Lambda_h) : 0.01106  tangent of horizontal tail sweep                        
            \eta_h : -0.7815  Lift efficiency (diff between sectional and actual lift)
          l_{fuse} : -1.912   Fuselage length                                         

