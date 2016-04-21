Beginning signomial solve.
Solving took 4 GP solves and 0.266 seconds.

Cost
----
 1.954e+04 [N] 

Free Variables
--------------
          AR : 5.776             Wing aspect ratio                      
     C_{D_w} : 0.05881           Drag coefficient                       
     C_{L_w} : 0.3962            Lift coefficient (wing)                
  C_{L_{aw}} : 3.962             Lift curve slope (wing)                
    D_{wing} : 1.592e+04  [N]    Wing drag                              
         L_w : 1.072e+05  [N]    Wing lift                              
 L_{max_{w}} : 7.854e+05  [N]    Maximum load                           
         S_w : 24.74      [m**2] Wing area                              
           W : 1.072e+05  [N]    Aircraft weight                        
    W_{wing} : 7243       [N]    Wing weight                            
    \alpha_w : 0.1               Wing angle of attack                   
     \lambda : 0.2               Wing taper ratio                       
      \tau_w : 0.15              Wing thickness/chord ratio             
         b_w : 11.95      [m]    Wing span                              
    c_{root} : 3.449      [m]    Wing root chord                        
     c_{tip} : 0.6898     [m]    Wing tip chord                         
         e_w : 0.9812            Oswald efficiency factor               
f(\lambda_w) : 0.003309          Empirical efficiency function of taper 
         p_w : 1.4               Substituted variable = 1 + 2*taper     
         q_w : 1.2               Substituted variable = 1 + taper       
                                                                        
     WingBox |                                                          
     I_{cap} : 1.335e-05         Non-dim spar cap area moment of inertia
         M_r : 2.646e+05  [N]    Root moment per root chord             
     W_{cap} : 4558       [N]    Weight of spar caps                    
     W_{web} : 615.5      [N]    Weight of shear web                    
         \nu : 0.8612            Dummy variable = $(t^2 + t + 1)/(t+1)$ 
     t_{cap} : 0.002928          Non-dim. spar cap thickness            
     t_{web} : 0.001757          Non-dim. shear web thickness           

Constants
---------
       C_{D_{0_w}} : 0.05                Wing parasitic drag coefficient                         
      C_{L_{wmax}} : 2.5                 Lift coefficient (wing)                                 
        V_{\infty} : 240       [m/s]     Freestream velocity                                     
            V_{ne} : 144       [m/s]     Never exceed velocity                                   
               W_0 : 1e+05     [N]       Weight excluding wing                                   
    \alpha_{max,w} : 0.1                 Max angle of attack                                     
            \eta_w : 0.97                Lift efficiency (diff between sectional and actual lift)
              \rho : 0.38      [kg/m**3] Air density (35,000 ft)                                 
            \rho_0 : 1.225     [kg/m**3] Air density (0 ft)                                      
     \tan(\Lambda) : 0.5774              tangent of wing sweep                                   
                                                                                                 
           WingBox |                                                                             
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
          N_{lift} : 0.2784   Wing loading multiplier                                 
                 g : 0.2674   Gravitational acceleration                              
        \rho_{cap} : 0.2356   Density of spar cap material                            
         f_{w,add} : 0.07641  Wing added weight fraction                              
               r_h : 0.03182  Fractional wing thickness at spar web                   
        \rho_{web} : 0.03182  Density of shear web material                           
                 w : -0.01092 Wingbox-width-to-chord ratio                            
\sigma_{max,shear} : -0.03182 Allowable shear stress                                  
      \sigma_{max} : -0.2465  Allowable tensile stress                                
                                                                                      
               W_0 : 1.134    Weight excluding wing                                   
       C_{D_{0_w}} : 0.6926   Wing parasitic drag coefficient                         
            V_{ne} : 0.5567   Never exceed velocity                                   
            \rho_0 : 0.2784   Air density (0 ft)                                      
      C_{L_{wmax}} : 0.2784   Lift coefficient (wing)                                 
     \tan(\Lambda) : 0.1751   tangent of wing sweep                                   
              \rho : -0.4012  Air density (35,000 ft)                                 
            \eta_w : -0.7003  Lift efficiency (diff between sectional and actual lift)
        V_{\infty} : -0.8023  Freestream velocity                                     
    \alpha_{max,w} : -0.9716  Max angle of attack                                     

