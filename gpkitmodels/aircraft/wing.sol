Beginning signomial solve.
Solving took 10 GP solves and 1.1 seconds.
Warning: Constraint [1.01*\bar{c}_{wing} >= 0.667*\...] is not tight because the left hand side evaluated to 7.91484339324 meter but the right hand side evaluated to 7.83647860716 meter (Allowable error: 1.0%, Actual error: 1%)


Cost
----
 1.582e+05 [N] 

Free Variables
--------------
            AR : 3.879             Wing aspect ratio                          
       C_{D_w} : 0.05998           Drag coefficient                           
       C_{L_w} : 0.3466            Lift coefficient (wing)                    
    C_{L_{aw}} : 3.466             Lift curve slope (wing)                    
      D_{wing} : 1.178e+05  [N]    Wing drag                                  
           L_w : 6.808e+05  [N]    Wing lift                                  
   L_{max_{w}} : 5.699e+06  [N]    Maximum load                               
           S_w : 179.5      [m**2] Wing area                                  
      V_{fuel} : 180.8      [m**3] Available fuel volume                      
             W : 6.808e+05  [N]    Aircraft weight                            
      W_{wing} : 8.078e+04  [N]    Wing weight                                
      \alpha_w : 0.1               Wing angle of attack                       
\bar{A}_{fuel} : 0.069             Non-dim. fuel area                         
\bar{c}_{wing} : 7.836      [m]    Mean aerodynamic chord (wing)              
       \lambda : 0.2               Wing taper ratio                           
        \tau_w : 0.15              Wing thickness/chord ratio                 
           b_w : 26.39      [m]    Wing span                                  
      c_{root} : 11.38      [m]    Wing root chord                            
       c_{tip} : 2.275      [m]    Wing tip chord                             
           e_w : 0.9874            Oswald efficiency factor                   
  f(\lambda_w) : 0.0033            Empirical efficiency function of taper     
           p_w : 1.4               Substituted variable = 1 + 2*taper         
           q_w : 1.2               Substituted variable = 1 + taper           
   y_{\bar{c}} : 7.539      [m]    Spanwise location of mean aerodynamic chord
                                                                              
       WingBox |                                                              
       I_{cap} : 6.019e-06         Non-dim spar cap area moment of inertia    
           M_r : 1.289e+06  [N]    Root moment per root chord                 
       W_{cap} : 4.784e+04  [N]    Weight of spar caps                        
       W_{web} : 9859       [N]    Weight of shear web                        
           \nu : 0.8612            Dummy variable = $(t^2 + t + 1)/(t+1)$     
       t_{cap} : 0.001288          Non-dim. spar cap thickness                
       t_{web} : 0.00118           Non-dim. shear web thickness               

Constants
---------
       C_{D_{0_w}} : 0.05                Wing parasitic drag coefficient                  
      C_{L_{wmax}} : 2.5                 Lift coefficient (wing)                          
        V_{\infty} : 240       [m/s]     Freestream velocity                              
            V_{ne} : 144       [m/s]     Never exceed velocity                            
               W_0 : 5e+05     [N]       Weight excluding wing                            
          W_{fuel} : 1e+05     [N]       Fuel weight                                      
    \alpha_{max,w} : 0.1                 Max angle of attack                              
     \cos(\Lambda) : 0.866               cosine of sweep angle                            
            \eta_w : 0.97                Lift efficiency (diff b/w sectional, actual lift)
              \rho : 0.38      [kg/m**3] Air density (35,000 ft)                          
            \rho_0 : 1.225     [kg/m**3] Air density (0 ft)                               
       \rho_{fuel} : 817       [kg/m**3] Density of fuel                                  
     \tan(\Lambda) : 0.5774              tangent of wing sweep                            
                 g : 9.81      [m/s**2]  Gravitational acceleration                       
                 w : 0.5                 Wingbox-width-to-chord ratio                     
                                                                                          
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
          N_{lift} : 0.4248   Wing loading multiplier                          
                 g : 0.4181   Gravitational acceleration                       
        \rho_{cap} : 0.3466   Density of spar cap material                     
         f_{w,add} : 0.1194   Wing added weight fraction                       
               r_h : 0.07143  Fractional wing thickness at spar web            
        \rho_{web} : 0.07143  Density of shear web material                    
\sigma_{max,shear} : -0.07143 Allowable shear stress                           
      \sigma_{max} : -0.3534  Allowable tensile stress                         
                                                                               
               W_0 : 1.008    Weight excluding wing                            
            V_{ne} : 0.8496   Never exceed velocity                            
       C_{D_{0_w}} : 0.6208   Wing parasitic drag coefficient                  
            \rho_0 : 0.4248   Air density (0 ft)                               
      C_{L_{wmax}} : 0.4248   Lift coefficient (wing)                          
                 g : 0.4181   Gravitational acceleration                       
          W_{fuel} : 0.2015   Fuel weight                                      
     \tan(\Lambda) : 0.1693   tangent of wing sweep                            
              \rho : -0.6271  Air density (35,000 ft)                          
            \eta_w : -0.6772  Lift efficiency (diff b/w sectional, actual lift)
    \alpha_{max,w} : -1.124   Max angle of attack                              
        V_{\infty} : -1.254   Freestream velocity                              

