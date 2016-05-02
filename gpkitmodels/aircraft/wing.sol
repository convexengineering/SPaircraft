Beginning signomial solve.
Solving took 10 GP solves and 1.16 seconds.
Warning: Constraint [1.01*\bar{c}_{wing} >= 0.667*\...] is not tight because the left hand side evaluated to 6.98747787557 meter but the right hand side evaluated to 6.9182949263 meter (Allowable error: 1.0%, Actual error: 1%)


Cost
----
 1.27e+05 [N] 

Free Variables
--------------
            AR : 4.046             Wing aspect ratio                          
       C_{D_w} : 0.0599            Drag coefficient                           
       C_{L_w} : 0.3523            Lift coefficient (wing)                    
    C_{L_{aw}} : 3.523             Lift curve slope (wing)                    
      D_{wing} : 9.566e+04  [N]    Wing drag                                  
           L_w : 5.627e+05  [N]    Wing lift                                  
   L_{max_{w}} : 4.634e+06  [N]    Maximum load                               
           S_w : 145.9      [m**2] Wing area                                  
             W : 5.627e+05  [N]    Aircraft weight                            
      W_{wing} : 6.275e+04  [N]    Wing weight                                
      \alpha_w : 0.1               Wing angle of attack                       
\bar{c}_{wing} : 6.918      [m]    Mean aerodynamic chord (wing)              
       \lambda : 0.2               Wing taper ratio                           
        \tau_w : 0.15              Wing thickness/chord ratio                 
           b_w : 24.3       [m]    Wing span                                  
      c_{root} : 10.04      [m]    Wing root chord                            
       c_{tip} : 2.009      [m]    Wing tip chord                             
           e_w : 0.9868            Oswald efficiency factor                   
  f(\lambda_w) : 0.0033            Empirical efficiency function of taper     
           p_w : 1.4               Substituted variable = 1 + 2*taper         
           q_w : 1.2               Substituted variable = 1 + taper           
   y_{\bar{c}} : 6.943      [m]    Spanwise location of mean aerodynamic chord
                                                                              
       WingBox |                                                              
       I_{cap} : 6.55e-06          Non-dim spar cap area moment of inertia    
           M_r : 1.094e+06  [N]    Root moment per root chord                 
       W_{cap} : 3.744e+04  [N]    Weight of spar caps                        
       W_{web} : 7383       [N]    Weight of shear web                        
           \nu : 0.8612            Dummy variable = $(t^2 + t + 1)/(t+1)$     
       t_{cap} : 0.001404          Non-dim. spar cap thickness                
       t_{web} : 0.001231          Non-dim. shear web thickness               

Constants
---------
       C_{D_{0_w}} : 0.05                Wing parasitic drag coefficient                  
      C_{L_{wmax}} : 2.5                 Lift coefficient (wing)                          
        V_{\infty} : 240       [m/s]     Freestream velocity                              
            V_{ne} : 144       [m/s]     Never exceed velocity                            
               W_0 : 5e+05     [N]       Weight excluding wing                            
    \alpha_{max,w} : 0.1                 Max angle of attack                              
            \eta_w : 0.97                Lift efficiency (diff b/w sectional, actual lift)
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
          N_{lift} : 0.4044   Wing loading multiplier                          
                 g : 0.3974   Gravitational acceleration                       
        \rho_{cap} : 0.3319   Density of spar cap material                     
         f_{w,add} : 0.1135   Wing added weight fraction                       
               r_h : 0.06546  Fractional wing thickness at spar web            
        \rho_{web} : 0.06546  Density of shear web material                    
\sigma_{max,shear} : -0.06546 Allowable shear stress                           
      \sigma_{max} : -0.339   Allowable tensile stress                         
                                                                               
               W_0 : 1.199    Weight excluding wing                            
            V_{ne} : 0.8089   Never exceed velocity                            
       C_{D_{0_w}} : 0.6286   Wing parasitic drag coefficient                  
            \rho_0 : 0.4044   Air density (0 ft)                               
      C_{L_{wmax}} : 0.4044   Lift coefficient (wing)                          
     \tan(\Lambda) : 0.1696   tangent of wing sweep                            
              \rho : -0.5961  Air density (35,000 ft)                          
            \eta_w : -0.6783  Lift efficiency (diff b/w sectional, actual lift)
    \alpha_{max,w} : -1.1     Max angle of attack                              
        V_{\infty} : -1.192   Freestream velocity                              

