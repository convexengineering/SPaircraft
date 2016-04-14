Beginning signomial solve.
Solving took 5 GP solves and 2.97 seconds.

Cost
----
 1.389e+04 [N] 

Free Variables
--------------
          AR : 12.04             Wing aspect ratio                      
     C_{D_w} : 0.0558            Drag coefficient                       
     C_{L_w} : 0.4593            Lift coefficient (wing)                
  C_{L_{aw}} : 4.593             Lift curve slope (wing)                
         D_h : 1.249e+04  [N]    Wing drag                              
           L : 1.028e+05  [N]    Lift                                   
     L_{max} : 1.193e+05  [N]    Maximum load                           
         S_w : 20.45      [m**2] Wing area                              
           W : 1.028e+05  [N]    Aircraft weight                        
         W_w : 2800       [N]    Wing weight                            
      \alpha : 0.1               Wing angle of attack                   
     \lambda : 0.2               Wing taper ratio                       
        \tau : 0.15              Wing thickness/chord ratio             
         b_w : 15.69      [m]    Wing span                              
    c_{root} : 2.172      [m]    Wing root chord                        
     c_{tip} : 0.4345     [m]    Wing tip chord                         
           e : 0.9618            Oswald efficiency factor               
f\(\lambda\) : 0.0033            Empirical efficiency function of taper 
           p : 1.4               Substituted variable = 1 + 2*taper     
           q : 1.2               Substituted variable = 1 + taper       
                                                                        
     WingBox |                                                          
     I_{cap} : 1.065e-05         Non-dim spar cap area moment of inertia
         M_r : 8.379e+04  [N]    Root moment per root chord             
     W_{cap} : 1877       [N]    Weight of spar caps                    
     W_{web} : 122.7      [N]    Weight of shear web                    
         \nu : 0.8612            Dummy variable = $(t^2 + t + 1)/(t+1)$ 
     t_{cap} : 0.002315          Non-dim. spar cap thickness            
     t_{web} : 0.0006728         Non-dim. shear web thickness           

Constants
---------
       C_{D_{0_w}} : 0.05                Wing parasitic drag coefficient                         
        V_{\infty} : 240       [m/s]     Freestream velocity                                     
            V_{ne} : 144       [m/s]     Never exceed velocity                                   
               W_0 : 1e+05     [N]       Weight excluding wing                                   
      \alpha_{max} : 0.1                 Max angle of attack                                     
              \eta : 0.97                Lift efficiency (diff between sectional and actual lift)
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
     N_{lift} : 0.135   Wing loading multiplier                                 
            g : 0.1306  Gravitational acceleration                              
   \rho_{cap} : 0.1226  Density of spar cap material                            
    f_{w,add} : 0.03732 Wing added weight fraction                              
 \sigma_{max} : -0.127  Allowable tensile stress                                
                                                                                
          W_0 : 1.065   Weight excluding wing                                   
  C_{D_{0_w}} : 0.8058  Wing parasitic drag coefficient                         
       V_{ne} : 0.2701  Never exceed velocity                                   
\tan(\Lambda) : 0.1666  tangent of wing sweep                                   
       \rho_0 : 0.135   Air density (0 ft)                                      
         \rho : -0.1959 Air density (35,000 ft)                                 
   V_{\infty} : -0.3919 Freestream velocity                                     
         \eta : -0.6664 Lift efficiency (diff between sectional and actual lift)
 \alpha_{max} : -0.7732 Max angle of attack                                     

