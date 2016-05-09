Beginning signomial solve.
Solving took 10 GP solves and 1.12 seconds.

Cost
----
 2.597e+04 [N] 

Free Variables
--------------
                 AR : 7.772 Wing aspect ratio 
            C_{D_w} : 0.01361 Drag coefficient  
        C_{D_{p_w}} : 0.005989 Wing parasitic drag coefficient
            C_{L_w} : 0.4259 Lift coefficient (wing)
         C_{L_{aw}} : 4.259 Lift curve slope (wing)
           D_{wing} : 2.597e+04  [N] Wing drag         
                L_w : 8.127e+05  [N] Wing lift         
        L_{max_{w}} : 5.536e+06  [N] Maximum load      
               Re_w : 3.508e+07 Cruise Reynolds number (wing)
                S_w : 174.4  [m**2] Wing area         
      V_{fuel, max} : 180.8  [m**3] Available fuel volume
                  W : 8.127e+05  [N] Aircraft weight   
           W_{wing} : 2.127e+05  [N] Wing weight       
           \alpha_w : 0.1   Wing angle of attack
\bar{A}_{fuel, max} : 0.069 Non-dim. fuel area
     \bar{c}_{wing} : 5.384  [m] Mean aerodynamic chord (wing)
            \lambda : 0.2   Wing taper ratio  
             \tau_w : 0.15  Wing thickness/chord ratio
                b_w : 36.81  [m] Wing span         
           c_{root} : 7.894  [m] Wing root chord   
            c_{tip} : 1.579  [m] Wing tip chord    
                e_w : 0.975 Oswald efficiency factor
       f(\lambda_w) : 0.0033 Empirical efficiency function of taper
                p_w : 1.4   Substituted variable = 1 + 2*taper
                q_w : 1.2   Substituted variable = 1 + taper
        y_{\bar{c}} : 10.52  [m] Spanwise location of mean aerodynamic chord
                            
            WingBox |                         
            I_{cap} : 2.417e-05 Non-dim spar cap area moment of inertia
                M_r : 2.51e+06  [N] Root moment per root chord
            W_{cap} : 1.386e+05  [N] Weight of spar caps
            W_{web} : 1.336e+04  [N] Weight of shear web
                \nu : 0.8612 Dummy variable = $(t^2 + t + 1)/(t+1)$
            t_{cap} : 0.005517 Non-dim. spar cap thickness
            t_{web} : 0.002364 Non-dim. shear web thickness

Constants
---------
      C_{L_{wmax}} : 2.5            Lift coefficient (wing)
        V_{\infty} : 240       [m/s] Freestream velocity   
            V_{ne} : 144       [m/s] Never exceed velocity 
               W_0 : 5e+05     [N]  Weight excluding wing 
          W_{fuel} : 1e+05     [N]  Fuel weight           
    \alpha_{max,w} : 0.1            Max angle of attack   
     \cos(\Lambda) : 0.866          cosine of sweep angle 
            \eta_w : 0.97           Lift efficiency (diff b/w sectional, actual lift)
               \mu : 1.4e-05   [N*s/m**2] Dynamic viscosity (35,000ft)
              \rho : 0.38      [kg/m**3] Air density (35,000 ft)
            \rho_0 : 1.225     [kg/m**3] Air density (0 ft)    
       \rho_{fuel} : 817       [kg/m**3] Density of fuel       
     \tan(\Lambda) : 0.5774         tangent of wing sweep 
                 g : 9.81      [m/s**2] Gravitational acceleration
                 w : 0.5            Wingbox-width-to-chord ratio
                                    
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
          N_{lift} : 0.486    Wing loading multiplier
                 g : 0.4472   Gravitational acceleration
        \rho_{cap} : 0.4079   Density of spar cap material
         f_{w,add} : 0.1278   Wing added weight fraction
               r_h : 0.03933  Fractional wing thickness at spar web
        \rho_{web} : 0.03933  Density of shear web material
                 w : -0.03882 Wingbox-width-to-chord ratio
\sigma_{max,shear} : -0.03933 Allowable shear stress
      \sigma_{max} : -0.4467  Allowable tensile stress
                              
               W_0 : 1.051    Weight excluding wing 
            V_{ne} : 0.9721   Never exceed velocity 
            \rho_0 : 0.486    Air density (0 ft)    
      C_{L_{wmax}} : 0.486    Lift coefficient (wing)
                 g : 0.4472   Gravitational acceleration
          W_{fuel} : 0.2102   Fuel weight           
     \tan(\Lambda) : 0.1171   tangent of wing sweep 
                 w : -0.03882 Wingbox-width-to-chord ratio
               \mu : -0.07545 Dynamic viscosity (35,000ft)
            \eta_w : -0.4685  Lift efficiency (diff b/w sectional, actual lift)
    \alpha_{max,w} : -0.594   Max angle of attack   
              \rho : -0.6331  Air density (35,000 ft)
        V_{\infty} : -1.342   Freestream velocity   

