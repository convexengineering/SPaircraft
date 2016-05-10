Beginning signomial solve.
Solving took 5 GP solves and 0.588 seconds.

Cost
----
 2.433e+04 [N] 

Free Variables
--------------
                 AR : 9.322 Wing aspect ratio 
            C_{D_w} : 0.0169 Drag coefficient  
        C_{D_{p_w}} : 0.005783 Wing parasitic drag coefficient
            C_{L_w} : 0.5619 Lift coefficient (wing)
         C_{L_{aw}} : 5.619 Lift curve slope (wing)
           D_{wing} : 2.433e+04  [N] Wing drag         
                L_w : 8.093e+05  [N] Wing lift         
        L_{max_{w}} : 4.485e+06  [N] Maximum load      
               Re_w : 2.81e+07 Cruise Reynolds number (wing)
                S_w : 141.2  [m**2] Wing area         
         V_{\infty} : 231.7  [m/s] Freestream velocity
      V_{fuel, max} : 180.8  [m**3] Available fuel volume
                  W : 8.093e+05  [N] Aircraft weight   
           W_{wing} : 2.093e+05  [N] Wing weight       
           \alpha_w : 0.1   Wing angle of attack
\bar{A}_{fuel, max} : 0.069 Non-dim. fuel area
     \bar{c}_{wing} : 4.469  [m] Mean aerodynamic chord (wing)
            \lambda : 0.2   Wing taper ratio  
             \tau_w : 0.15  Wing thickness/chord ratio
                b_w : 36.28  [m] Wing span         
           c_{root} : 6.487  [m] Wing root chord   
            c_{tip} : 1.297  [m] Wing tip chord    
                e_w : 0.9702 Oswald efficiency factor
       f(\lambda_w) : 0.0033 Empirical efficiency function of taper
                p_w : 1.4   Substituted variable = 1 + 2*taper
                q_w : 1.2   Substituted variable = 1 + taper
        y_{\bar{c}} : 10.37  [m] Spanwise location of mean aerodynamic chord
                            
            WingBox |                         
            I_{cap} : 3.477e-05 Non-dim spar cap area moment of inertia
                M_r : 2.439e+06  [N] Root moment per root chord
            W_{cap} : 1.388e+05  [N] Weight of spar caps
            W_{web} : 1.067e+04  [N] Weight of shear web
                \nu : 0.8612 Dummy variable = $(t^2 + t + 1)/(t+1)$
            t_{cap} : 0.008301 Non-dim. spar cap thickness
            t_{web} : 0.002836 Non-dim. shear web thickness

Constants
---------
      C_{L_{wmax}} : 2.5            Lift coefficient (wing)
                 M : 0.78           Cruise Mach number    
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
                 a : 297       [m/s] Speed of sound (35,000 ft)
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
          N_{lift} : 0.497    Wing loading multiplier
                 g : 0.4332   Gravitational acceleration
        \rho_{cap} : 0.4023   Density of spar cap material
         f_{w,add} : 0.1238   Wing added weight fraction
               r_h : 0.03093  Fractional wing thickness at spar web
        \rho_{web} : 0.03093  Density of shear web material
\sigma_{max,shear} : -0.03093 Allowable shear stress
                 w : -0.06373 Wingbox-width-to-chord ratio
      \sigma_{max} : -0.4661  Allowable tensile stress
                              
               W_0 : 1.035    Weight excluding wing 
            V_{ne} : 0.994    Never exceed velocity 
            \rho_0 : 0.497    Air density (0 ft)    
      C_{L_{wmax}} : 0.497    Lift coefficient (wing)
                 g : 0.4332   Gravitational acceleration
          W_{fuel} : 0.207    Fuel weight           
     \tan(\Lambda) : 0.1221   tangent of wing sweep 
               \mu : -0.05124 Dynamic viscosity (35,000ft)
                 w : -0.06373 Wingbox-width-to-chord ratio
            \eta_w : -0.2655  Lift efficiency (diff b/w sectional, actual lift)
    \alpha_{max,w} : -0.3482  Max angle of attack   
              \rho : -0.6242  Air density (35,000 ft)
                 a : -1.3     Speed of sound (35,000 ft)
                 M : -1.523   Cruise Mach number    

