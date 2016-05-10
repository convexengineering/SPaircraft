Beginning signomial solve.
Solving took 5 GP solves and 0.584 seconds.

Cost
----
 2.315e+04 [N] 

Free Variables
--------------
                 AR : 9.985 Wing aspect ratio 
            C_{D_w} : 0.01685 Drag coefficient  
        C_{D_{p_w}} : 0.005742 Wing parasitic drag coefficient
            C_{L_w} : 0.5808 Lift coefficient (wing)
         C_{L_{aw}} : 5.808 Lift curve slope (wing)
           D_{wing} : 2.315e+04  [N] Wing drag         
                L_w : 7.979e+05  [N] Wing lift         
        L_{max_{w}} : 3.986e+06  [N] Maximum load      
               Re_w : 2.652e+07 Cruise Reynolds number (wing)
                S_w : 125.5  [m**2] Wing area         
      V_{fuel, max} : 180.8  [m**3] Available fuel volume
                  W : 7.979e+05  [N] Aircraft weight   
           W_{wing} : 1.979e+05  [N] Wing weight       
           \alpha_w : 0.1   Wing angle of attack
\bar{A}_{fuel, max} : 0.069 Non-dim. fuel area
     \bar{c}_{wing} : 4.071  [m] Mean aerodynamic chord (wing)
            \lambda : 0.2   Wing taper ratio  
             \tau_w : 0.15  Wing thickness/chord ratio
                b_w : 35.4   [m] Wing span         
           c_{root} : 5.909  [m] Wing root chord   
            c_{tip} : 1.182  [m] Wing tip chord    
                e_w : 0.9681 Oswald efficiency factor
       f(\lambda_w) : 0.0033 Empirical efficiency function of taper
                p_w : 1.4   Substituted variable = 1 + 2*taper
                q_w : 1.2   Substituted variable = 1 + taper
        y_{\bar{c}} : 10.12  [m] Spanwise location of mean aerodynamic chord
                            
            WingBox |                         
            I_{cap} : 3.989e-05 Non-dim spar cap area moment of inertia
                M_r : 2.321e+06  [N] Root moment per root chord
            W_{cap} : 1.321e+05  [N] Weight of spar caps
            W_{web} : 9251   [N] Weight of shear web
                \nu : 0.8612 Dummy variable = $(t^2 + t + 1)/(t+1)$
            t_{cap} : 0.009758 Non-dim. spar cap thickness
            t_{web} : 0.003038 Non-dim. shear web thickness

Constants
---------
      C_{L_{wmax}} : 2.5            Lift coefficient (wing)
                 M : 0.8            Mach number           
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
          N_{lift} : 0.479    Wing loading multiplier
                 g : 0.4045   Gravitational acceleration
        \rho_{cap} : 0.378    Density of spar cap material
         f_{w,add} : 0.1156   Wing added weight fraction
               r_h : 0.02647  Fractional wing thickness at spar web
        \rho_{web} : 0.02647  Density of shear web material
\sigma_{max,shear} : -0.02647 Allowable shear stress
                 w : -0.07454 Wingbox-width-to-chord ratio
      \sigma_{max} : -0.4525  Allowable tensile stress
                              
               W_0 : 1.022    Weight excluding wing 
            V_{ne} : 0.958    Never exceed velocity 
            \rho_0 : 0.479    Air density (0 ft)    
      C_{L_{wmax}} : 0.479    Lift coefficient (wing)
                 g : 0.4045   Gravitational acceleration
          W_{fuel} : 0.2044   Fuel weight           
     \tan(\Lambda) : 0.1094   tangent of wing sweep 
               \mu : -0.04828 Dynamic viscosity (35,000ft)
                 w : -0.07454 Wingbox-width-to-chord ratio
                 M : -0.2101  Mach number           
            \eta_w : -0.2277  Lift efficiency (diff b/w sectional, actual lift)
    \alpha_{max,w} : -0.2946  Max angle of attack   
              \rho : -0.5825  Air density (35,000 ft)
        V_{\infty} : -1.213   Freestream velocity   

