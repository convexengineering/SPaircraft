Beginning signomial solve.
Solving took 4 GP solves and 0.444 seconds.

Cost
----
 1809 [N] 

Free Variables
--------------
           A_{fan} : 2.405  [m**2] Engine reference area                          
            A_{vt} : 1.948      Vertical tail aspect ratio                     
       C_{D_{vis}} : 0.00508      Viscous drag coefficient                       
        C_{L_{vt}} : 0.4537      Vertical tail lift coefficient                 
            D_{vt} : 1390  [N] Vertical tail viscous drag, cruise             
            D_{wm} : 3609  [N] Engine out windmill drag                       
      L_{max_{vt}} : 1.737e+06  [N] Maximum load for structural sizing             
       L_{v_{max}} : 8.686e+05  [N] Maximum load for structural sizing             
            L_{vt} : 3.581e+04  [N] Vertical tail lift in engine out               
           Re_{vt} : 2.591e+07      Vertical tail reynolds number, cruise          
                 S : 52.61  [m**2] Vertical tail reference area (full)            
            S_{vt} : 26.3  [m**2] Vertical tail reference area (half)            
        W_{struct} : 1.675e+04  [N] Full span weight                               
            W_{vt} : 8374  [N] Vertical tail weight                           
 \Delta x_{lead_v} : 15.21  [m] Distance from CG to vertical tail leading edge 
\Delta x_{trail_v} : 21  [m] Distance from CG to vertical tail trailing edge
      \bar{c}_{vt} : 4.079  [m] Vertical tail mean aero chord                  
      \lambda_{vt} : 0.27      Vertical tail taper ratio                      
         \tau_{vt} : 0.1406      Vertical tail thickness/chord ratio            
                 b : 14.32  [m] Vertical tail full span                        
            b_{vt} : 7.158  [m] Vertical tail half span                        
     c_{root_{vt}} : 5.787  [m] Vertical tail root chord                       
      c_{tip_{vt}} : 1.562  [m] Vertical tail tip chord                        
            l_{vt} : 17.88  [m] Vertical tail moment arm                       
            p_{vt} : 1.54      Substituted variable = 1 + 2*taper             
            q_{vt} : 1.27      Substituted variable = 1 + taper               
  z_{\bar{c}_{vt}} : 1.968  [m] Vertical location of mean aerodynamic chord    
                             
           WingBox |                                                        
                 A : 3.896      Aspect ratio                                   
           I_{cap} : 7.292e-06      Non-dim spar cap area moment of inertia        
               M_r : 4.343e+05  [N] Root moment per root chord                     
           W_{cap} : 1.02e+04  [N] Weight of spar caps                            
        W_{struct} : 1.675e+04  [N] Structural weight                              
           W_{web} : 1766  [N] Weight of shear web                            
               \nu : 0.8327      Dummy variable = $(t^2 + t + 1)/(t+1)$         
           t_{cap} : 0.001794      Non-dim. spar cap thickness                    
           t_{web} : 0.001473      Non-dim. shear web thickness                   

Constants
---------
        C_{D_{wm}} : 0.5            Windmill drag coefficient
      C_{L_{vmax}} : 2.6            Max lift coefficient  
               T_e : 1.29e+05  [N]  Thrust per engine at takeoff
               V_1 : 70        [m/s] Minimum takeoff velocity
        V_{\infty} : 234       [m/s] Cruise velocity       
            V_{ne} : 144       [m/s] Never exceed velocity 
               \mu : 1.4e-05   [N*s/m**2] Dynamic viscosity (35,000 ft)
            \rho_c : 0.38      [kg/m**3] Air density (35,000ft)
         \rho_{TO} : 1.225     [kg/m**3] Air density (SL))     
\tan(\Lambda_{vt}) : 0.8391         Tangent of leading edge sweep (40 deg)
        c_{l_{vt}} : 0.5            Sectional lift force coefficient (engine out)
           d_{fan} : 1.75      [m]  Fan diameter          
               e_v : 0.8            Span efficiency of vertical tail
          l_{fuse} : 39        [m]  Length of fuselage    
            x_{CG} : 18        [m]  x-location of CG      
           y_{eng} : 4.83      [m]  Engine moment arm     
                                    
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
          N_{lift} : 0.2373   Wing loading multiplier
                 g : 0.2315   Gravitational acceleration
        \rho_{cap} : 0.1973   Density of spar cap material
         f_{w,add} : 0.06613  Wing added weight fraction
               r_h : 0.03417  Fractional wing thickness at spar web
        \rho_{web} : 0.03417  Density of shear web material
\sigma_{max,shear} : -0.03417 Allowable shear stress
      \sigma_{max} : -0.2031  Allowable tensile stress
                              
        V_{\infty} : 1.525    Cruise velocity       
           y_{eng} : 1.215    Engine moment arm     
               T_e : 1.182    Thrust per engine at takeoff
            \rho_c : 0.7563   Air density (35,000ft)
            V_{ne} : 0.4745   Never exceed velocity 
      C_{L_{vmax}} : 0.2373   Max lift coefficient  
           d_{fan} : 0.06617  Fan diameter          
        C_{D_{wm}} : 0.03308  Windmill drag coefficient
               \mu : 0.01229  Dynamic viscosity (35,000 ft)
\tan(\Lambda_{vt}) : -0.1122  Tangent of leading edge sweep (40 deg)
               e_v : -0.1126  Span efficiency of vertical tail
         \rho_{TO} : -0.9451  Air density (SL))     
        c_{l_{vt}} : -1.103   Sectional lift force coefficient (engine out)
          l_{fuse} : -1.427   Length of fuselage    
               V_1 : -2.365   Minimum takeoff velocity

