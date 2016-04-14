Beginning signomial solve.
Solving took 5 GP solves and 4.43 seconds.

Cost
----
 5519 [N] 

Free Variables
--------------
         A_{fan} : 2.405      [m**2] Engine reference area                          
          A_{vt} : 1.095             Vertical tail aspect ratio                     
     C_{D_{vis}} : 0.00513           Viscous drag coefficient                       
      C_{L_{vt}} : 0.4231            Vertical tail lift coefficient                 
          D_{vt} : 2114       [N]    Vertical tail viscous drag, cruise             
          D_{wm} : 3112       [N]    Engine out windmill drag                       
         L_{max} : 2.615e+06  [N]    Maximum load for structural sizing             
         L_{max} : 2.615e+06  [N]    Maximum wing load                              
     L_{v_{max}} : 1.308e+06  [N]    Maximum load for structural sizing             
          L_{vt} : 4.336e+04  [N]    Vertical tail lift in engine out               
            Re_c : 4.241e+07         Vertical tail reynolds number, cruise          
               S : 79.2       [m**2] Reference area                                 
               S : 79.2       [m**2] Vertical tail reference area (full)            
          S_{vt} : 39.6       [m**2] Vertical tail ref. area (half)                 
      W_{struct} : 1.362e+04  [N]    Full span weight                               
          W_{vt} : 6812       [N]    Vertical tail weight                           
 \Delta x_{lead} : 11.53      [m]    Distance from CG to vertical tail leading edge 
\Delta x_{trail} : 21         [m]    Distance from CG to vertical tail trailing edge
         \bar{c} : 6.677      [m]    Vertical tail mean aero chord                  
    \lambda_{vt} : 0.27              Vertical tail taper ratio                      
            \tau : 0.15              Vertical tail thickness/chord ratio            
            \tau : 0.15              Thickness to chord ratio                       
               b : 13.17      [m]    Vertical tail full span                        
               b : 13.17      [m]    Span                                           
          b_{vt} : 6.584      [m]    Vertical tail half span                        
        c_{root} : 9.472      [m]    Vertical tail root chord                       
         c_{tip} : 2.557      [m]    Vertical tail tip chord                        
          l_{vt} : 14.72      [m]    Vertical tail moment arm                       
               p : 1.54              Substituted variable = 1 + 2*taper             
               q : 1.27              Substituted variable = 1 + taper               
     z_{\bar{c}} : 1.81       [m]    Vertical location of mean aerodynamic chord    
                                                                                    
         WingBox |                                                                  
               A : 2.189             Aspect ratio                                   
         I_{cap} : 2.457e-06         Non-dim spar cap area moment of inertia        
         L_{max} : 2.615e+06  [N]    Maximum wing load                              
             M_r : 3.674e+05  [N]    Root moment per root chord                     
               S : 79.2       [m**2] Reference area                                 
         W_{cap} : 7285       [N]    Weight of spar caps                            
      W_{struct} : 1.362e+04  [N]    Structural weight                              
         W_{web} : 2445       [N]    Weight of shear web                            
             \nu : 0.8327            Dummy variable = $(t^2 + t + 1)/(t+1)$         
            \tau : 0.15              Thickness to chord ratio                       
               b : 13.17      [m]    Span                                           
               p : 1.54              Substituted variable = 1 + 2*taper             
               q : 1.27              Substituted variable = 1 + taper               
         t_{cap} : 0.00052           Non-dim. spar cap thickness                    
         t_{web} : 0.0007758         Non-dim. shear web thickness                   

Constants
---------
        C_{D_{wm}} : 0.5                  Windmill drag coefficient                    
      C_{L_{vmax}} : 2.6                  Max lift coefficient                         
               T_e : 1.29e+05  [N]        Thrust per engine at takeoff                 
               V_1 : 65        [m/s]      Minimum takeoff velocity                     
               V_c : 234       [m/s]      Cruise velocity                              
            V_{ne} : 144       [m/s]      Never exceed velocity                        
               \mu : 1.4e-05   [N*s/m**2] Dynamic viscosity (35,000ft)                 
            \rho_c : 0.38      [kg/m**3]  Air density (35,000ft)                       
         \rho_{TO} : 1.225     [kg/m**3]  Air density (SL))                            
\tan(\Lambda_{LE}) : 0.8391               Tangent of leading edge sweep (40 deg)       
        c_{l_{vt}} : 0.5                  Sectional lift force coefficient (engine out)
           d_{fan} : 1.75      [m]        Fan diameter                                 
                 e : 0.8                  Span efficiency of vertical tail             
          l_{fuse} : 39        [m]        Length of fuselage                           
            x_{CG} : 18        [m]        x-location of CG                             
           y_{eng} : 4.83      [m]        Engine moment arm                            
                                                                                       
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
          N_{lift} : 0.6206  Wing loading multiplier                      
                 g : 0.6171  Gravitational acceleration                   
        \rho_{cap} : 0.462   Density of spar cap material                 
         f_{w,add} : 0.1763  Wing added weight fraction                   
               r_h : 0.1551  Fractional wing thickness at spar web        
        \rho_{web} : 0.1551  Density of shear web material                
\sigma_{max,shear} : -0.1551 Allowable shear stress                       
      \sigma_{max} : -0.4655 Allowable tensile stress                     
                                                                          
            x_{CG} : 2.368   x-location of CG                             
           y_{eng} : 1.66    Engine moment arm                            
               T_e : 1.62    Thrust per engine at takeoff                 
            V_{ne} : 1.241   Never exceed velocity                        
               V_c : 0.7599  Cruise velocity                              
      C_{L_{vmax}} : 0.6206  Max lift coefficient                         
            \rho_c : 0.377   Air density (35,000ft)                       
           d_{fan} : 0.07819 Fan diameter                                 
        C_{D_{wm}} : 0.03909 Windmill drag coefficient                    
\tan(\Lambda_{LE}) : -0.1713 Tangent of leading edge sweep (40 deg)       
                 e : -0.2552 Span efficiency of vertical tail             
         \rho_{TO} : -0.9999 Air density (SL))                            
        c_{l_{vt}} : -1.404  Sectional lift force coefficient (engine out)
          l_{fuse} : -2.368  Length of fuselage                           
               V_1 : -3.241  Minimum takeoff velocity                     

