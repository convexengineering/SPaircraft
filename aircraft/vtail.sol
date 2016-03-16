Beginning signomial solve.
Solving took 4 GP solves and 0.49 seconds.

Cost
----
 3770 [N] 

Free Variables
--------------
         VerticalTail |                                                                  
               A_{vt} : 1.255             Vertical tail aspect ratio                     
          C_{D_{vis}} : 0.005139          Viscous drag coefficient                       
           C_{L_{vt}} : 0.4316            Vertical tail lift coefficient                 
              D_{vis} : 1965       [N]    Vertical tail viscous drag, cruise             
               D_{wm} : 3112       [N]    Engine out windmill drag                       
              L_{max} : 1.214e+06  [N]    Maximum wing load                              
               L_{vt} : 4.105e+04  [N]    Vertical tail lift in engine out               
                 Re_c : 3.816e+07         Vertical tail reynolds number, cruise          
                    S : 73.51      [m**2] Reference area                                 
               S_{vt} : 36.75      [m**2] Vertical tail ref. area (half span)            
           W_{struct} : 7220       [N]    Full span weight                               
               W_{vt} : 3610       [N]    Vertical tail weight                           
      \Delta x_{lead} : 12.48      [m]    Distance from CG to vertical tail leading edge 
     \Delta x_{trail} : 21         [m]    Distance from CG to vertical tail trailing edge
              \bar{c} : 6.008      [m]    Vertical tail mean aero chord                  
              \lambda : 0.27              Vertical tail taper ratio                      
                 \tau : 0.15              Thickness to chord ratio                       
                    b : 13.58      [m]    Span                                           
               b_{vt} : 6.791      [m]    Vertical tail half span                        
             c_{root} : 8.523      [m]    Vertical tail root chord                       
              c_{tip} : 2.301      [m]    Vertical tail tip chord                        
               l_{vt} : 15.55      [m]    Vertical tail moment arm                       
                    p : 1.54              Substituted variable = 1 + 2*taper             
                    q : 1.27              Substituted variable = 1 + taper               
          z_{\bar{c}} : 1.867      [m]    Vertical location of mean aerodynamic chord    
                                                                                         
WingBox, VerticalTail |                                                                  
                    A : 2.51              Aspect ratio                                   
              I_{cap} : 1.614e-06         Non-dim spar cap area moment of inertia        
                  M_r : 1.954e+05  [N]    Root moment per root chord                     
              W_{cap} : 3986       [N]    Weight of spar caps                            
              W_{web} : 1170       [N]    Weight of shear web                            
                  \nu : 0.8327            Dummy variable = $(t^2 + t + 1)/(t+1)$         
              t_{cap} : 0.0003408         Non-dim. spar cap thickness                    
              t_{web} : 0.0004447         Non-dim. shear web thickness                   

Constants
---------
         VerticalTail |                                                                   
              A_{eng} : 2.405     [m**2]     Engine reference area                        
           C_{D_{wm}} : 0.5                  Windmill drag coefficient                    
         C_{L_{vmax}} : 2.6                  Max lift coefficient                         
             L_{fuse} : 39        [m]        Length of fuselage                           
                  T_e : 1.29e+05  [N]        Thrust per engine at takeoff                 
                  V_1 : 65        [m/s]      Minimum takeoff velocity                     
                  V_c : 234       [m/s]      Cruise velocity                              
               V_{ne} : 144       [m/s]      Never exceed velocity                        
                  \mu : 1.4e-05   [N*s/m**2] Dynamic viscosity (35,000ft)                 
               \rho_c : 0.38      [kg/m**3]  Air density (35,000ft)                       
            \rho_{TO} : 1.225     [kg/m**3]  Air density (SL))                            
   \tan(\Lambda_{LE}) : 0.8391               Tangent of leading edge sweep (40 deg)       
           c_{l_{vt}} : 0.5                  Sectional lift force coefficient (engine out)
                    e : 0.8                  Span efficiency of vertical tail             
                  l_e : 4.83      [m]        Engine moment arm                            
               x_{CG} : 18        [m]        x-location of CG                             
                                                                                          
WingBox, VerticalTail |                                                                   
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
WingBox, VerticalTail |                                                      
             N_{lift} : 0.4806  Wing loading multiplier                      
                    g : 0.4788  Gravitational acceleration                   
           \rho_{cap} : 0.3701  Density of spar cap material                 
            f_{w,add} : 0.1368  Wing added weight fraction                   
                  r_h : 0.1087  Fractional wing thickness at spar web        
           \rho_{web} : 0.1087  Density of shear web material                
   \sigma_{max,shear} : -0.1087 Allowable shear stress                       
         \sigma_{max} : -0.372  Allowable tensile stress                     
                                                                             
         VerticalTail |                                                      
                  l_e : 1.497   Engine moment arm                            
                  T_e : 1.462   Thrust per engine at takeoff                 
                  V_c : 1.034   Cruise velocity                              
               V_{ne} : 0.9612  Never exceed velocity                        
               \rho_c : 0.5123  Air density (35,000ft)                       
         C_{L_{vmax}} : 0.4806  Max lift coefficient                         
           C_{D_{wm}} : 0.03526 Windmill drag coefficient                    
              A_{eng} : 0.03526 Engine reference area                        
   \tan(\Lambda_{LE}) : -0.1509 Tangent of leading edge sweep (40 deg)       
                    e : -0.2048 Span efficiency of vertical tail             
            \rho_{TO} : -0.9809 Air density (SL))                            
           c_{l_{vt}} : -1.292  Sectional lift force coefficient (engine out)
             L_{fuse} : -2.021  Length of fuselage                           
                  V_1 : -2.923  Minimum takeoff velocity                     

