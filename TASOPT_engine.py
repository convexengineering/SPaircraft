"SP Implementation of the TASOPT engine model"
import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, units, Vectorize, SignomialEquality
##from gpkit.nomials import SignomialEquality
from gpkit.constraints.tight import Tight as TCS
from collections import defaultdict
from gpkit.small_scripts import mag

#Cp and gamma values estimated from https://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html

#goption sets the gamma value
#Cp and gamma values estimated from https://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html
goption = 1

if goption == 1:
    fgamma = 1.401
    lpcgamma = 1.398
    hpcgamma = 1.354
    ccgamma = 1.313    #gamma value of air @ 1400 K
    lptgamma = 1.354
    hptgamma = 1.318
if goption == 0:
    fgamma = 1.401
    lpcgamma = 1.398
    hpcgamma = 1.354
    ccgamma = 1.313    #gamma value of air @ 1400 K
    lptgamma = 1.3060
    hptgamma = 1.2987

#Fan
faneta = .8948
fexp1 = (fgamma - 1)/(faneta * fgamma)

#LPC
LPCeta = .88
lpcexp1 = (lpcgamma - 1)/(LPCeta * lpcgamma)

#HPC
HPCeta = .87
hpcexp1 = (hpcgamma - 1)/(HPCeta * hpcgamma)

#combustor cooling exponents
ccexp1 = ccgamma/(1 - ccgamma)
ccexp2 = -ccgamma/(1 - ccgamma)

#Turbines
#LPT
LPTeta = .889
lptexp1 = lptgamma * LPTeta / (lptgamma - 1)

#HPT
HPTeta = .899
hptexp1 = hptgamma * HPTeta / (hptgamma - 1)

#Exhaust and Thrust
#station 8, fan exit
sta8gamma = 1.4
fanexexp = (sta8gamma - 1)/ sta8gamma

#station 6, turbine exit
sta6gamma = 1.387
turbexexp = (sta6gamma - 1) / sta6gamma

class Engine(Model):
    """
    Tasopt engine model
    """
    def setup(self, res7, cooling, N, state):
        """
        setup method for the engine model

        Arguments
        ---------
        res7 - 0 = Thrust constrained engine, 1 = burner exit temp/turbine entry temp constrained engine
        cooling - True = use cooling flow and mixing model, False = neglect cooling flow and mixing
        N - the number of flight segments
        """
        self.compressor = Compressor()
        self.combustor = Combustor()
        self.turbine = Turbine()
        self.fanmap = FanMap()
        self.lpcmap = LPCMap()
        self.hpcmap = HPCMap()
        self.thrust = Thrust()
        self.sizing = Sizing()
        self.constants = EngineConstants()
        self.state = state

        with Vectorize(N):
            #-------------------Specified Thrust or Tt4-----------------------
            #variables for the thrust constraint
            Fspec = Variable('F_{spec}', 'N', 'Specified Total Thrust')
##            Tt4spec = Variable('T_{t_{4spec}}', 'K', 'Specified Combustor Exit (Station 4) Stagnation Temperature')
##            self.state= TestState()
            self.engineP = self.dynamic(self.state, res7)
            
        models = [self.compressor , self. combustor, self. turbine, self. thrust, self.fanmap, self.lpcmap, self.hpcmap, self.sizing, self.state, self.engineP]

        #declare variables
        #max variables for the engine weight model
##        mtotmax = Variable('\dot{m}_{tot_{max}}', 'kg/s', 'Max Fan Mass Flow')
##        mCoremax = Variable('\dot{m}_{core_{max}}', 'kg/s', 'Max Core Mass Flow')
##        pifmax = Variable('\\pi_{f_{max}}', '-', 'Max Fan Pressure Ratio')
##        pilcmax = Variable('\\pi_{lc_{max}}', '-', 'Max LPC Pressure Ratio')
##        pihcmax = Variable('\\pi_{hc_{max}}', '-', 'Max HPC Pressure Ratio')

        #engine weight
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')

        #make the constraints
        constraints = []

        with SignomialsEnabled():

            weight = [
##                W_engine >= ((mtotmax/(self.constants['\\alpha_{+1_{max}}']*mCoremax)*mCoremax)*.0984)*(1684.5+17.7*(pifmax*pilcmax*pihcmax)/30+1662.2*(self.constants['\\alpha_{max}']/5)**1.2)*units('m/s'),
                W_engine >= ((self.engineP['m_{total}']/(self.engineP['alphap1']*self.engineP['m_{core}'])*self.engineP['m_{core}'])*.0984)*(1684.5+17.7*(self.engineP['\pi_f']*self.engineP['\pi_{lc}']*self.engineP['\pi_{hc}'])/30+1662.2*(self.engineP['\\alpha']/5)**1.2)*self.engineP['dum2'],

##                mtotmax >= self.engineP['m_{total}'],
##                mCoremax >= self.engineP['m_{core}'],
                  ]

            fmix = [
                #compute f with mixing
                TCS([self.combustor['eta_{B}'] * self.engineP['f'] * self.combustor['h_f'] >= (1-self.combustor['\\alpha_c'])*self.engineP['h_{t_4}']-(1-self.combustor['\\alpha_c'])*self.engineP['h_{t_3}']+self.combustor['Cp_{fuel}']*self.engineP['f']*(self.engineP['T_{t_4}']-self.combustor['T_{t_f}'])]),
                #compute Tt41...mixing causes a temperature drop
                #had to include Tt4 here to prevent it from being pushed down to zero
                SignomialEquality(self.engineP['h_{t_{4.1}}']*self.engineP['fp1'], ((1-self.combustor['\\alpha_c']+self.engineP['f'])*self.engineP['h_{t_4}'] + self.combustor['\\alpha_c']*self.engineP['h_{t_3}'])),

                self.engineP['P_{t_4}'] == self.combustor['\pi_{b}'] * self.engineP['P_{t_3}'],   #B.145
                ]

            fnomix = [
                #only if mixing = false
                #compute f without mixing, overestimation if there is cooling
                TCS([self.combustor['eta_{B}'] * self.engineP['f'] * self.combustor['h_f'] + self.engineP['h_{t_3}'] >= self.engineP['h_{t_4}']]),

                self.engineP['P_{t_4}'] == self.combustor['\pi_{b}'] * self.engineP['P_{t_3}'],   #B.145
                ]

            shaftpower = [
                #HPT shafter power balance
                #SIGNOMIAL   
                SignomialEquality(self.constants['M_{takeoff}']*self.turbine['\eta_{HPshaft}']*(1+self.engineP['f'])*(self.engineP['h_{t_{4.1}}']-self.engineP['h_{t_{4.5}}']), self.engineP['h_{t_3}'] - self.engineP['h_{t_{2.5}}']),    #B.161

                #LPT shaft power balance
                #SIGNOMIAL  
                SignomialEquality(self.constants['M_{takeoff}']*self.turbine['\eta_{LPshaft}']*(1+self.engineP['f'])*
                (self.engineP['h_{t_{4.9}}'] - self.engineP['h_{t_{4.5}}']),-((self.engineP['h_{t_{2.5}}']-self.engineP['h_{t_{1.8}}'])+self.engineP['alphap1']*(self.engineP['h_{t_2.1}'] - self.engineP['h_{t_2}']))),    #B.165
                ]

            hptexit = [
                #HPT Exit states (station 4.5)
                self.engineP['P_{t_{4.5}}'] == self.engineP['\pi_{HPT}'] * self.engineP['P_{t_{4.1}}'],
                self.engineP['\pi_{HPT}'] == (self.engineP['T_{t_{4.5}}']/self.engineP['T_{t_{4.1}}'])**(hptexp1),      #turbine efficiency is 0.9
                ]

            fanmap = [
                self.engineP['\pi_f']*(1.7/self.fanmap['\pi_{f_D}']) == (1.05*self.engineP['N_f']**.0871)**10,
                (self.engineP['\pi_f']*(1.7/self.fanmap['\pi_{f_D}']))**(.1) <= 1.1*(1.06 * (self.engineP['m_{tild_f}'])**0.137),
                (self.engineP['\pi_f']*(1.7/self.fanmap['\pi_{f_D}']))**(.1) >= .9*(1.06 * (self.engineP['m_{tild_f}'])**0.137),

                #define mbar
                self.engineP['m_{f}'] == self.engineP['m_{fan}']*((self.engineP['T_{t_2}']/self.constants['T_{ref}'])**.5)/(self.engineP['P_{t_2}']/self.constants['P_{ref}']),    #B.280

                self.engineP['\pi_f'] >= 1,

##                pifmax >= self.engineP['\pi_f'],
                ]

            lpcmap = [
                self.engineP['\pi_{lc}']*(26/self.lpcmap['\pi_{lc_D}']) == (1.38 * (self.engineP['N_1'])**0.566)**10,
                self.engineP['\pi_{lc}']*(26/self.lpcmap['\pi_{lc_D}']) <= 1.1*(1.38 * (self.engineP['m_{tild_lc}'])**0.122)**10,
                self.engineP['\pi_{lc}']*(26/self.lpcmap['\pi_{lc_D}']) >= .9*(1.38 * (self.engineP['m_{tild_lc}'])**0.122)**10,

                self.engineP['\pi_{lc}'] >= 1,

##                pilcmax >= self.engineP['\pi_{lc}'],

                #define mbar..technially not needed b/c constrained in res 2 and/or 3
                self.engineP['m_{lc}'] == self.engineP['m_{core}']*((self.engineP['T_{t_2}']/self.constants['T_{ref}'])**.5)/(self.engineP['P_{t_2}']/self.constants['P_{ref}']),    #B.280
                ]

            hpcmap = [
                self.engineP['\pi_{hc}']*(26/self.hpcmap['\pi_{hc_D}']) == (1.35 * (self.engineP['N_2'])**0.566)**10,
                self.engineP['\pi_{hc}']*(26/self.hpcmap['\pi_{hc_D}']) >= .9*(1.38 * (self.engineP['m_{tild_hc}'])**0.122)**10,
                self.engineP['\pi_{hc}']*(26/self.hpcmap['\pi_{hc_D}']) <= 1.1*(1.38 * (self.engineP['m_{tild_hc}'])**0.122)**10,

                self.engineP['m_{hc}'] == self.engineP['m_{core}']*((self.engineP['T_{t_{2.5}}']/self.constants['T_{ref}'])**.5)/(self.engineP['P_{t_{2.5}}']/self.constants['P_{ref}']),    #B.280

                self.engineP['\pi_{hc}'] >= 1,

##                pihcmax >= self.engineP['\pi_{hc}']
                ]
 
            thrust = [
                self.engineP['P_{t_8}'] == self.engineP['P_{t_7}'], #B.179
                self.engineP['T_{t_8}'] == self.engineP['T_{t_7}'], #B.180

                self.engineP['P_{t_6}'] == self.engineP['P_{t_5}'], #B.183
                self.engineP['T_{t_6}'] == self.engineP['T_{t_5}'], #B.184

                TCS([self.engineP['F_6']/(self.constants['M_{takeoff}']*self.engineP['m_{core}']) + (self.engineP['f']+1)*self.state['V'] <= (self.engineP['fp1'])*self.engineP['u_6']]),

                #ISP
                self.engineP['I_{sp}'] == self.engineP['F_{sp}']*self.state['a']*(self.engineP['alphap1'])/(self.engineP['f']*self.constants['g']),  #B.192
                ]

            res1 = [
                #residual 1 Fan/LPC speed
                self.engineP['N_f']*self.sizing['G_f'] == self.engineP['N_1'],
                self.engineP['N_1'] <= 1.1,
                self.engineP['N_2'] <= 1.1,
                ]

                #note residuals 2 and 3 differ from TASOPT, by replacing mhc with mlc
                #in residual 4 I was able to remove the LPC/HPC mass flow equality
                #in residual 6 which allows for convergence

            res2 = [
                #residual 2 HPT mass flow
                self.sizing['m_{htD}'] == (self.engineP['fp1'])*self.engineP['m_{hc}']*self.constants['M_{takeoff}']*
                (self.engineP['P_{t_{2.5}}']/self.engineP['P_{t_{4.1}}'])*
                (self.engineP['T_{t_{4.1}}']/self.engineP['T_{t_{2.5}}'])**.5,
                ]

            res3 = [
                #residual 3 LPT mass flow
                (self.engineP['fp1'])*self.engineP['m_{lc}']*self.constants['M_{takeoff}']*
                (self.engineP['P_{t_{1.8}}']/self.engineP['P_{t_{4.5}}'])*
                (self.engineP['T_{t_{4.5}}']/self.engineP['T_{t_{1.8}}'])**.5
                == self.sizing['m_{ltD}'],
                ]

            res4 = [
                #residual 4
                (self.engineP['P_{7}']/self.engineP['P_{t_7}']) == (self.engineP['T_{7}']/self.engineP['T_{t_7}'])**(3.5),
                (self.engineP['T_{7}']/self.engineP['T_{t_7}'])**-1 >= 1 + .2 * self.engineP['M_7']**2,
                ]

            res5 = [
                #residual 5 core nozzle mass flow
                (self.engineP['P_{5}']/self.engineP['P_{t_5}']) == (self.engineP['T_{5}']/self.engineP['T_{t_5}'])**(3.583979),
                (self.engineP['T_{5}']/self.engineP['T_{t_5}'])**-1 >= 1 + .2 * self.engineP['M_5']**2,
                ]


            massflux = [
                #compute core mass flux
                self.constants['M_{takeoff}'] * self.engineP['m_{core}'] == self.engineP['\\rho_5'] * self.sizing['A_5'] * self.engineP['u_5']/(self.engineP['fp1']),

                #compute fan mas flow
                self.engineP['m_{fan}'] == self.engineP['\\rho_7']*self.sizing['A_7']*self.engineP['u_7'],

                self.engineP['m_{total}'] >= self.engineP['m_{fan}'] + self.engineP['m_{core}'],
                ]

            #component area sizing
            fanarea = [
                #fan area
                self.engineP['P_2'] == self.engineP['P_{t_2}']*(self.compressor['hold_{2}'])**(-3.512),
                self.engineP['T_2'] == self.engineP['T_{t_2}'] * self.compressor['hold_{2}']**-1,
                self.sizing['A_2'] == self.engineP['m_{fan}']/(self.engineP['\\rho_2']*self.engineP['u_2']),     #B.198
                ]

            HPCarea = [
                #HPC area
                self.engineP['P_{2.5}'] == self.engineP['P_{t_{2.5}}']*(self.compressor['hold_{2.5}'])**(-3.824857),
                self.engineP['T_{2.5}'] == self.engineP['T_{t_{2.5}}'] * self.compressor['hold_{2.5}']**-1,
                self.sizing['A_{2.5}'] == self.engineP['m_{core}']/(self.engineP['\\rho_2.5']*self.engineP['u_{2.5}']),     #B.203
                ]

            onDest = [
                #estimate relevant on design values
                self.sizing['m_{htD}'] <= 1.2*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1400.0/288)**.5)/(1498/101.325),
                self.sizing['m_{htD}'] >= .8*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}']*((1400.0/288)**.5)/(1498/101.325),
                self.sizing['m_{ltD}'] <= 1.2*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1144.8/288)**.5)/(788.5/101.325),
                self.sizing['m_{ltD}'] >= .8*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}']*((1144.8/288)**.5)/(788.5/101.325),
                self.lpcmap['m_{lc_D}'] >= .8*self.sizing['m_{coreD}']*((294.5/288)**.5)/(84.25/101.325),
                self.lpcmap['m_{lc_D}'] <= 1.2*self.sizing['m_{coreD}'] *((294.5/288)**.5)/(84.25/101.325),
                self.hpcmap['m_{hc_D}'] >= .8*self.sizing['m_{coreD}'] *((482.7/288)**.5)/(399.682/101.325),
                self.hpcmap['m_{hc_D}'] <= 1.2*self.sizing['m_{coreD}'] *((482.7/288)**.5)/(399.682/101.325),
                self.fanmap['\\bar{m}_{fan_{D}}'] >= .8 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}'] *((250.0/288)**.5)/(50/101.325),
                self.fanmap['\\bar{m}_{fan_{D}}'] <= 1.2 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}']* ((250.0/288)**.5)/(50/101.325),
                ]
##908
##            print Pt0, Tt0, Pt21, Tt21, Pt25, Tt25, Pt3, Tt3, Pt3, Tt41, Pt45, Tt45

        if res7 == 0:
            res7list = [
                #residual 7
                #option #1, constrain the engine's thrust
                self.engineP['F'] == Fspec,
                ]
    
        if res7 == 1:
            if cooling == True:
                res7list = [
                    #residual 7
                    #option #2 constrain the burner exit temperature
                    self.engineP['T_{t_{4.1}}'] == Tt4spec,  #B.265
                    ]
            if cooling == False:
                res7list = [
                    #residual 7
                    #option #2 constrain the burner exit temperature
                    self.engineP['T_{t_4}'] == Tt4spec,  #B.265
                    ]

        if cooling == True:
            constraints = [weight, fmix, shaftpower, hptexit, fanmap, lpcmap, hpcmap, thrust, res1, res2, res3, res4, res5, massflux, fanarea, HPCarea, onDest, res7list]

        if cooling == False:
            constraints = [weight, fmix, shaftpower, hptexit, fanmap, lpcmap, hpcmap, thrust, res1, res2, res3, res4, res5, massflux, fanarea, HPCarea, onDest, res7list]

        #NEED TO ACCOUNT FOR THE NO MIXING CASE
        
        return models, constraints

    def dynamic(self, state, res7):
        """
        creates an instance of the engine performance model
        """
        return EnginePerformance(self, state, res7)

class EnginePerformance(Model):
    """
    Engine performance model
    """
    def setup(self, engine, state, res7, **kwargs):

        #create the subcomponent performance models
        self.compP = engine.compressor.dynamic(engine.constants, state)
        self.combP = engine.combustor.dynamic(engine.constants, state)
        self.turbineP = engine.turbine.dynamic(engine.constants)
        self.thrustP = engine.thrust.dynamic(engine.constants, state)
        self.fanmapP = engine.fanmap.dynamic(engine.constants)
        self.lpcmapP = engine.lpcmap.dynamic(engine.constants)
        self.hpcmapP = engine.hpcmap.dynamic(engine.constants)
        self.sizingP = engine.sizing.dynamic(engine.constants, engine.compressor, engine.fanmap, engine.lpcmap, engine.hpcmap, state, res7)

        models = [self.compP, self.combP, self.turbineP, self.thrustP, self.fanmapP, self.lpcmapP, self.hpcmapP, self.sizingP]
    
        return models

class EngineConstants(Model):
    """
    Class of constants used in the engine model
    """
    def setup(self):
        #-----------------------air properties------------------
        #ambient
        R = Variable('R', 287, 'J/kg/K', 'R')

        #gravity
        g = Variable('g', 9.81, 'm/(s^2)', 'Gravitational Acceleration')
        
        #-------------------------reference temp and pressure--------------------
        Tref = Variable('T_{ref}', 'K', 'Reference Stagnation Temperature')
        Pref = Variable('P_{ref}', 'kPa', 'Reference Stagnation Pressure')

        #---------------------------efficiencies & takeoffs-----------------------
        Mtakeoff = Variable('M_{takeoff}', '-', '1 Minus Percent mass flow loss for de-ice, pressurization, etc.')

        #------------------By-Pass Ratio (BPR)----------------------------
##        alphap1 = Variable('alphap1', '-', '1 plus BPR')
##        alphap1max = Variable('\\alpha_{+1_{max}}', '-', '1 Plus Max BPR')
        alphamax = Variable('\\alpha_{max}', '-', 'Max BPR')

class Compressor(Model):
    """"
    Compressor model
    """
    def setup(self):
        #define new variables
        #fan, LPC, HPC
        Cp1 = Variable('Cp_{1}', 1008, 'J/kg/K', "Cp Value for Air at 350K")#gamma = 1.398
        Cp2 = Variable('Cp_{2}', 1099, 'J/kg/K', "Cp Value for Air at 800K") #gamma = 1.354

        hold25 = Variable('hold_{2.5}', '-', '1+(gamma-1)/2 * M_2.5**2')
        hold2 = Variable('hold_{2}', '-', '1+(gamma-1)/2 * M_2**2')
        c1 = Variable('c1', '-', 'Constant in Stagnation Eqn')

        #-------------------------diffuser pressure ratios--------------------------
        pid = Variable('\pi_{d}', '-', 'Diffuser Pressure Ratio')
        pifn = Variable('\pi_{fn}', '-', 'Fan Duct Pressure Loss Ratio')

        gammaAir = Variable('gamma_{air}', 1.4, '-', 'Specific Heat Ratio for Ambient Air')
        Cpair = Variable('Cp_{air}', 1003, 'J/kg/K', "Cp Value for Air at 250K")

    def dynamic(self, engine, state):
        """
        creates an instance of the compressor performance model
        """
        return CompressorPerformance(self, engine, state) 

class CompressorPerformance(Model):
    """
    combustor perfomrance constraints
    """
    def setup(self, comp, engine, state):
        self.comp = comp
        self.engine = engine
        
        #define new variables
        #--------------------------free stream stagnation states--------------------------
        Pt0 = Variable('P_{t_0}', 'kPa', 'Free Stream Stagnation Pressure')
        Tt0 = Variable('T_{t_0}', 'K', 'Free Stream Stagnation Temperature')
        ht0 = Variable('h_{t_0}', 'J/kg', 'Free Stream Stagnation Enthalpy')

        #--------------------------diffuser exit stagnation states------------------------
        Pt18 = Variable('P_{t_{1.8}}', 'kPa', 'Stagnation Pressure at the Diffuser Exit (1.8)')
        Tt18 = Variable('T_{t_{1.8}}', 'K', 'Stagnation Temperature at the Diffuser Exit (1.8)')
        ht18 = Variable('h_{t_{1.8}}', 'J/kg', 'Stagnation Enthalpy at the Diffuser Exit (1.8)')
        
        #--------------------------fan inlet (station 2) stagnation states---------------------------
        Pt2 = Variable('P_{t_2}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2)')
        Tt2 = Variable('T_{t_2}', 'K', 'Stagnation Temperature at the Fan Inlet (2)')
        ht2 = Variable('h_{t_2}', 'J/kg', 'Stagnation Enthalpy at the Fan Inlet (2)')

        #--------------------------fan exit (station 2.1) stagnation states---------------------
        Pt21 = Variable('P_{t_2.1}', 'kPa', 'Stagnation Pressure at the Fan Exit (2.1)')
        Tt21 = Variable('T_{t_2.1}', 'K', 'Stagnation Temperature at the Fan Exit (2.1)')
        ht21 = Variable('h_{t_2.1}', 'J/kg', 'Stagnation Enthalpy at the Fan Exit (2.1)')

        #-------------------------LPC exit (station 2.5) stagnation states-------------------
        Pt25 = Variable('P_{t_{2.5}}', 'kPa', 'Stagnation Pressure at the LPC Exit (2.5)')
        Tt25 = Variable('T_{t_{2.5}}', 'K', 'Stagnation Temperature at the LPC Exit (2.5)')
        ht25 = Variable('h_{t_{2.5}}', 'J/kg', 'Stagnation Enthalpy at the LPC Exit (2.5)')

        #--------------------------HPC exit stagnation states (station 3)---------------------
        Pt3 = Variable('P_{t_3}', 'kPa', 'Stagnation Pressure at the HPC Exit (3)')
        Tt3 = Variable('T_{t_3}', 'K', 'Stagnation Temperature at the HPC Exit (3)')
        ht3 = Variable('h_{t_3}', 'J/kg', 'Stagnation Enthalpy at the HPC Exit (3)')

        #---------------------------fan nozzle exit (station 7) stagnation states---------------
        Pt7 = Variable('P_{t_7}', 'kPa', 'Stagnation Pressure at the Fan Nozzle Exit (7)')
        Tt7 = Variable('T_{t_7}', 'K', 'Stagnation Temperature at the Fan Nozzle Exit (7)')
        ht7 = Variable('h_{t_7}', 'J/kg', 'Stagnation Enthalpy at the Fan Nozzle Exit (7)')

        #------------------------turbo machinery pressure ratios--------------
        pif = Variable('\pi_f', '-', 'Fan Pressure Ratio')
        pilc = Variable('\pi_{lc}', '-', 'LPC Pressure Ratio')
        pihc = Variable('\pi_{hc}', '-', 'HPC Pressure Ratio')

        diffuser = [
            #free stream stagnation values
            Pt0 == state["P_{atm}"] / (self.comp['c1'] ** -3.5), #https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
            Tt0 == state["T_{atm}"] / (self.comp['c1']) ** (-1),             #https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
            ht0 == self.comp['Cp_{air}'] * Tt0,

            #diffuser exit stagnation values (station 1.8)
            Pt18 == self.comp['\pi_{d}'] * Pt0,  #B.113
            Tt18 == Tt0,        #B.114
            ht18 == ht0,        #B.115
            ]

        fan = [
            #fan inlet constraints (station 2)
            Tt2 == Tt18,    #B.120
            ht2 == ht18,    #B.121
            Pt2 == Pt18,
                        
            #fan exit constraints (station 2.1)
            Pt21 == pif * Pt2,  #16.50
            Tt21 == Tt2 * pif ** (fexp1),   #16.50
            ht21 == self.comp['Cp_{air}'] * Tt21,   #16.50
                       
            #fan nozzle exit (station 7)
            Pt7 == self.comp['\pi_{fn}'] * Pt21,     #B.125
            Tt7 == Tt21,    #B.126
            ht7 == ht21,    #B.127
            ]

        lpc = [
            #LPC exit (station 2.5)
            Pt25 == pilc * pif * Pt2,
            Tt25 == Tt2 * (pif*pilc) ** (lpcexp1),
            ht25 == Tt25 * self.comp['Cp_{1}'],
            ]

        hpc = [
            Pt3 == pihc * Pt25,
            Tt3 == Tt25 * pihc ** (hpcexp1),
            ht3 == self.comp['Cp_{2}'] * Tt3
            ]
        
        return diffuser, fan, lpc, hpc

class Combustor(Model):
    """"
    Combustor model
    """
    def setup(self):
        #define new variables
        Cpc = Variable('Cp_c', 1216, 'J/kg/K', "Cp Value for Fuel/Air Mix in Combustor") #1400K, gamma equals 1.312
        Cpfuel = Variable('Cp_{fuel}', 2010, 'J/kg/K', 'Specific Heat Capacity of Kerosene (~Jet Fuel)')
        hf = Variable('h_f', 43.003, 'MJ/kg', 'Heat of Combustion of Jet Fuel')     #http://hypeRbook.com/facts/2003/EvelynGofman.shtml...prob need a better source
        #43.003
        #40.8

        #-------------------------diffuser pressure ratios--------------------------
        pib = Variable('\pi_{b}', '-', 'Burner Pressure Ratio')

        #---------------------------efficiencies & takeoffs-----------------------
        etaB = Variable('eta_{B}', '-', 'Burner Efficiency')

        #------------------------Variables for cooling flow model---------------------------
        #cooling flow bypass ratio
        ac = Variable('\\alpha_c', '-', 'Total Cooling Flow Bypass Ratio')
        #variables for cooling flow velocity
        ruc = Variable('r_{uc}', '-', 'User Specified Cooling Flow Velocity Ratio')

        hold4a = Variable('hold_{4a}', '-', '1+(gamma-1)/2 * M_4a**2')

        Ttf = Variable('T_{t_f}', 'K', 'Incoming Fuel Total Temperature')

    def dynamic(self, engine, state):
        """
        creates an instance of the fan map performance model
        """
        return CombustorPerformance(self, engine, state) 

class CombustorPerformance(Model):
    """
    combustor perfomrance constraints
    """
    def setup(self, combustor, engine, state, mixing = True):
        self.combustor = combustor
        self.engine = engine

        #define new variables
        #--------------------------combustor exit (station 4) stagnation states------------------
        Pt4 = Variable('P_{t_4}', 'kPa', 'Stagnation Pressure at the Combustor Exit (4)')
        ht4 = Variable('h_{t_4}', 'J/kg', 'Stagnation Enthalpy at the Combustor Exit (4)')
        Tt4 = Variable('T_{t_4}', 'K', 'Combustor Exit (Station 4) Stagnation Temperature')

        #--------------------High Pressure Turbine inlet state variables (station 4.1)-------------------------
        Pt41 = Variable('P_{t_{4.1}}', 'kPa', 'Stagnation Pressure at the Turbine Inlet (4.1)')
        Tt41 = Variable('T_{t_{4.1}}', 'K', 'Stagnation Temperature at the Turbine Inlet (4.1)')
        ht41 = Variable('h_{t_{4.1}}', 'J/kg', 'Stagnation Enthalpy at the Turbine Inlet (4.1)')
        u41 = Variable('u_{4.1}', 'm/s', 'Flow Velocity at Station 4.1')
        T41 = Variable('T_{4.1}', 'K', 'Static Temperature at the Turbine Inlet (4.1)')

        #------------------------Variables for cooling flow model---------------------------
        #define the f plus one variable, limits the number of signomials
        fp1 = Variable('fp1', '-', 'f + 1')
        #variables for station 4a
        u4a = Variable('u_{4a}', 'm/s', 'Flow Velocity at Station 4a')
        M4a = Variable('M_{4a}', '-', 'User Specified Station 4a Mach #')
        P4a = Variable('P_{4a}', 'kPa', 'Static Pressure at Station 4a (4a)')
        uc = Variable('u_c', 'm/s', 'Cooling Airflow Speed at Station 4a')

        #---------------------------fuel flow fraction f--------------------------------
        f = Variable('f', '-', 'Fuel Air Mass Flow Fraction')
        fp1 = Variable('fp1', '-', 'f + 1')

        #make the constraints
        constraints = []        

        with SignomialsEnabled():
            #combustor constraints
            constraints.extend([
                #flow through combustor
                ht4 == self.combustor['Cp_c'] * Tt4,

                #compute the station 4.1 enthalpy
                ht41 == self.combustor['Cp_c'] * Tt41,

                #making f+1 GP compatible --> needed for convergence
                SignomialEquality(fp1,f+1),

                #investigate doing this with a substitution
                M4a == .1025,
                ])
                     
            #mixing constraints
            if mixing == True:
                constraints.extend([
                    fp1*u41 == (u4a*(fp1)*self.combustor['\\alpha_c']*uc)**.5,
                    #this is a stagnation relation...need to fix it to not be signomial
                    SignomialEquality(T41, Tt41-.5*(u41**2)/self.combustor['Cp_c']),
                    
                    #here we assume no pressure loss in mixing so P41=P4a
                    Pt41 == P4a*(Tt41/T41)**(ccexp1),
                    #compute station 4a quantities, assumes a gamma value of 1.313 (air @ 1400K)
                    u4a == M4a*((1.313*self.engine['R']*Tt4)**.5)/self.combustor['hold_{4a}'],
                    uc == self.combustor['r_{uc}']*u4a,
                    P4a == Pt4*self.combustor['hold_{4a}']**(ccexp2),
                    ])
            #combustor constraints with no mixing
            else:
                constraints.extend([
                    Pt41 == Pt4,
                    Tt41 == Tt4,
                    ])
                
        return constraints

class Turbine(Model):
    """"
    Turbine model
    """
    def setup(self):
        #define new variables
         #turbines
##        Cpt1 =Variable('Cp_t1', 1190, 'J/kg/K', "Cp Value for Combustion Products in HP Turbine") #1300K gamma = 1.318
##        Cpt2 =Variable('Cp_t2', 1099, 'J/kg/K', "Cp Value for Combustion Products in LP Turbine") #800K gamma = 1.354
        Cpt1 =Variable('Cp_t1', 1280, 'J/kg/K', "Cp Value for Combustion Products in HP Turbine") #1300K gamma = 1.318
        Cpt2 =Variable('Cp_t2', 1184, 'J/kg/K', "Cp Value for Combustion Products in LP Turbine") #800K gamma = 1.354

        #-------------------------diffuser pressure ratios--------------------------
        pitn = Variable('\\pi_{tn}', '-', 'Turbine Nozzle Pressure Ratio')

        #---------------------------efficiencies & takeoffs-----------------------
        etaHPshaft = Variable('\eta_{HPshaft}', '-', 'Power Transmission Efficiency of High Pressure Shaft, Smears in Losses for Electrical Power')
        etaLPshaft = Variable('\eta_{LPshaft}', '-', 'Power Transmission Efficiency of Low Pressure Shaft, Smeras in Losses for Electrical Power')

    def dynamic(self, engine):
        """
        creates an instance of the fan map performance model
        """
        return TurbinePerformance(self, engine)

class TurbinePerformance(Model):
    """
    combustor perfomrance constraints
    """
    def setup(self, turbine, engine):
        self.turbine = turbine
        self.engine = engine

        #define new variables
        #------------------LPT inlet stagnation states (station 4.5)------------------
        ht45 = Variable('h_{t_{4.5}}', 'J/kg', 'Stagnation Enthalpy at the HPT Exit (4.5)')
        Pt45 = Variable('P_{t_{4.5}}', 'kPa', 'Stagnation Pressure at the HPT Exit (4.5)')
        Tt45 = Variable('T_{t_{4.5}}', 'K', 'Stagnation Temperature at the HPT Exit (4.5)')

        #-------------------HPT exit (station 4.5) stagnation states--------------------
        Pt49 = Variable('P_{t_{4.9}}', 'kPa', 'Stagnation Pressure at the HPTExit (49)')
        Tt49 = Variable('T_{t_{4.9}}', 'K', 'Stagnation Temperature at the HPT Exit (49)')
        ht49 = Variable('h_{t_{4.9}}', 'J/kg', 'Stagnation Enthalpy at the HPT Exit (49)')

        #------------------------turbo machinery pressure ratios--------------
        pihpt = Variable('\pi_{HPT}', '-', 'HPT Pressure Ratio')
        pilpt = Variable('\pi_{LPT}', '-', 'LPT Pressure Ratio')

        #------------------turbine nozzle exit stagnation states (station 5)------------
        Pt5 = Variable('P_{t_5}', 'kPa', 'Stagnation Pressure at the Turbine Nozzle Exit (5)')
        Tt5 = Variable('T_{t_5}', 'K', 'Stagnation Temperature at the Turbine Nozzle Exit (5)')
        ht5 = Variable('h_{t_5}', 'J/kg', 'Stagnation Enthalpy at the Turbine Nozzle Exit (5)')

        #make the constraints
        constraints = []

        with SignomialsEnabled():
            #turbine constraints
            constraints.extend([
                #HPT Exit states (station 4.5)
                ht45 == self.turbine['Cp_t1'] * Tt45,

                #LPT Exit States
                Pt49 == pilpt * Pt45,
                pilpt == (Tt49/Tt45)**(lptexp1),    #turbine efficiency is 0.9
                ht49 == self.turbine['Cp_t2'] * Tt49,

                #turbine nozzle exit states
                Pt5 == self.turbine['\\pi_{tn}'] * Pt49, #B.167
                Tt5 == Tt49,    #B.168
                ht5 == ht49     #B.169
                ])

        return constraints

class FanMap(Model):
    """"
    Fan map model
    """
    def setup(self):
        #define new variables
        #------------------Fan map variables----------------
##        mFanD = Variable('m_{fan_D}', 'kg/s', 'Fan On-Design Mass Flow')
        mFanBarD = Variable('\\bar{m}_{fan_{D}}', 'kg/s', 'Fan On-Design Corrected Mass Flow')
        piFanD = Variable('\pi_{f_D}', '-', 'On-Design Pressure Ratio')

    def dynamic(self, engine):
        """
        creates an instance of the fan map performance model
        """
        return FanMapPerformance(self, engine)   

class FanMapPerformance(Model):
    """
    Fan map perfomrance constraints
    """
    def setup(self, fanmap, engine):
        self.fanmap = fanmap
        self.engine = engine

        #define new variables
        #-----------------------Fan Map Variables--------------------
        #Mass Flow Variables
        mf = Variable('m_{f}', 'kg/s', 'Fan Corrected Mass Flow')
        mtildf = Variable('m_{tild_f}', '-', 'Fan Normalized Mass Flow')
        
        #----------------------Compressor Speeds--------------------
        #Speed Variables...by setting the design speed to be 1 since only ratios are
        #imporant I was able to drop out all the other speeds
        Nf = Variable('N_f', '-', 'Fan Speed')

        #make the constraints
        constraints = []
        
        #fan map
        constraints.extend([
            #define mtild
            mtildf == mf/self.fanmap['\\bar{m}_{fan_{D}}'],   #B.282
            ])

        return constraints

class LPCMap(Model):
    """"
    LPC map model
    """
    def setup(self):
        #define new variables
        #-----------------LPC map variables-------------------
        mlcD = Variable('m_{lc_D}', 'kg/s', 'On Design LPC Corrected Mass Flow')
        pilcD = Variable('\pi_{lc_D}', '-', 'LPC On-Design Pressure Ratio')

    def dynamic(self, engine):
        """
        creates an instance of the HPC map performance model
        """
        return LPCMapPerformance(self, engine)   
        
class LPCMapPerformance(Model):
    """
    LPC map perfomrance constraints
    """
    def setup(self, lpcmap, engine):
        self.lpcmap = lpcmap
        self.engine = engine

        #define new variables
        #-------------------------LPC Map Variables-------------------------
        #Mass Flow Variables
        mlc = Variable('m_{lc}', 'kg/s', 'LPC Corrected Mass Flow')
        mtildlc = Variable('m_{tild_lc}', '-', 'LPC Normalized Mass Flow')

        #----------------------Compressor Speeds--------------------
        #Speed Variables...by setting the design speed to be 1 since only ratios are
        #imporant I was able to drop out all the other speeds
        N1 = Variable('N_1', '-', 'LPC Speed')

        #make the constraints
        constraints = []

        #LPC map
        constraints.extend([
            #define mtild
            mtildlc == mlc/self.lpcmap['m_{lc_D}'],   #B.282
        ])

        return constraints

class HPCMap(Model):
    """"
    HPC map model
    """
    def setup(self):
        #define new variables
        #-----------------HPC map variables-----------
        mhcD = Variable('m_{hc_D}', 'kg/s', 'On Design HPC Corrected Mass Flow')
        pihcD = Variable('\pi_{hc_D}', '-', 'HPC On-Design Pressure Ratio')
        mhcD = Variable('m_{hc_D}', 'kg/s', 'On Design HPC Corrected Mass Flow')

    def dynamic(self, engine):
        """
        creates an instance of the HPC map performance model
        """
        return HPCMapPerformance(self, engine)   
        
class HPCMapPerformance(Model):
    """
    HPC map perfomrance constraints
    """
    def setup(self, hpcmap, engine):
        self.hpcmap = hpcmap
        self.engine = engine

        #define new variables
        #--------------------------HPC Map Variables------------------
        #Mass Flow Variables
        mhc = Variable('m_{hc}', 'kg/s', 'HPC Corrected Mass Flow')
        mtildhc = Variable('m_{tild_hc}', '-', 'HPC Normalized Mass Flow')

        #----------------------Compressor Speeds--------------------
        #Speed Variables...by setting the design speed to be 1 since only ratios are
        #imporant I was able to drop out all the other speeds
        N2 = Variable('N_2', '-', 'HPC Speed')

        #make the constraints
        constraints = []

        #HPC map
        constraints.extend([
            #define mtild
            mtildhc == mhc/self.hpcmap['m_{hc_D}'],   #B.282
            ])

        return constraints

class Thrust(Model):
    """"
    thrust sizing model
    """
    def setup(self):
        #define new variables
        #fan and exhaust
        Cptex =Variable('Cp_tex', 1029, 'J/kg/K', "Cp Value for Combustion Products at Core Exhaust") #500K, gamma = 1.387
        Cpfanex = Variable('Cp_fex', 1005, 'J/kg/K', "Cp Value for Air at 300K") #gamma =1.4        #heat of combustion of jet fuel

    def dynamic(self, engine, state):
        """
        creates an instance of the thrust performance model
        """
        return ThrustPerformance(self, engine, state)

class ThrustPerformance(Model):
    """
    thrust performacne model
    """
    def setup(self, thrust, engine, state):
        self.thrust = thrust
        self.engine = engine
        
        #define new variables
        #------------------fan exhaust (station 8) statge variables------------
        P8 = Variable('P_8', 'kPa', 'Fan Exhaust Static Pressure')
        Pt8 = Variable('P_{t_8}', 'kPa', 'Fan Exhaust Stagnation Pressure')
        ht8 = Variable('h_{t_8}', 'J/kg', 'Fan Exhaust Stagnation Enthalpy')
        h8 = Variable('h_8', 'J/kg', 'Fan Exhasut Static Enthalpy')
        Tt8 = Variable('T_{t_8}', 'K', 'Fan Exhaust Stagnation Temperature (8)')
        T8 = Variable('T_{8}', 'K', 'Fan Exhaust Sttic Temperature (8)')

        #-----------------core exhaust (station 6) state variables-------------
        P6 = Variable('P_6', 'kPa', 'Core Exhaust Static Pressure')
        Pt6 = Variable('P_{t_6}', 'kPa', 'Core Exhaust Stagnation Pressure')
        Tt6 = Variable('T_{t_6}', 'K', 'Core Exhaust Stagnation Temperature (6)')
        T6 = Variable('T_{6}', 'K', 'Core Exhaust Static Temperature (6)')
        ht6 = Variable('h_{t_6}', 'J/kg', 'Core Exhaust Stagnation Enthalpy')
        h6 = Variable('h_6', 'J/kg', 'Core Exhasut Static Enthalpy')

        #thrust variables
        F8 = Variable('F_8', 'N', 'Fan Thrust')
        F6 = Variable('F_6', 'N', 'Core Thrust')
        F = Variable('F', 'N', 'Total Thrust')
        Fsp = Variable('F_{sp}', '-', 'Specific Net Thrust')
        Isp = Variable('I_{sp}', 's', 'Specific Impulse')
        TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumption')

        #exhaust speeds
        u6 = Variable('u_6', 'm/s', 'Core Exhaust Velocity')
        u8 = Variable('u_8', 'm/s', 'Fan Exhaust Velocity')

        #mass flows
        mCore = Variable('m_{core}', 'kg/s', 'Core Mass Flow')
        mFan = Variable('m_{fan}', 'kg/s', 'Fan Mass Flow')

        #------------------By-Pass Ratio (BPR)----------------------------
        alpha = Variable('\\alpha', '-', 'By Pass Ratio')
        alphap1 = Variable('alphap1', '-', '1 plus BPR')

        hold = Variable('hold', '-', 'unecessary hold var')

        #constraints
        constraints = []

        with SignomialsEnabled():

            #exhaust and thrust constraints
            constraints.extend([
                P8 == state["P_{atm}"],
                h8 == self.thrust['Cp_fex'] * T8,
                TCS([u8**2 + 2*h8 <= 2*ht8]),
##               SignomialEquality(u8**2 + 2*h8, 2*ht8),
                (P8/Pt8)**(fanexexp) == T8/Tt8,
                ht8 == self.thrust['Cp_fex'] * Tt8,
                
                #core exhaust
                P6 == state["P_{atm}"],   #B.4.11 intro
 
                (P6/Pt6)**(turbexexp) == T6/Tt6,
                TCS([u6**2 + 2*h6 <= 2*ht6]),
##                SignomialEquality(u6**2 + 2*h6, 2*ht6),
                h6 == self.thrust['Cp_tex'] * T6,
                ht6 == self.thrust['Cp_tex'] * Tt6,

                u6 >= state['V'],
                u8 >= state['V'],
                F6 >= .01*units('N'),
                F8 >= .01*units('N'),

                #constrain the new BPR
                alpha == mFan / mCore,
                hold == alphap1,
                SignomialEquality(hold, alpha + 1),
                alpha <= 5.1362,

                #overall thrust values
                TCS([F8/(alpha * mCore) + state['V'] <= u8]),  #B.188
                
                #SIGNOMIAL
                TCS([F <= F6 + F8]),
##                SignomialEquality(F, F6 + F8),

                Fsp == F/((alphap1)*mCore*state['a']),   #B.191

                F >= .1*units('N'),

                #TSFC
                TSFC == 1/Isp,

##                self.engine['\\alpha_{max}'] >= alpha,
##                self.engine['\\alpha_{+1_{max}}'] >= self.engine['alphap1'],
                ])

        return constraints

class Sizing(Model):
    """"
    engine sizing model
    """
    def setup(self):
        #define new variables
        #gear ratio, set to 1 if no gearing present
        Gf = Variable('G_f', '', 'Gear Ratio Between Fan and LPC')

        mhtD = Variable('m_{htD}', 'kg/s', 'Design HPT Corrected Mass Flow (see B.225)')
        mltD = Variable('m_{ltD}', 'kg/s', 'Design LPT Corrected Mass Flow (see B.226)')

        #on design by pass ratio
        alpha_OD = Variable('\\alpha_{OD}', '-', 'By Pass Ratio')

        #-------------------------Areas------------------------
        A2 = Variable('A_2', 'm^2', 'Fan Area')
        A25 = Variable('A_{2.5}', 'm^2', 'HPC Area')
        A5 = Variable('A_5', 'm^2', 'Core Exhaust Nozzle Area')
        A7 = Variable('A_7', 'm^2', 'Fan Exhaust Nozzle Area')

        mCoreD = Variable('m_{coreD}', 'kg/s', 'Estimated on Design Mass Flow')  

        #constraints
        constraints = []

        constraints.extend([
            #-------------------------Areas------------------------
            A5 + A7 <= A2,
            ])

        return constraints

    def dynamic(self, engine, compressor, fanmap, lpcmap, hpcmap, state, res7):
        """
        creates an instance of the engine sizing performance model
        """
        return SizingPerformance(self, engine, compressor, fanmap, lpcmap, hpcmap, state, res7)

class SizingPerformance(Model):
    """
    engine sizing perofrmance model
    """
    def setup(self, sizing, engine, compressor, fanmap, lpcmap, hpcmap, state, res7, cooling = True):
        self.sizing = sizing
        self.engine = engine
        self.compressor = compressor
        self.fanmap = fanmap
        self.lpcmap = lpcmap
        self.hpcmap = hpcmap

        #new variables
        #exhaust mach numbers
        a5 = Variable('a_5', 'm/s', 'Speed of Sound at Station 5')
        a7 = Variable('a_7', 'm/s', 'Speed of Sound at Station 7')
        
        #mass flows
        mtot = Variable('m_{total}', 'kg/s', 'Total Engine Mass Flux')
        
        #-------------------fan face variables---------------------
        rho2 = Variable('\\rho_2', 'kg/m^3', 'Air Static Density at Fan Face')
        T2 = Variable('T_2', 'K', 'Air Static Temperature at Fan Face')
        P2 = Variable('P_2', 'kPa', 'Air Static Pressure at Fan Face')
        u2 = Variable('u_2', 'm/s', 'Air Speed at Fan Face')
        h2 = Variable('h_{2}', 'J/kg', 'Static Enthalpy at the Fan Inlet (2)')
        M2 = Variable('M_2', '-', 'Fan Face/LPC Face Axial Mach Number')

        #------------------HPC face variables---------------------
        rho25 = Variable('\\rho_2.5', 'kg/m^3', 'Static Air Density at HPC Face')
        T25 = Variable('T_{2.5}', 'K', 'Static Air Temperature at HPC Face')
        P25 = Variable('P_{2.5}', 'kPa', 'Static Air Pressure at HPC Face')
        u25 = Variable('u_{2.5}', 'm/s', 'Air Speed at HPC Face')
        h25 = Variable('h_{2.5}', 'J/kg', 'Static Enthalpy at the LPC Exit (2.5)')
        M25 = Variable('M_{2.5}', '-', 'HPC Face Axial Mach Number')

        #fan exhuast states
        P7 = Variable('P_{7}', 'kPa', 'Fan Exhaust Static Pressure (7)')
        T7 = Variable('T_{7}', 'K', 'Static Temperature at the Fan Nozzle Exit (7)')
        u7 = Variable('u_7', 'm/s', 'Station 7 Exhaust Velocity')
        M7 = Variable('M_7', '-', 'Station 7 Mach Number')
        rho7 = Variable('\\rho_7', 'kg/m^3', 'Air Static Density at Fam Exhaust Exit (7)')

        #core exhaust states
        P5 = Variable('P_{5}', 'kPa', 'Core Exhaust Static Pressure (5)')
        T5 = Variable('T_{5}', 'K', 'Static Temperature at the Turbine Nozzle Exit (5)')
        M5 = Variable('M_5', '-', 'Station 5 Mach Number')
        u5 = Variable('u_5', 'm/s', 'Station 5 Exhaust Velocity')
        rho5 = Variable('\\rho_5', 'kg/m^3', 'Air Static Density at Core Exhaust Exit (5)')

        #dummy vairables for pint purposes
        dum = Variable("dum", 781, 'J/kg/K')
        dum2 = Variable('dum2', 1, 'm/s')

        #constraints
        constraints = []
 
        #sizing constraints that hold the engine together
        constraints.extend([
            #residual 4
            P7 >= state["P_{atm}"],
            M7 <= 1,
            u7 >= state['V'],
            a7 == (1.4*self.engine['R']*T7)**.5,
            a7*M7 == u7,
            rho7 == P7/(self.engine['R']*T7),
            
            #residual 5 core nozzle mass flow
            P5 >= state["P_{atm}"],
            M5 <= 1,
            u5 >= state['V'],
            a5 == (1.387*self.engine['R']*T5)**.5,
            a5*M5 == u5,
            rho5 == P5/(self.engine['R']*T5),
            
            #component area sizing
            #fan area
            h2 == self.compressor['Cp_{1}'] * T2,
            rho2 == P2/(self.engine['R'] * T2),  #B.196
            u2 == M2*(self.compressor['Cp_{1}']*self.engine['R']*T2/(dum))**.5,  #B.197

            #HPC area
 
            h25 == self.compressor['Cp_{2}'] * T25,
            rho25 == P25/(self.engine['R']*T25),
            u25 == M25*(self.compressor['Cp_{2}']*self.engine['R']*T25/(dum))**.5,   #B.202
        ])

        return constraints

class TestState(Model):
    """
    state class only to be used for testing purposes
    """
    def setup(self):
        #define variables
        p_atm = Variable("P_{atm}", "kPa", "air pressure")
##        R_atm = Variable("R_{atm}", "J/mol/K", "air specific heating value")
        TH = 5.257386998354459 #(g*M_atm/R_atm/L_atm).value
##        rho = Variable('\\\rho', 'kg/m^3', 'Density of air')
        T_atm = Variable("T_{atm}", "K", "air temperature")
  

##        mu  = Variable('\\mu', 'kg/(m*s)', 'Dynamic viscosity')

        V = Variable('V', 'kts', 'Aircraft Flight Speed')
        a = Variable('a', 'm/s', 'Speed of Sound')
##        h = Variable('h', 'm', 'Segment Altitude [meters]')
##        hft = Variable('hft', 'feet', 'Segment Altitude [feet]')
        R = Variable('R', 287, 'J/kg/K', 'Air Specific Heat')
        gamma = Variable('\\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        M = Variable('M', '-', 'Mach Number')

        #make constraints
        constraints = []

        constraints.extend([
            #INVESTIGATE RHO
##            rho == .38*units('kg/m^3'),
            
            V == M * a,
            a  == (gamma * R * T_atm)**.5,
            ])

        return constraints

    def bound_all_variables(self, model, eps=1e-30, lower=None, upper=None):
        "Returns model with additional constraints bounding all free variables"
        lb = lower if lower else eps
        ub = upper if upper else 1/eps
        constraints = []
        freevks = tuple(vk for vk in model.varkeys if "value" not in vk.descr)
        for varkey in freevks:
            units = varkey.descr.get("units", 1)
            varub = Variable('varub', ub, units)
            varlb = Variable('varls', lb, units)
            constraints.append([varub >= Variable(**varkey.descr),
                                Variable(**varkey.descr) >= varlb])
        m = Model(model.cost, [constraints, model], model.substitutions)
        m.bound_all = {"lb": lb, "ub": ub, "varkeys": freevks}
        return m

    # pylint: disable=too-many-locals
    def determine_unbounded_variables(self, model, verbosity=4,
                                      eps=1e-30, lower=None, upper=None, **kwargs):
        "Returns labeled dictionary of unbounded variables."
        m = self.bound_all_variables(model, eps, lower, upper)
        sol = m.localsolve(solver='mosek', verbosity=4, iteration_limit = 100, **kwargs)
        solhold = sol
        lam = sol["sensitivities"]["la"][1:]
        out = defaultdict(list)
        for i, varkey in enumerate(m.bound_all["varkeys"]):
            lam_gt, lam_lt = lam[2*i], lam[2*i+1]
            if abs(lam_gt) >= 1e-7:  # arbitrary threshold
                out["sensitive to upper bound"].append(varkey)
            if abs(lam_lt) >= 1e-7:  # arbitrary threshold
                out["sensitive to lower bound"].append(varkey)
            value = mag(sol["variables"][varkey])
            distance_below = np.log(value/m.bound_all["lb"])
            distance_above = np.log(m.bound_all["ub"]/value)
            if distance_below <= 3:  # arbitrary threshold
                out["value near lower bound"].append(varkey)
            elif distance_above <= 3:  # arbitrary threshold
                out["value near upper bound"].append(varkey)
        return out, solhold

class TestMission(Model):
    """
    place holder of a mission calss
    """
    def setup(self, engine):
        M2 = .8025
        M25 = .6
        M4a = .1025
        Mexit = 1
        M0 = .8025

        toclimb = [
            engine['F_{spec}'][0] == 94.971 * units('kN'),
##            engine['F_{spec}'][1] == 63.184 * units('kN'),
            engine['F_{spec}'][1] == 30.109* units('kN'),
            engine['F_{spec}'][2] == 22.182 * units('kN'),

#CHECK THE TT4SPEC

##            engine['T_{t_{4spec}}'] [0]== 1400*units('K'),
##            engine['T_{t_{4spec}}'][1] == 1400*units('K'),

            engine.state['P_{atm}'][1] == 23.84*units('kPa'),    #36K feet
            engine.state["T_{atm}"][1] == 218*units('K'),
            engine.state['M'][1] == M0,
            engine.engineP['M_2'][1] == M2,
            engine.engineP['M_{2.5}'][1] == M25,
            engine.compressor['hold_{2}'] == 1+.5*(1.398-1)*M2**2,
            engine.compressor['hold_{2.5}'] == 1+.5*(1.354-1)*M25**2,
            engine.compressor['c1'] == 1+.5*(.401)*M0**2,
            ]

        M2 = .8
        M25 = .6
        M4a = .1025
        Mexit = 1
        M0 = .8

        cruise = [
            engine.state['P_{atm}'][2] == 23.84*units('kPa'),    #36K feet
            engine.state["T_{atm}"][2] == 218*units('K'),
            engine.state['M'][2] == M0,
            engine.engineP['M_2'][2] == M2,
            engine.engineP['M_{2.5}'][2] == M25,
            ]

        M2 = .223
        M25 = .6
        M4a = .1025
        Mexit = 1
        M0 = .223

        rotation = [
            engine.state['P_{atm}'][0] == 101.325*units('kPa'),    #36K feet
            engine.state["T_{atm}"][0] == 291*units('K'),
            engine.state['M'][0] == M0,
            engine.engineP['M_2'][0] == M2,
            engine.engineP['M_{2.5}'][0] == M25,
            ]

        M2 = .4596
        M25 = .6
        M4a = .1025
        Mexit = 1
        M0 = .4596

##        climb = [
##            engine.state['P_{atm}'][1] == 72.4*units('kPa'),    #36K feet
##            engine.state["T_{atm}"][1] == 270*units('K'),
##            engine.state['M'][1] == M0,
##            engine.engineP['M_2'][1] == M2,
##            engine.engineP['M_{2.5}'][1] == M25,
##            ]

        return rotation, toclimb, cruise

if __name__ == "__main__":
    engine = Engine(0, True, 3)
    mission = TestMission(engine)

    M4a = .1025
    fan = 1.685
    lpc  = 4.744
    hpc = 3.75
 
    substitutions = {
            '\\pi_{tn}': .989,
            '\pi_{b}': .94,
            '\pi_{d}': .998,
            '\pi_{fn}': .98,
            'T_{ref}': 288.15,
            'P_{ref}': 101.325,
            '\eta_{HPshaft}': .97,
            '\eta_{LPshaft}': .97,
            'eta_{B}': .9827,

            '\pi_{f_D}': fan,
            '\pi_{hc_D}': hpc,
            '\pi_{lc_D}': lpc,

##            '\\alpha_{OD}': 5.105,

##            'M_{4a}': M4a,
            'hold_{4a}': 1+.5*(1.313-1)*M4a**2,#sol('hold_{4a}'),
            'r_{uc}': .5,
            '\\alpha_c': .19036,
            'T_{t_f}': 435,

            'M_{takeoff}': .972,

            'G_f': 1,

            'h_f': 43.003,

            'Cp_t1': 1280,
            'Cp_t2': 1184,
            'Cp_c': 1216,
           }

    m = Model((10*engine.engineP.thrustP['TSFC'][2]+engine.engineP.thrustP['TSFC'][1]+engine.engineP.thrustP['TSFC'][0]) * (engine['W_{engine}'] * units('1/hr/N'))**.00001, [engine, mission], substitutions)
##    for posy in m.sp().gp().posynomials:
##            print posy.str_without(["models"])
    m.substitutions.update(substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 1)
##    ts = TestState()
##    bounds, sol = ts.determine_unbounded_variables(m)

    rotationerror = 100*(mag(sol('TSFC')[0]) - .48434)/.48434
##    climberror = 100*(mag(sol('TSFC')[1]) - 0.57297)/0.57297
    tocerror = 100*(mag(sol('TSFC')[1]) - .65290)/.65290
    cruiseerror = 100*(mag(sol('TSFC')[2]) - .64009)/.64009

    print rotationerror, tocerror, cruiseerror

    print 100*(mag(sol('A_2').to('m^2'))-1.6026)/1.6026
    print 100*(mag(sol('A_7').to('m^2'))-.7423)/.7423
    print 100*(mag(sol('A_5').to('m^2'))-.2262)/.2262


    Pt0 = 50
    Tt0 = 250
    Pt3 = Pt0*fan*lpc*hpc
    Pt21 = fan * Pt0
    Pt25 = Pt0 *fan* lpc

    Tt21 = Tt0 * (fan)**(.4/(1.4*.91))
    Tt25 = Tt21 * (lpc)**(.4/(1.4*.90))
    Tt3 = Tt25 * (hpc)**(.4/(1.4*.89))

    Tt41 = 1400

    Tt45 = Tt41 - (Tt3 - Tt25)

    Tt49 = Tt45 - (Tt25 - Tt21)

    piHPT = (Tt45/Tt41)**(.9121*1.4/.4)

    piLPT = (Tt49/Tt45)**(.9228*1.4/.4)

    Pt45 = piHPT * Pt3
    print Pt0, Tt0, Pt21, Tt21, Pt25, Tt25, Pt3, Tt3, Pt3, Tt41, Pt45, Tt45
