from gpkit.shortcuts import *

class  B_Range(Model):

	"""Breguet Range Model

	Arguments
	--------------------
	TSFC_min 
		Minimum TSFC value, assumed constant
	MTOW 
		Max Take Off Weight
	W_oew 
		Operating empty weight
	LoverD_max
		Maximum Lift to Drag ratio
	M_max
		Maximum Mach number 

	"""
	def setup(self, TSFC_min=0.307, MTOW=10000, W_oew=7000, LoverD_max=15, M_max=0.78):
		#store attributes for later external use
		self.TSFC_min, self.MTOW, self.W_oew, self.LoverD_max, self.M_max = TSFC_min, MTOW, W_oew, LoverD_max, M_max
		TSFC_min = Var('TSFC_{min}', TSFC_min, "lb/lbf/hr")
		MTOW = Var('MTOW', MTOW, "lbf")
		W_oew = Var('W_{oew}', W_oew, "lbf")
		LoverD_max = Var('LoverD_{max}', LoverD_max)
		M_max = Var('M_max', M_max)
	
		#Constants
		g = Var('g', 9.81, "m/s^2","gravity")
		a0 = Var('a0', 340.29, "m/s", "speed of sound at sea level")

		#Free Variables
		R = Var('R', "nautical_miles", "range")
		M = Var('M', "Mach number")
		LoverD = Var('LoverD', "life to drag ratio")
		TSFC = Var('TSFC', "lb/lbf/hr", "thrust specific fuel consuption")
		W_init = Var('W_init', "lbf", "initial weight")
		W_fuel = Var('W_fuel', "lbf", "fuel weight")
		z_bre = Var('z_{bre}', "Breguet parameter")
		t = Var('t', "hr", "time")

		#Constraint Equations

		weight_eq = (W_init >= W_end + W_fuel)
		MTOW_eq = (W_init <= MTOW)
		LD_eq = (LoverD <= LoverD_max)
		TSFC_eq = (TSFC >= TSFC_min)
		M_eq = (M <= M_max)
		t_eq = (t <= T/M/a0)
		zbre_eq = (z_bre >= t*TSFC*g/LoverD)
		breguet_eq = (W_fuel/W_end >= te_exp_minus1(z_bre, nterm=3))

		#Maximum time
		return 1/t, [weight_eq, MTOW_eq, LD_eq, TSFC_eq, M_eq, t_eq, zbre_eq, breguet_eq]
