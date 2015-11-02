from gpkit.shortcuts import *
from gpkit import Variable, Model, units
import numpy.testing as npt

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

                #Set up Model Equations
		objective = 1/t #Maximize time
		constraints = [
				W_init >= W_oew + W_fuel,
				W_init <= MTOW,
				LoverD <= LoverD_max,
				TSFC >= TSFC_min,
				M <= M_max,
				t <= T/M/a0,
				z_bre >= t*TSFC*g/LoverD,
				W_fuel/W_end >= te_exp_minus1(z_bre, nterm=3)
				]
		return objective, constraints
							
		if __name__ == "__main__":
			m = Breguet_Range(0.307, 10000, 7000, 15, 0.78)
			sol = m.solve() 
