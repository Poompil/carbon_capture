import csv, numpy
import math
from scipy.optimize import minimize

global RH_Hi_Matrix
global RH_Lo_Matrix
global TempHi_Matrix 
global TempLo_Matrix 
global cost_of_electricity 
global EER 
global fan_efficiency 
global capital_cost_per_frontal_area 
global packing_and_fluid_distr_cost 
global annual_operation_fraction 
global capital_charge_factor 
global maintenance_and_operation 
global packing_efficiency 
global mass_density_of_CO2_400PPMV 
global specific_packing_area 
global mass_transfer_coefficient 
global Air_Density_Matrix 
global HeatCapacity_Matrix 

RH_Hi_Matrix = numpy.genfromtxt('Mecca_RH_Hi.csv', delimiter=',')
RH_Lo_Matrix = numpy.genfromtxt('Mecca_RH_Lo.csv', delimiter=',')
TempHi_Matrix = numpy.genfromtxt('Mecca_TempHi.csv', delimiter=',')
TempLo_Matrix = numpy.genfromtxt('Mecca_TempLo.csv', delimiter=',')
cost_of_electricity = 1.39E-08	#$/J MECCA
#cost_of_electricity = 3.89E-09	#$/J TEHRAN
#cost_of_electricity = 3.35E-08	#$/J EL PASO http://www.epelectric.com/files/html/Rates/NM_Residential/Residential_Rate_No_1.pdf
#cost_of_electricity = 1.83E-08	#$/J WINDHOEK
EER = 9

fan_efficiency = 0.56
capital_cost_per_frontal_area = 3700.0 #$/m^2
packing_and_fluid_distr_cost = 250.0 #$/m^3
annual_operation_fraction = 2.70E+07 #s/yr
capital_charge_factor = 0.15
maintenance_and_operation = 0.05
packing_efficiency = 0.8
mass_density_of_CO2_400PPMV = 7.30E-4 #kg/m^3
specific_packing_area = 210 #m^2/m^3
mass_transfer_coefficient = 1.50E-03 #m/s

Air_Density_Matrix = numpy.genfromtxt('Temp_vs_density_air.csv', delimiter=',')
HeatCapacity_Matrix = numpy.genfromtxt('Temp_vs_heatcapacity.csv', delimiter=',')

total_cost_array = numpy.zeros(365)
fraction_of_the_day_requiring_cooling_array = numpy.zeros(365)

depth = 5.0	#m
velocity = 3.06618752	#m/s

parameters_to_optimize = [depth, velocity] 
constraints = [(5,20), (1,5)]

def total_cost_average_function(parameters_to_optimize):

	depth = parameters_to_optimize[0]
	velocity = parameters_to_optimize[1]

	total_cost_array = numpy.zeros(365)
	fraction_of_the_day_requiring_cooling_array = numpy.zeros(365)
	electricity_cost_saved_day = 0.0
	electricity_cost_saved_night = 0.0

	for i in range(len(RH_Hi_Matrix)):
		fraction_of_the_day_requiring_cooling = 0
		electricity_cost_saved_day = 0
		electricity_cost_saved_night = 0
		#Daytime
		inlet_dry_bulb_day = TempHi_Matrix[i]
		inlet_RH_day = RH_Lo_Matrix[i]
		
		#Nightime	
		inlet_dry_bulb_night = TempLo_Matrix[i]	
		inlet_RH_night = RH_Hi_Matrix[i]
		
		if inlet_dry_bulb_day > 27.5:
			fraction_of_the_day_requiring_cooling += 0.5
		if inlet_dry_bulb_night > 27.5:
			fraction_of_the_day_requiring_cooling += 0.5
		


		#Daytime
		if fraction_of_the_day_requiring_cooling > 0:
			
			saturation_efficiency = 1-0.024488*velocity #Based on CelDEK relationships: https://www.munters.com/en/munters/products/coolers--humidifiers/celdek-7060-15-evaporative-cooling-pad/
			
			inlet_wet_bulb_day = (inlet_dry_bulb_day*math.atan(0.151977*((inlet_RH_day+8.313659)**(0.5)))) + (math.atan(inlet_dry_bulb_day+inlet_RH_day)) - (math.atan(inlet_RH_day-1.676331)) + (0.00391838*((inlet_RH_day)**(1.5))*math.atan(0.023101*inlet_RH_day)) - 4.6686035
			outlet_dry_bulb_day_NaOH_Solution = inlet_dry_bulb_day-(saturation_efficiency*0.89*(inlet_dry_bulb_day-inlet_wet_bulb_day))
			
			density_of_air_day = numpy.interp(outlet_dry_bulb_day_NaOH_Solution, Air_Density_Matrix[:,0], Air_Density_Matrix[:,1])
			heat_capacity_day = numpy.interp(outlet_dry_bulb_day_NaOH_Solution, HeatCapacity_Matrix[:,0], HeatCapacity_Matrix[:,1])
			
			sensible_heat_energy_day = density_of_air_day*velocity*heat_capacity_day*(inlet_dry_bulb_day-outlet_dry_bulb_day_NaOH_Solution)*annual_operation_fraction
			coefficient_of_performance = EER/3.41214
			exergy_correction = math.tanh(0.147221949*(inlet_dry_bulb_day-outlet_dry_bulb_day_NaOH_Solution))
			electric_energy_saved_day = sensible_heat_energy_day/coefficient_of_performance*exergy_correction
			electricity_cost_saved_day = electric_energy_saved_day*0.5*cost_of_electricity

		#Nightime	
		if fraction_of_the_day_requiring_cooling > 0.5:

			saturation_efficiency = 1-0.024488*velocity
			
			inlet_wet_bulb_night = (inlet_dry_bulb_night*math.atan(0.151977*((inlet_RH_night+8.313659)**(0.5)))) + (math.atan(inlet_dry_bulb_night+inlet_RH_night)) - (math.atan(inlet_RH_night-1.676331)) + (0.00391838*((inlet_RH_night)**(1.5))*math.atan(0.023101*inlet_RH_night)) - 4.6686035
			outlet_dry_bulb_night_NaOH_Solution = inlet_dry_bulb_night-(saturation_efficiency*0.89*(inlet_dry_bulb_night-inlet_wet_bulb_night))
			
			density_of_air_night = numpy.interp(outlet_dry_bulb_night_NaOH_Solution, Air_Density_Matrix[:,0], Air_Density_Matrix[:,1])
			heat_capacity_night = numpy.interp(outlet_dry_bulb_night_NaOH_Solution, HeatCapacity_Matrix[:,0], HeatCapacity_Matrix[:,1])
			
			sensible_heat_energy_night = density_of_air_night*velocity*heat_capacity_night*(inlet_dry_bulb_night-outlet_dry_bulb_night_NaOH_Solution)*annual_operation_fraction
			coefficient_of_performance = EER/3.41214
			exergy_correction = math.tanh(0.147221949*(inlet_dry_bulb_night-outlet_dry_bulb_night_NaOH_Solution))
			electric_energy_saved_night = sensible_heat_energy_night/coefficient_of_performance*exergy_correction
			electricity_cost_saved_night = electric_energy_saved_night*0.5*cost_of_electricity

		
		CO2_capture_flux = annual_operation_fraction*mass_density_of_CO2_400PPMV*velocity*(1-math.exp(-packing_efficiency*specific_packing_area*depth*mass_transfer_coefficient/velocity))
		pressure_drop_in_packing = 7.4*depth*(velocity**2.14)
		energy_for_fans = annual_operation_fraction*pressure_drop_in_packing*velocity/fan_efficiency
		capital_cost = capital_cost_per_frontal_area + packing_and_fluid_distr_cost*depth
		operating_cost = energy_for_fans*cost_of_electricity + maintenance_and_operation*capital_cost - electricity_cost_saved_day - electricity_cost_saved_night
		total_cost = (operating_cost + capital_charge_factor*capital_cost)/CO2_capture_flux*1000

		total_cost_array[i] = total_cost
		fraction_of_the_day_requiring_cooling_array[i] = fraction_of_the_day_requiring_cooling
		#print i+1, ",", total_cost
		#print "Day ", i+1
		#print "Temp Day: ", inlet_dry_bulb_day, "	Temp Night: " , inlet_dry_bulb_night
		#print "RH Day: ", inlet_RH_day, "	RH Night: " , inlet_RH_night
		#print "Fraction of this day cooling was used: ", fraction_of_the_day_requiring_cooling
		#print "Cost of CO2 capture today: ", total_cost, "$/ton"
		#print "--------------------------------------------------------------"

	print "Average cost of CO2 capture: " , numpy.average(total_cost_array)
	print "Average fraction of cooling required: " , numpy.average(fraction_of_the_day_requiring_cooling_array)
	return numpy.average(total_cost_array)

#optimization methods that work with bounds:
#'L-BFGS-B'
#'TNC'
#'SLSQP'
optimizied_values = minimize(total_cost_average_function, parameters_to_optimize, bounds = constraints, method='L-BFGS-B')
average_cost_at_optimized_values = total_cost_average_function(parameters_to_optimize)
print optimizied_values
print average_cost_at_optimized_values
	
