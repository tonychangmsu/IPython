#Title: BIOE 591 Homework 9
#Author: Tony Chang
#Abstract: Assignment 9 for ecological modeling covering land use and climate

#7. LAND USE AND CLIMATE
#
#Q: Suppose that 20% of the land area of Earth is deforested and the area subsequently desertifies. By about how much would Earth's average surface temperature change?

'''
Variables:
a: albedo of Earth (entire)
a = solar flux reflected from Earth to space/ solar flux incident on Earth
R_s: albedo of Earth (surface)
R_s = solar flux reflected from Earth's surface to the atmosphere/ solar flux incident on Earth's surface
R_a: albedo of Earth (atmosphere)
A_s: absorption of surface
A_a: absorption of atmosphere

T_a: transmission coefficient of atmosphere
T_a = 1 - R_a - A_a

a = R_a + (((1 - R_a -A_a)**2) * R_s)/(1 - (R_a * R_s))

f: total fraction of the incoming flux that is absorbed in the atmosphere

f = A_a * (1 + ((1 - R_a - A_a)*R_s)/(1 - R_a*R_s))
'''

a = 0.3
f = 86/343 

#Direct measurements of the other variables

A_a = 0.23
R_s = 0.16
R_a = 0.26

R_s_forest = 0.15
R_s_desert = 0.25

#total surface area of the earth that is forest is 6%
F_s = 0.06
delta_R_s = F_s*(R_s_desert-R_s_forest)
delta_a = 0.26*delta_R_s

a_new = a + delta_a

#using equations from Problem III.6
def T_n(n, T_o):
	return(((n+1)**(1/4))*T_o)
	
#determine a 2 level atmosphere energy balance
omega = 1372 #radiation from sun (W/m2)
sigma = 5.67e-8 #Stefan-Boltzman constant (J/m2 sec K4)
#a = 0.3 #Earth's albedo (%) ###(a=0.39  Hartmann 2005; Freedman and Kaufmann 2002)
a = a_new
n = 2 #number of layers represented in atmosphere
 
T_o = ((omega*(1-a))/(4*sigma))**(1/4)
T_s = T_n(n, T_o)
print("The estimated mean surface temperature of Earth based on the initial model is: %0.1f K \n"%(T_s))

#adding more parameters to account for the overestimate.
#Considerations are: 1. energy absorption by atmosphere, 2. latent heat flux, 3. Narrow band allowed to penetrate atmosphere

F_w = 20 #portion of IR emitted from the surface that is radiated directly to space
F_s = 86 #portion of the solar flux absorbed in the atmosphere
F_s = 0.55*(omega/4) - a*(omega/4)
F_e = 80 #flux of latent heat leaving Earth's surface
F_c = 17 #flux of convective heat leaving Earth's surface

W = (a*(omega/4)) + (sigma*T_o**4) + F_w - (omega/4)
W = 0 #initially ignoring the waste heat
T_o_hat = (((omega*(1-a))/(4*sigma)) - (F_w/sigma))**(1/4)
T_1_hat = (((2*sigma*(T_o_hat**4)) - (0.5*F_e) - (0.7*F_s))/sigma)**(1/4)
T_s_hat = (((2*sigma*(T_1_hat**4)) - (sigma*(T_o_hat**4)) + F_w - F_c - (0.5*F_e) - (0.3*F_s))/sigma)**(1/4)

print(' T_o: %0.1f K \n T_1: %0.1f K \n T_s: %0.1f K \n' %(T_o_hat,T_1_hat,T_s_hat))

'''
Exercise 4: Climate is affected by the $CO_2$ content of the atmosphere (see Problem III.8). By roughly what percentage will the atmospheric concentration of $CO_2$ immediately increase if the deforestation and subsequent burning of the cleared vegetation occur so rapidly that ocean uptake of $CO_2$ can be ignored? See Appendix, Section XII.2; assume that only tropical and temperate forests are cut and that for every 3 $km^2$ of tropical deforestation, there is 1 $km^2$ of temperate deforestation.
'''

#unclear question statement, does this mean that deforestation equals the npp from the ocean? If so, can we just subtract ocean npp from the total respired carbon per year? 

#total tropical forest area 
A_trf = 24.5e12 #m^2
npp_trf_by_A = 0.83 #kg(C)/m^2/yr
npp_trf = A_trf * npp_trf_by_A #kg/yr

#total temperate forest area
A_tef = 12.0e12 #m^2
npp_tef_by_A = 0.56 #kg(C)/m^2/yr
npp_tef = A_tef * npp_tef_by_A

#total ocean area
A_ocean = 332e12
npp_ocean_by_A = 0.057
npp_ocean = A_ocean * npp_ocean_by_A

#the question first requires solving how much tropical and temperate forest area would have to be deforested to equal the npp of the ocean, where tropical to temperate forest deforestation occurs at a 3:1 ratio

A_deforest = npp_ocean/(3 * npp_trf_by_A + npp_tef_by_A)

#solve for the total mass of plant burned
trop_bio_area = 18.8 #kg(c)/m2
trop_biomass = 3*A_deforest*trop_bio_area

temp_bio_area = 14.6 #kg(c)/m2
temp_biomass = A_deforest*temp_bio_area

total_biomass_lost = trop_biomass+temp_biomass

#total CO2 in atmosphere = 735e12 kg(C)
C_a = 735e12

#assume that the total water content in wood is roughly 40%, then there is only 60% dry carbon content
#if 1kg burning wood resulted in a 1.9 kg CO2 release?

CO2_burned_wood = total_biomass_lost * .6 * 1.9

percent_change = CO2_burned_wood/C_a

print('The total percent change of burning %0.1e kg(C) of deforested wood would result in a %0.1f%% increase in atmospheric CO2'%(total_biomass_lost,percent_change*100))


# import numpy as np

# area = np.array([24.5,12,12,8,15,9,8,18,24,14,2,2.5,332,.4,26.6,.6,1.4])
# npp = np.array([0.83,0.56,0.36,0.27,0.32,0.23,0.065,0.032,0.015,0.29,1.13,0.23,0.057,0.23,0.16,0.9,0.81])
# total_npp = np.sum(area*1e12*npp)
# area[0] = area[0]-(3*A_deforest/1e12)
# area[1] = area[1]-(A_deforest/1e12)
# new_total_npp = np.sum(area*1e12*npp)
# total_C_notsequestered = total_npp-new_total_npp #kg(C)

# #now determine the total amount of CO2 emitted per year
# #from appendix XIII.1
# C_emitted = (50+20+10+5.3+.2+.1+.1)*1e12
# new_C_emitted = C_emitted + total_C_notsequestered



#solution is 70, which would require an increase of 514.5e12 kg(c)%