
#BIOE 591: Ecological Modeling Homework Assignment 1
Author: Tony Chang

Date: 09/08/14

Abstract: Solutions to BIOE1 assignment 1 written in Python 3.4

Dependencies: Python 3.4

###Pre-exercise question: *The Greens We Eat*

    What fraction of the total annual plant growth on Earth was eaten by humans
in 1983?


    #constants
    daily_cal_human = 2.5e6 # average daily human consumption in calories per day
    population = 4.7e9 # total number of humans in 1983
    J_per_cal = 4.18 # conversion of J/cal
    npp_g = 7.5e16 # g(C)/yr amount of C assimilated by plants per year
    e_per_biomass = 1.6e4 # J/g(biomass dry) amount of energy per gram dry biomass
    days_year = 365


    #we estimate biomass grown by plants with total glucose produced given known carbon assimilation
    #write a function to calculate the carbon within glucose
    def biomassGlucose(C,H,O):
    	mC = 12.011; mH = 1.0079; mO = 15.999 #mass of C, H, and O per mol
    	return((C*mC)/(C*mC+H*mH+O*mO))	


    glucose = biomassGlucose(6, 12, 6) 
    npp_J = npp_g*ec/glucose #npp in units J/year
    daily_J_human = daily_cal_human * J_per_cal
    annual_J_human = daily_J_human * days_year
    pop_J = annual_J_human * population
    f = pop_J/npp_J
    print("%.1f percent plant material (glucose) are consumed by humans annually"%(f*100))

    0.6 percent plant material (glucose) are consumed by humans annually


My intial guess of human consumption of total annual plant growth was 5%. It
turns out, that based on
the estimates with assumptions of a 2,500 kcal diet per human and known net
primary productivity that only 0.6% of plant material are consumed.

--------------

###Exercise 1
_A formula representing the approximate chemical composition of typical dry
freshly photosynthesized biomass is $H_{2960} O_{1480} C_{1480} N_{160} P_{18}
S_{10}$, where each subscript denotes the relative number of atoms of that
elemental type. If this is more precise representation is used instead of $C_6
H_{12} O_6$, recalculate the fraction, f._

Begin by writing a function for the carbon content in dry biomass and compare
that estimate for carbon content in glucose.


    #function to calculate the carbon content in dry biomass (g(C)/g(biomass))
    def biomassDry(H,O,C,N,P,S):	
    	mC = 12.011; mH = 1.0079; mO = 15.999 #mass of C, H, and O per mol
    	mS = 32.065; mP = 30.974; mN = 14.007 #mass of S, P, and N per mol
    	return((C*mC)/(C*mC+H*mH+O*mO+N*mN+P*mP+S*mS))	


    drymass = biomassDry(2960,1480,1480,160,18,10) 
    print("Carbon content in glucose: %.2f" %glucose)
    print("Carbon content in dry biomass: %.2f" %drymass)
    npp_J_ex1 = npp_g*ec/drymass #npp in units J/year
    ex1_f = pop_J/npp_J_ex1
    print("%.1f percent plant material consumed by humans annually"%(ex1_f*100))

    Carbon content in glucose: 0.40
    Carbon content in dry biomass: 0.37
    0.6 percent plant material consumed by humans annually


Using the more precise estimate of photosynthesized biomass do not make a major
change in end estimate of annual human consumption of plant material. The new
_f_ estimate is still 0.6% for 1983.

-------------

###Exercise 2

_The production of animal-derived foods, such as beef, eggs, fish, and milk,
requires the production of plants as fodder. To produce 1 J of energy in the
form of beef requires about 8 J of energy in the form of grains, while for
poultry about 3 J of energy from grains are required. These represent extremes.
The production of 1 J of other animal-derived foods requires very roughly 5 J of
plant matter. Estimate how much meat you eat per year and use this to work out
the following: If all of Earth's people ate a diet like yours, approximately
what would the fraction f be? (Hint: if you start with an estimate of the mass
of meat you eat, you will have to assume something about its water content. You
may assume fresh meat has about the same water content as fresh vegetation
roughly 70%)_

For this exercise, I make some general assumptions regarding my diet. I think on
a daily basis, I eat about 1 lbs. of meat type product a day. The proportions of
those meat products are about $\frac{1}{10}$ beef product, $\frac{1}{10}$
poultry product, and $\frac{4}{5}$ other meat product. The rest of my diet, I
assume I eat about 4 lbs of plant material a day. Using the values, I found the
amount of energy (kcal) derived from beef, poultry, and other (mostly eggs and
cheese) per every 100g from the USDA <https://ndb.nal.usda.gov/ndb/>.


    grams_per_lbs = 453.592
    water_content = 0.7 #approximate water content in meats
    lbs_meat = 1 *(1-water_content)
    my_meat_cons = lbs_meat * grams_per_lbs  #assume I eat 1 lbs (453.592g) of meat a day
    pr_beef = 0.1; pr_poul = 0.1; pr_other = 0.8 # proportions of beef, poultry, and other
    my_bf = my_meat_cons * pr_beef
    my_pf = my_meat_cons * pr_poul
    my_of = my_meat_cons * pr_other
    
    #USDA derived energy values for meat products in kcal/100g
    beef_ec = 270/100 #Used USDA #23573 Beef, ground, 80% lean meat / 20% fat, patty, cooked, broiled
    poultry_ec = 172/100 #Used USDA #05324, Chicken patty, frozen, cooked
    other_ec = (196/100 + 406/100)/2 #Used USDA #01128, Egg, whole, cooked, fried and #01009, Cheese, cheddar (averaged)


    #perform a small check first to make sure energy estimates seems reasonable
    cal_per_kcal = 1000
    beef_cal = my_bf * beef_ec * cal_per_kcal
    poul_cal = my_pf * poultry_ec * cal_per_kcal
    other_cal = my_of * other_ec * cal_per_kcal
    total_cal_meat = (beef_cal + poul_cal + other_cal) 
    print(total_cal_meat)

    387821.16000000003


So, from these numbers I get about 387.8 Calories from meat, which seems pretty
reasonable.


    total_J_meat = total_cal_meat * J_per_cal #convert to J
    lbs_plant = 4
    my_plant_cons = lbs_plant * grams_per_lbs
    my_dry_plant = my_plant_cons * (1 - water_content)
    my_J_plant = my_dry_plant * e_per_biomass 
    print(my_J_plant/J_per_cal)
    print(my_J_plant/J_per_cal+total_cal_meat)

    2083484.7846889955
    2471305.9446889954


I get about 2083.5 Calories from plant material and a total diet of 2471.3
Calories daily, so this sounds like a reasonable diet for the rest of the world,
especially following exercise 1? This diet is probably not representative, most
Americans eat 3000+ Calories per day and many people in 3rd world countries eat
much less.


    #convert the meat energy to equivalent in plant material
    bf_to_plant = 8 # beef to grain energy cost
    poul_to_plant = 3 # poultry energy cost
    oth_to_plant = 5 # other meat energy cost
    beef_plant_equi_J = beef_cal * J_per_cal * bf_to_plant
    poul_plant_equi_J = poul_cal * J_per_cal * poul_to_plant
    oth_plant_equi_J = other_cal * J_per_cal * oth_to_plant
    total_meat_plant_J = beef_plant_equi_J + poul_plant_equi_J + oth_plant_equi_J
    my_annual_J = (total_meat_plant_J + my_J_plant) * days_year
    human_annual_J = my_annual_J * population
    f_ex2 = human_annual_J/npp_J
    print("%.1f percent plant material consumed by humans annually (considering meat products)"%(f_ex2*100))

    1.0 percent plant material consumed by humans annually (considering meat products)


The _f_ of plant material needed by humans increased to 1.0%. So if every person
on the planet had a diet like mine (in 1983), the plant growth requirements
almost double.

-----------

###Exercise 3

_About what fraction of Earth's current npp would we need to consume if we
derived all the energy we now (1980) get from fossil fuel from biomass instead?
What does your answer tell you about the wisdom of replacing fossil fuels with
biomass? What ecological problems would you anticipate this might cause?_

For this exercise I looked up some type of published source of world wide fossil
fuel use. I found a report from World Watch Vital Signs (2003)
<https://www.worldwatch.org/brain/media/pdf/pubs/vs/2003_fossil_fuel.pdf>, which
is founded by Lester Brown a US environmental analyst. Although this is an
independent research institute, the estimates seem to be calculated from
reasonable sources (BP, DOE, IEA, IGU, and LBL). Using this as a source, from
1980, the world used 1814, 2972, and 1304 million tonnes of oil equivalent for
coal, oil, and natural gas respectively. We can convert this to energy (J)
assuming one toe (tonne of oil equivalent) to release 41.868 GJ.


    coal = 1814e6 #tonnes of oil equivalent
    oil = 2972e6 
    gas = 1304e6
    # fossil fuel taken from World Watch estimates
    f_J = 41.868e9 #joules in a tonne oil equivalent
    total_fos = (coal+oil+gas)*f_J
    print("%.1e joules of energy consumed in 1980." %(total_fos))

    2.5e+20 joules of energy consumed in 1980.


This seems reasonable since the world used 143,851 terawatt hours of energy in
2013 <http://en.wikipedia.org/wiki/World_energy_consumption>, which comes out to
5.2e20 joules. Seems like we doubled energy use since 1980.


    f_fos = total_fos/npp_J
    print("%.1f percent annual plant growth material equivalent used by fossil fuels"%(f_fos*100))

    8.5 percent annual plant growth material equivalent used by fossil fuels


If we were to replace all our fossil fuels with annual plant growth, we would
have to use 8.5% of the global NPP. So each year we would be using a substantial
portion of the global standing biomass just to power our needs. After 10-15
years we would run out of all our plant material, meaning we would lose a major
carbon sink, intensifying the carbon concentration in the atmosphere.
Ecologically, this much loss in plant material would mean loss of vital habitat
for wildlife/biodiversity, water filtration loss from wetlands, erosion control
from riparian plants, and a host of other ecosystem services.

----------

###Exercise 4

_If the human population continues to grow at about 2%/yr, in what year will
humans be eating Earth's current rate of npp?_

For this exercise, we have to write a compound growth and linear growth function
to determine the total Earth population for a given year projection.


    def compoundGrowth(p, yrs, growth_rate): 
        #given number of years, calculates total population at a 2% growth rate
        cP = p
        for i in range(yrs):
            pop_year = cP * growth_rate
            cP += pop_year
        return(cP)
    
    def linearGrowth(p, yrs, growth_rate):
        return(p+(p*growth_rate*yrs))
    
    #figure out the population for 2014 to see which growth rate model is reasonable
    growth_rate = 0.02
    currentyear = 2014
    dy = currentyear-1983
    pop_cg = compoundGrowth(population, dy, growth_rate)
    pop_lg = linearGrowth(population, dy, growth_rate)
    print("For compound growth the population estimate for year %i is %i individuals" %(currentyear, pop_cg))
    print("For linear growth the population estimate for year %i is %i individuals" %(currentyear, pop_lg))

    For compound growth the population estimate for year 2014 is 8683667434 individuals
    For linear growth the population estimate for year 2014 is 7614000000 individuals


Seems like the compound rate is a bit of an over estimate since we currently
have 7.19 billion individuals. We'll try both to see what's the difference.


    #Compound case
    annual_pop_J = pop_J
    yr = 0
    while (annual_pop_J<npp_J):
        annual_pop = compoundGrowth(population, yr, growth_rate)
        annual_pop_J = annual_pop * annual_J_human 
        yr+=1
    print("Total human population by year %i: %e" %(yr+1983,annual_pop))
    print("Annual human consumption by year %i: %.e J: " %(yr+1983,annual_pop_J))
    print("Annual NPP: %e"%(npp_J))
    print('With a compound growth rate of 2%% population per year, \nhumans would consume all plant material in %i years.'%(yr))


    Total human population by year 2243: 7.934916e+11
    Annual human consumption by year 2243: 3e+21 J: 
    Annual NPP: 2.999830e+21
    With a compound growth rate of 2% population per year, 
    humans would consume all plant material in 260 years.



    annual_pop_J = pop_J
    yr = 0
    while (annual_pop_J<npp_J):
        annual_pop = linearGrowth(population, yr, growth_rate)
        annual_pop_J = annual_pop * annual_J_human 
        yr+=1
    print("Total human population by year %i: %e" %(yr+1983,annual_pop))
    print("Annual human consumption by year %i: %.e J: " %(yr+1983,annual_pop_J))
    print("Annual NPP: %e"%(npp_J))
    print('With a linear growth rate of 2%% population per year, \nhumans would consume all plant material in %i years.'%(yr))


    Total human population by year 10301: 7.864980e+11
    Annual human consumption by year 10301: 3e+21 J: 
    Annual NPP: 2.999830e+21
    With a linear growth rate of 2% population per year, 
    humans would consume all plant material in 8318 years.


So, by year 2243 at the current human plant consumption rate and compound
population growth rate, we would consume all the global NPP with a total
population of ~793 billion individuals. With the linear growth model, it would
take us 8318 years (year 10,301) which is much more optimistic. I might go with
the more optimistic estimate, as population growth has been reported to slow
down globally and seems unlikely that humans will consume all of the NPP in the
near future.
