### parameters
# implementation parameters and age ranges
ages<-30			# maximum age
years<-250			# simulation years
skip<-200			# transient years to be skipped
N.ini<-1e8			# initial numbers at age 1

# biological parameters
r0<-21.77072		# maximum age 1 recruits per kg SSB
Bhalf<-365426284	# half-saturation constant in recruitment ***tuned***
Linf<-196			# VB: asymptotic length
a0<- -0.937			# VB: age at length zero
k<-0.0555			# VB: steepness
theta<-8.10e-6		# weight-length: condition
g<-3.01				# weight-length: allometric exponent
# age-specific biological parameters from 2004-2013 (AFWG2014)
M<-c(1.48936,0.61272,0.36018,0.26108,0.22300,0.20447,0.20000) # natural mortality 
mat<-c(0,0,0.0002,0.0029,0.0565,0.3240,0.6361,0.8451,0.9473,0.9838,0.9944,0.9964,1) # maturity-at-age
sf<-0.4858775		# # recruitment noise at age 3 from the last 30 years (AFWG2014)

# fishing parameters
LF50<-61.4806 		# length at 50% selectivity estimated from length-at-age and selection pattern (F in 2004-2013, AFWG2014)
SF<-0.1222			# parameter related to steepness, as above
b<-0.75				# concentration profile parameter
q<-2.827164e-06		# scaling parameter relating to catchability and density ***tuned***

# employment parameters
d.sea<-0.054		# FTE/vessel day
d.max<-30000		# maximum available offshore FTE
d.shr<-0.000004		# FTE/kg

# economic parameters
price.sea<-13.13	# landing price NOK/kg
price.shr<-17.0		# selling price NOK/kg

s.sea<-1078000			# employment cost sea NOK/FTE
s.shr<-348000			# employment cost shore NOK/FTE
fixedC.sea<-351123000	# fixed costs sea NOK (=average per unit * #units)
fixedC.shr<-103246800	# fixed costs shore NOK
varC.sea<-65000 		# variable costs NOK/vessel day
scale.catch<-0.53 		# percentage of total codfish catch that is cod
