## Lab 10: age-structured model

########################################################################################
## PART 1: Two AGE GROUP MODEL
########################################################################################

library(deSolve)

# sample code for the 2-age group model
SIR2ageGrs <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dSC = nu - SC * (betaCC * IC + betaCA * IA) - muC * SC - lC * SC
    dIC = SC * (betaCC * IC + betaCA * IA) - gamma * IC - muC * IC - lC * IC
    dSA = lC * SC - SA * (betaAC * IC + betaAA * IA) - muA * SA
    dIA = lC * IC + SA * (betaAC * IC + betaAA * IA) - gamma * IA - muA * IA
    
    list(c(dSC, dIC, dSA, dIA))
  })
}



## parameters/initial conditions 
## note all time units are in year here
betaCC = 1100; # transmission rate among children per year
betaCA = betaAC = 110;  # transmission rate between children and adults per year
betaAA = 220; # transmission rate among adults per year
gamma = 1 / (14/365); # recovery rate: 1/ infectious period, infectious period = 2 weeks  

lC = 0.066667; # rate of leaving the children group: 1/15 per year

muC = 0; # mortality rate of children
muA = 0.0166667; # mortality rate of adults:  1/(75-15);  mean life expectancy=75yr; subtract the 15 yrs in youth with 0 death rate

# fraction of population in the two age groups
nC = 0.2;  # % children = 15/75; 
nA = 0.8; # % adults = 1-nC;
nu = 1/75; # birth rate, per year

# initial susceptibility 
SC0 = .1; # initial susceptibility in children 10%
SA0 =.1; # initial susceptibility in adults 10%
# initial infectious:
IC0 = .0001; # initial infection rate in children 
IA0 = .0001; # initial infection rate in adults
# initial recovered:
RC0 = nC - IC0 - SC0; 
RA0 = nA - IA0 - SA0; 

######  CODE ON YOUR OWN  ###### 
## to run the model, you will have to put all pieces (parameters, initial conditions, etc.) together. 
## [Q1] Run the model using parameters/initial conditions listed in the next slide for 100 years.  
## Plot the fraction of susceptibles (Si/ni) for the two groups. what do you find, which group has a higher susceptibility? Why? (1pt)


## SAMPLE CODE TO NORMALIZE THE PROPORTIONS RELATIVE TO THE SPECIFIC GROUP (VS THE ENTIRE POPULATION)
if(F){
  ## proportion susceptible in the two groups, i.e. Si/ni
  fSC=sim[,'SC']/nC; # note: it is normalized by population size in each group
  fSA=sim[,'SA']/nA;
  
  ## proportion infectious in the two groups, i.e. Ii/ni
  fIC=sim[,'IC']/nC;
  fIA=sim[,'IA']/nA;
  
  
  ## plot # people infected
  totI=rowSums(sim[,c('IC','IA')]);  # total infectious
  totS=rowSums(sim[,c('SC','SA')]);  # total susecptible
  
}





# [Q2] What is R0 for the entire population? What would the average age of infection be, given this R0 value and a life span of 75 yr? (0.5pt)
# [Hint: use the same method we used for risk structured model, ignore death rate and aging rate]
# [Q3] What is the force of infection for each of the two groups, when the disease / your simulation reaches equilibrium? (0.5pt)[Hint: What is the force of infection?]
# [Q4] What is the force of infection for each of the two groups, given the parameters and initial state variables at the beginning of the simulation? (0.5pt)[Hint: What is the force of infection?]


########################################################################################
## PART 2: VACCINATION
########################################################################################
# [Q5] Code a model with 2 age groups and vaccination at birth using the equations in the slide 12
# (Submit the finished code in your lab report.)  





# [Q6] Then test your model with p=0.5 and the same initial conditions/parameters as in Part 1. 
# Plot the fraction of recovered/immune (Ri/ni) for the two age groups.   (0.5 pt)






####################
## TEST DIFFERENT VACCINATION RATES
## use the R Shiny App: ShinyApp_Vac.R
