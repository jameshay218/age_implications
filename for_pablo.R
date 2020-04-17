## Functions we need
setwd("~/Documents/Github/age_implications/")

source("auxiliary_funcs.R")
source("sir_functions.R")

## Libraries we need
library(socialmixr)
library(pracma)
library(deSolve)
library(ggplot2)
library(tidyverse)

## Demography data
## Age distribution of China
age_distribution <- read_csv("data/population_demography_10year.csv")

## Population of Wuhan
N_tot <- 11080000

## POLYMOD contact matrix
## DEAR PABLO - THIS IS WHERE WE READ IN THE POLYMOD SOCIAL MIXING DATA. HAVE A LOOK AT IT
data("polymod")
View(polymod$contacts)

## Proportion of poulation in each age class
N_props1 <- c(0.118574496482027, 0.115752336341461, 0.128634781768929, 0.15898446700289, 
              0.150148368364849, 0.154367713154031, 0.105371761585134, 0.0496726586040798, 
              0.0184934166965997)
age_dat <- data.frame(lower.age.limit=seq(0,80,by=10),population=N_props1)
polymod_c <- contact_matrix(polymod,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                            missing.contact.age = "sample",
                            missing.participant.age = "remove")
C <- polymod_c$matrix
row.names(C) <- colnames(C)

## Create polymod with no school contacts
## ie. filter out all school contacts
contacts <- polymod$contacts
contacts_no_schools <- contacts %>% filter(cnt_school == 0)# & cnt_leisure == 0)
polymod1 <- polymod
polymod1$contacts <- contacts_no_schools
polymod_c_no_schools <- contact_matrix(polymod1,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                                       missing.contact.age = "sample",
                                       missing.participant.age = "remove")
C_no_schools <- polymod_c_no_schools$matrix
row.names(C_no_schools) <- colnames(C_no_schools)

## Proportion of population in each age group
age_distribution <- age_distribution %>% mutate(propn = number/sum(number))
age_distribution$age_group <- factor(age_distribution$age_group, levels=unique(age_distribution$age_group))
age_distribution <- age_distribution %>% mutate(wuhan=N_tot*propn)

## Solve base model for comparison
N_props <- age_distribution %>% pull(propn)
N_age_classes <- length(N_props1)

## Ignore this next parameter - not used
N_immunity_classes <- 1

## Number of people in each age group and age class
N <- matrix(N_props*N_tot,ncol=1,nrow=N_age_classes)

## Relative susceptibility of each age group - IGNORE, NOT USED
alphas <- c(1)

R0 <- 2
gamma <- 5

## IGNORE - NOT USED
prop_immune <- rep(0, N_age_classes)
alphas <- c(1)

## This is used to scale relative transmissibility by age
beta_scales <- rep(1, N_age_classes)

N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
N_props_long <- N_props
N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)

beta_scales <- rep(1,9)

## Normalise the contact matrix
C_use <- setup_C(C, N, beta_scales)
## Find the beta parameter needed to get desired R0, given this contact matrix and gamma
beta_par <- get_beta(C,age_dat$population,gamma, R0)

## Solve model
y_base <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(1,365,by=1),
                       alphas=alphas, age_seed=4,immunity_seed=1)

## Do the same thing but after school contacts removed
C_use <- setup_C(C_no_schools, N, beta_scales)
y_noschools <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(1,365,by=1),
                            alphas=alphas, age_seed=4,immunity_seed=1)


plot(y_base,type='l',col="red",ylim=c(0,1),ylab="Final cumulative incidence", xlab="Age group")
## Results with all school contacts removed
lines(y_noschools,col="green")
