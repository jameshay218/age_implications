## Functions we need
source("auxiliary_funcs.R")
source("sir_functions.R")

## Libraries we need
library(socialmixr)
library(pracma)
library(deSolve)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(doParallel)

rerun <- TRUE
## Set up parallelisation
n_clusters <- 8
#registerDoParallel(cores=n_clusters)
cl <- makeCluster(n_clusters)
registerDoParallel(cl)

## Demography data
## Distribution of cases
china_cases <- read_csv("data/case_data_2020.02.11.csv")

## Age distribution of China
age_distribution <- read_csv("data/population_demography_10year.csv")

## Population of Wuhan
N_tot <- 11080000

## POLYMOD contact matrix
data("polymod")
N_props1 <- c(0.118574496482027, 0.115752336341461, 0.128634781768929, 0.15898446700289, 
              0.150148368364849, 0.154367713154031, 0.105371761585134, 0.0496726586040798, 
              0.0184934166965997)
age_dat <- data.frame(lower.age.limit=seq(0,80,by=10),population=N_props1)
polymod_c <- contact_matrix(polymod,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                            missing.contact.age = "sample",
                            missing.participant.age = "remove")
C <- polymod_c$matrix
row.names(C) <- colnames(C)

## Create polymod with No school contacts
contacts <- polymod$contacts
contacts_no_schools <- contacts %>% filter(cnt_school == 0)# & cnt_leisure == 0)
polymod1 <- polymod
polymod1$contacts <- contacts_no_schools
polymod_c_no_schools <- contact_matrix(polymod1,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                            missing.contact.age = "sample",
                            missing.participant.age = "remove")
C_no_schools <- polymod_c_no_schools$matrix
row.names(C_no_schools) <- colnames(C_no_schools)

######################
## FIGURE 1
######################
C_melted <- reshape2::melt(C)
C_melted$Scenario <- "POLYMOD"

C_melted_noschool <- reshape2::melt(C_no_schools)
C_melted_noschool$Scenario <- "No school contacts"

C_all <- bind_rows(C_melted, C_melted_noschool)
C_all$Scenario <- factor(C_all$Scenario, levels=c("POLYMOD","No school contacts"))

p_all_C <- ggplot(C_all) + 
  geom_tile(aes(x=Var1,y=contact.age.group,fill=value)) + 
  scale_fill_gradient2(low="#0072B2",high="#D55E00",mid="#F0E442", midpoint=6) +
  theme_pubr() +
  ylab("Age group of contact") +
  xlab("Age group") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(legend.position="none",
        axis.text=element_text(size=6),
        axis.title=element_text(size=10),
        plot.margin=unit(c(0.2,0.25,0.2,0.2),"cm")) +
  facet_wrap(~Scenario)

png("figS1.png",width=8,height=3,res=300,units="in")
p_all_C
dev.off()

p_C <- ggplot(C_melted) + 
  geom_tile(aes(x=Var1,y=contact.age.group,fill=value)) + 
  scale_fill_gradient2(low="#0072B2",high="#D55E00",mid="#F0E442", midpoint=6) +
  theme_pubr() +
  ylab("Age group of contact") +
  xlab("Age group") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(legend.position="none",
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        plot.margin=unit(c(0.2,0.25,0.2,0.2),"cm"))


## Proportion of population in each age group
age_distribution <- age_distribution %>% mutate(propn = number/sum(number))
age_distribution$age_group <- factor(age_distribution$age_group, levels=unique(age_distribution$age_group))
age_distribution <- age_distribution %>% mutate(wuhan=N_tot*propn)
age_dist_p <- age_distribution %>% ggplot() + 
  geom_bar(aes(x=age_group,y=propn),stat="identity",fill="grey40") +
  theme_pubr() + 
  xlab("Age group") +
  ylab("Proportion of population")+ 
  coord_flip() +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=10),
        axis.text.x=element_text(size=7),
        plot.margin=unit(c(0.2,0.4,0.2,0.2),"cm")) +
  labs(tag="A")

## Observed proportion of each age group infected
china_cases <- china_cases %>% left_join(age_distribution)
china_cases <- china_cases %>% mutate(incidence=confirmed_cases*0.75/(propn*58500000))

case_dist_p <- china_cases %>% ggplot() + 
  geom_bar(aes(x=age_group,y=confirmed_cases),stat="identity", fill="grey40") +
  theme_pubr() + 
  xlab("Age group") +
  ylab("Confirmed cases")+ 
  coord_flip() +
  scale_y_continuous(expand=c(0,0),limits=c(0,12000)) +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=10),
        axis.text.x=element_text(size=7),
        plot.margin=unit(c(0.2,0.4,0.2,0.2),"cm"))+
  labs(tag="B")



## Solve base model for comparison
N_props <- age_distribution %>% pull(propn)
N_age_classes <- length(N_props1)
N_immunity_classes <- 1

## Number of people in each age group and age class
N <- matrix(N_props*N_tot,ncol=1,nrow=N_age_classes)

## Relative susceptibility of each age group
alphas <- c(1)

R0 <- 2
gamma <- 5

prop_immune <- rep(0, N_age_classes)
beta_scales <- rep(1, N_age_classes)
alphas <- c(1)

N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
N_props_long <- N_props
N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)

beta_scales <- rep(1,9)

C_use <- setup_C(C, N, beta_scales)
beta_par <- get_beta(C,age_dat$population,gamma, R0)
y_base <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(1,365,by=1),
                       alphas=alphas, age_seed=4,immunity_seed=1)

C_use <- setup_C(C_no_schools, N, beta_scales)
y_noschools <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(1,365,by=1),
                           alphas=alphas, age_seed=4,immunity_seed=1)

######################
## FIGURE 2
######################
####################################
## GRID SEARCH OF PARAMETER SPACE
####################################
to_test <- expand.grid(transmissibility=seq(0.03,0.2,by=0.001),
                       susceptibility_scale=1,
                       relative_susceptibility=c(0),
                       minimum_transmissibility_scale=c(0,0.00001,0.0001,0.001,0.01),
                       transmissibility_scale=seq(0.1,0.3,by=0.1),
                       transmissibility_midpoint=seq(40,70,by=5))

## Split up to send to different cores
to_test$cluster <- n_clusters
cluster_indices <- rep(1:n_clusters, each=nrow(to_test)/n_clusters)
to_test$cluster[1:length(cluster_indices)] <- cluster_indices

prop_immune <- rep(0, N_age_classes)
ages <- seq(0,80,by=10)
N_immunity_classes <- 2

## Solve model overall whole grid search
if(rerun){
  res <- foreach(i=1:n_clusters,.packages=c("pracma", "deSolve", "foreach"), .combine=rbind) %dopar% {
  #for(i in 1:n_clusters){
  ## Get subset of parameter combinations to try on this core
  subset_pars <- to_test[to_test$cluster == i,]
  
  ## Iterate through each row
  tmp_res <- foreach(j=1:nrow(subset_pars),.combine=rbind) %do% {
    #print(j)
    beta_par <- subset_pars$transmissibility[j]
    alpha_immune <- subset_pars$relative_susceptibility[j]
    A <- subset_pars$minimum_transmissibility_scale[j]
    k <- subset_pars$transmissibility_scale[j]
    x0 <- subset_pars$transmissibility_midpoint[j]
    
    alphas <- c(1, alpha_immune)
    alphas <- 1
    N_immunity_classes <- 1

    N_props_long <- c(N_props, N_props*0)
    N <- matrix(N_props_long*N_tot, ncol=N_immunity_classes, nrow=N_age_classes)

    beta_scales1 <- sigmoid_func(seq(0,80,by=10),k,x0,A)
    
    
    C_use <- setup_C(C, N, beta_scales1)
    y <- epi_ode_size(C_use, beta_par, gamma, N, alphas, seq(0,365,by=1),age_seed=9,immunity_seed = 1)
    
    y[!is.finite(y)] <- 0
    
    ## Get number infected in each age and immunity class, then remnormalise to get average incidence
    ## per age class
    numbers <- rowSums(y*N)
    popn_size <- rowSums(N)
    overall_i <- sum(numbers)
    overall_n <- sum(popn_size)
    
    numbers <- rowSums(y*N)/rowSums(N)
    tmp <- c(numbers, overall_i/overall_n)
    names(tmp) <- c(colnames(C),"overall_ar")
    tmp
  }
  final <- cbind(subset_pars, tmp_res)
}
  write_csv(res, "res.csv")
} else {
  res <- read_csv("res.csv")
}
res2 <- res %>% pivot_longer(cols=c("[0,10)", "[10,20)", "[20,30)", "[30,40)", "[40,50)", "[50,60)", 
                                    "[60,70)", "[70,80)", "80+"))
ggplot(res2[res2$transmissibility_scale==0.1 & res2$minimum_transmissibility_scale < 0.5,]) + 
  geom_line(aes(x=name,y=value, col=transmissibility, group=transmissibility)) +# facet_wrap(~relative_susceptibility) +
  scale_color_gradient2(low="blue",mid="orange",high="red",midpoint=0.09) +
  facet_grid(minimum_transmissibility_scale~transmissibility_midpoint)

sigmoid_res <- get_sigmoid_func_outputs(ages, 
                                        unique(to_test[,c("minimum_transmissibility_scale",
                                                          "transmissibility_midpoint","transmissibility_scale")]))
ggplot(sigmoid_res[sigmoid_res$minimum_transmissibility_scale < 0.5 &
                     sigmoid_res$transmissibility_scale == 0.1,]) + 
  geom_line(aes(x=age_group,y=beta,col=transmissibility_scale, group=transmissibility_scale)) + 
  facet_grid(minimum_transmissibility_scale~transmissibility_midpoint)


final_ar <- res2[res2$transmissibility_scale==0.1 & 
                   res2$minimum_transmissibility_scale == 0.01 &
                   res2$transmissibility_midpoint == 60 &
                   res2$transmissibility == 0.16,] %>% select(name, value)


C <- contact_matrix(polymod,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                    missing.contact.age = "remove",
                    missing.participant.age = "remove")
C <- C$matrix
row.names(C) <- colnames(C)

beta_par <- get_beta_vector(C,age_dat$population,gamma, R0,beta_scales)

ar_base <- data.frame(age_group=china_cases$age_group, ar=y_base[,1], Scenario="POLYMOD")
ar_noschool <- data.frame(age_group=china_cases$age_group, ar=y_noschools[,1], Scenario="No school contacts")
ars <- bind_rows(ar_base, ar_noschool)

beta_scales <- sigmoid_func(seq(0,80,by=10),0.1, 60, 0.01)
C_use1 <- setup_C(C, N, beta_scales)
y_model <- epi_ode_size(C_use1, 0.16, gamma, N,ts=seq(1,365,by=1),
                        alphas=alphas, age_seed=4,immunity_seed=1)
ar_model <- data.frame(age_group=china_cases$age_group, ar=y_model[,1],Scenario="Age-dependent transmissibility")
ars <- bind_rows(ars, ar_model)

y_scale <- ar_model %>% filter(age_group == "80+") %>% summarise(max(ar)) %>% pull(`max(ar)`)
y_scale <- as.numeric(y_scale/china_cases[china_cases$age_group == "80+","incidence"])


## Fig 2 inset, baseline AR
ars$Scenario <- factor(ars$Scenario, 
                       levels=c("POLYMOD","No school contacts","Age-dependent transmissibility"))
inset_p <- ggplot(ars) +
  geom_bar(data=china_cases, aes(x=age_group,y=incidence*y_scale),stat="identity",fill="grey40") +
  geom_line(data=ars,aes(x=as.numeric(age_group),y=ar, col=Scenario),size=0.5) +
  geom_point(aes(x=age_group,y=ar, col=Scenario),size=1) +
  xlab("Age group") +
  scale_y_continuous(expand=c(0,0), limits=c(0,1), "Predicted cumulative incidence (line)", breaks=seq(0,1,by=0.2),
                     sec.axis=sec_axis(~.*(1/y_scale), name="Observed incidence (bars)*", breaks=seq(0,0.002,by=0.0002))) +
  scale_color_manual(values=c("#D55E00","#0072B2","#009E73")) +
  theme_pubr() +
  theme(axis.text=element_text(size=8),
        axis.text.x=element_text(size=8),
        axis.title=element_text(size=10),
        legend.title = element_blank(),
        legend.position=c(0.5,1.1),
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.background=element_blank(),
        legend.text=element_text(size=6),
        plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))

inset_p
y <- sigmoid_func(0:80,0.1, 60, 0.01)
age_dist_trans <- data.frame(beta=y,age_group=0:80)
beta_par <- get_beta(C,age_dat$population,gamma, R0)
age_dist_trans$age_group <- as.numeric(age_dist_trans$age_group)
transmissibility_dist <- ggplot(age_dist_trans) +  
  geom_hline(yintercept=beta_par, col="#0072B2",linetype="dashed", size=0.5) +
  geom_hline(yintercept=max(age_dist_trans$beta)*0.16,col="#D55E00", 
             linetype="dashed", size=0.5) +
  geom_line(aes(x=as.numeric(age_group),y=beta*0.16),size=0.5,col="grey40") +
  xlab("Age group") +  
  ylab("Transmissibility") +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.15),breaks=seq(0,0.15,by=0.01)) +
  scale_x_continuous(expand=c(0,0)) +
  theme_pubr() +
  theme(axis.text=element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.title=element_text(size=10),
        legend.title = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))
transmissibility_dist

## Blue line is R0 of 2 in polymod mixing world
## Orange line is transmissibility in 

rhs <- inset_p/transmissibility_dist + plot_layout(heights=c(1.2,1))
#png("figS1.png",width=10,height=6,res=300,units="in")
(((age_dist_p | case_dist_p)/p_C) | rhs) + 
  plot_layout(widths=c(1,1.2))
#dev.off()


#pdf("figS1.pdf",width=10,height=6)
(((age_dist_p | case_dist_p)/p_C) | rhs) + 
  plot_layout(widths=c(1,1.2))
#dev.off()

