
t_cutoff_min <- 20
t_cutoff_max <- 100

## Run baseline matrix with no school contacts and increasing proportion of  removed adult contacts
all_dats <- NULL
actual_r0s <- NULL
thin_propns <- seq(0,1,by=0.05)
R0 <- 2
gamma <- 5

C <- contact_matrix(polymod,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                    missing.contact.age = "remove",
                    missing.participant.age = "remove")
C <- C$matrix
row.names(C) <- colnames(C)
beta_par <- get_beta(C,age_dat$population,gamma, R0)

all_C_thinned_base <- NULL

for(i in seq_along(thin_propns)){
  thin_propn <- thin_propns[i]
  ## Solve base model for comparison
  N_props <- age_distribution %>% pull(propn)
  N_age_classes <- length(N_props)
  N_immunity_classes <- 1
  
  ## Number of people in each age group and age class
  N <- matrix(N_props*N_tot,ncol=1,nrow=N_age_classes)
  
  ## Relative susceptibility of each age group
  alphas <- c(1)
  
  prop_immune <- rep(0, N_age_classes)
  beta_scales <- rep(1, N_age_classes)
  alphas <- c(1)
  
  N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
  N_props_long <- N_props
  N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)
  
  beta_scales <- rep(1,9)
  
  contacts <- polymod$contacts
  contacts_no_schools <- contacts #%>% filter(cnt_school == 0)
  polymod1 <- polymod
  polymod1$contacts <- contacts_no_schools
  
  polymod_thinned <- thin_polymod(polymod1, t_cutoff_min, t_cutoff_max, thin_propn, TRUE)
  polymod_c_thin <- contact_matrix(polymod_thinned,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                                   missing.contact.age = "remove",
                                   missing.participant.age = "remove")
  C_thin <- polymod_c_thin$matrix
  row.names(C_thin) <- colnames(C_thin)
  
  tmp_C <- reshape2::melt(C_thin)
  colnames(tmp_C) <- c("age_group","contact_age_group","value")
  tmp_C$thin_propn <- thin_propn
  all_C_thinned_base <- rbind(all_C_thinned_base, tmp_C)
  
  C_use <- setup_C(C_thin, N, beta_scales)
  
  C_use_explicit <- setup_C_explicit(C_thin, N, beta_scales)
  actual_r0 <- 2
  #actual_r0 <- max(eigen(C_use_explicit)$values)*beta_par*gamma
  
  y <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(0,2000,by=1),alphas=alphas, 
                    age_seeds=1:9,immunity_seed=1,return_peak=TRUE)
  t_peak <- y[[2]]
  y <- y[[1]]
  
  dat <- data.frame(age_group=china_cases$age_group, ar=y[,1], thin_propn=thin_propn,r0=actual_r0,
                    peak_time=t_peak)
  all_dats <- bind_rows(all_dats, dat)
}

beta_scales1 <- sigmoid_func(seq(0,80,by=10),0.1, 60, 0.01)

## Run baseline matrix with no school contacts and increasing proportion of  removed adult contacts
all_dats2 <- NULL
actual_r0s <- NULL
R0 <- 2

C <- contact_matrix(polymod,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                    missing.contact.age = "remove",
                    missing.participant.age = "remove")
C <- C$matrix
row.names(C) <- colnames(C)
beta_par <- get_beta_vector(C,age_dat$population,gamma, R0,beta_scales1)

all_C_thinned <- NULL

for(i in seq_along(thin_propns)){
  thin_propn <- thin_propns[i]
  ## Solve base model for comparison
  N_props <- age_distribution %>% pull(propn)
  N_age_classes <- length(N_props)
  N_immunity_classes <- 1
  
  ## Number of people in each age group and age class
  N <- matrix(N_props*N_tot,ncol=1,nrow=N_age_classes)
  
  ## Relative susceptibility of each age group
  prop_immune <- rep(0, N_age_classes)
  alphas <- c(1)
  
  N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
  N_props_long <- N_props
  N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)
  
  
  contacts <- polymod$contacts
  contacts_no_schools <- contacts# %>% filter(cnt_school == 0)
  polymod1 <- polymod
  polymod1$contacts <- contacts_no_schools
  
  polymod_thinned <- thin_polymod(polymod1, t_cutoff_min, t_cutoff_max, thin_propn, TRUE)
  polymod_c_thin <- contact_matrix(polymod_thinned,survey.pop=age_dat,age.limits = seq(0,80,by=10),
                                   symmetric=TRUE,
                                   missing.contact.age = "remove",
                                   missing.participant.age = "remove")
  C_thin <- polymod_c_thin$matrix
  row.names(C_thin) <- colnames(C_thin)
  tmp_C <- reshape2::melt(C_thin)
  colnames(tmp_C) <- c("age_group","contact_age_group","value")
  tmp_C$thin_propn <- thin_propn
  all_C_thinned <- rbind(all_C_thinned, tmp_C)
  
  C_use <- setup_C(C_thin, N, beta_scales1)
  
  C_use_explicit <- setup_C_explicit(C_thin, N, beta_scales1)
  #actual_r0 <- max(eigen(C_use_explicit)$values)*beta_par*gamma
  actual_r0 <- 2
  y <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(0,2000,by=1),alphas=alphas, 
                    age_seeds=1:9,immunity_seed=1,return_peak=TRUE)
  t_peak <- y[[2]]
  y <- y[[1]]
  
  dat <- data.frame(age_group=china_cases$age_group, ar=y[,1], thin_propn=thin_propn,r0=actual_r0,
                    peak_time=t_peak)
  all_dats2 <- bind_rows(all_dats2, dat)
}
C <- contact_matrix(polymod,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                    missing.contact.age = "remove",
                    missing.participant.age = "remove")
C <- C$matrix
row.names(C) <- colnames(C)
beta_scales_none <- rep(1,9)
C_use <- setup_C(C, N, beta_scales_none)
beta_par <- get_beta(C, N, gamma, 2)
y <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(0,2000,by=1),alphas=alphas, 
                  age_seeds=1:9,immunity_seed=1, return_peak=TRUE)
t_peak <- y[[2]]
y <- y[[1]]
dat <- data.frame(age_group=china_cases$age_group, ar=y[,1], thin_propn=thin_propn,r0=actual_r0)


r0_labels <- unique(all_dats[all_dats$age_group == "80+" & all_dats$r0 > 1.1,c("ar","r0")])
colnames(all_dats)[3] <- "Proportion of adult-adult contacts removed"
p3 <- ggplot(all_dats) + 
  geom_line(aes(x=as.numeric(age_group),y=ar,group=`Proportion of adult-adult contacts removed`,
                col=`Proportion of adult-adult contacts removed`),size=0.4) + 
  geom_line(data=dat,aes(x=as.numeric(age_group),y=ar),size=0.5,col="grey40") +
  geom_point(aes(x=as.numeric(age_group),y=ar,
                 group=`Proportion of adult-adult contacts removed`,
                 col=`Proportion of adult-adult contacts removed`),size=0.5) + 
  geom_point(data=dat,aes(x=as.numeric(age_group),y=ar),size=0.4,col="grey40") +
  #geom_text_repel(data=r0_labels, aes(label=signif(r0,3), y=ar, x= 10)) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(breaks=seq(1,9,by=1), labels=rownames(C)) +
  scale_color_gradientn(colors=c("#0072B2","#56B4E9","#009E73","#E69F00","#D55E00")) +
  theme_pubr() +
  guides(colour = guide_colorbar(title.position = "top", title.hjust = 0.5,barwidth = 8,barheight=1)) +
  theme(legend.position=c(0.7,0.95),
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.title=element_text(size=6),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8))  +
  ylab("Final cumulative incidence") +
  xlab("Age group")
C_use <- setup_C(C, N, beta_scales1)

C <- contact_matrix(polymod,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                    missing.contact.age = "remove",
                    missing.participant.age = "remove")
C <- C$matrix
row.names(C) <- colnames(C)
beta_par <- get_beta_vector(C,age_dat$population,gamma, R0,beta_scales1)
y <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(0,2000,by=1),alphas=alphas, age_seeds=1:9,immunity_seed=1)
dat1 <- data.frame(age_group=china_cases$age_group, ar=y[,1], thin_propn=thin_propn,r0=actual_r0)

r0_labels2 <- unique(all_dats2[all_dats2$age_group == "80+" & all_dats2$r0 > 1.1,c("ar","r0")])
p4 <- ggplot(all_dats2) + 
  geom_line(aes(x=as.numeric(age_group),y=ar,group=thin_propn,col=thin_propn),size=0.4) + 
  geom_line(data=dat1,aes(x=as.numeric(age_group),y=ar),size=0.4,col="grey40") +
  geom_point(aes(x=as.numeric(age_group),y=ar,group=thin_propn,col=thin_propn),size=0.5) + 
  geom_point(data=dat1,aes(x=as.numeric(age_group),y=ar),size=0.5,col="grey40") +
  #geom_text_repel(data=r0_labels2, aes(label=signif(r0,3), y=ar, x= 10)) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(breaks=seq(1,9,by=1), labels=rownames(C)) +
  scale_color_gradientn(colors=c("#0072B2","#56B4E9","#009E73","#E69F00","#D55E00")) +
  #scale_color_gradientn(low="#0072B2",high="red",mid="green",midpoint=0.3)+
  theme_pubr() +
  guides(colour = guide_colorbar(title.position = "top", title.hjust = 0.5,barwidth = 10)) +
  theme(legend.position="none",
        legend.text=element_text(size=6),
        legend.title=element_text(size=6),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8)
  ) +
  ylab("") +
  xlab("Age group")

#png("fig2.png",width=8,height=3,res=300,units="in")
png("tmp_workers.png",width=8, height=5, res=300,units="in")
(p3+p1)/(p4 + p2)
dev.off()

pdf("tmp_workers.pdf",width=8, height=5)
(p3+p1)/(p4 + p2)
dev.off()
