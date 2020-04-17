library(ggrepel)

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
  contacts_no_schools <- contacts %>% filter(cnt_school == 0)
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
  contacts_no_schools <- contacts %>% filter(cnt_school == 0)
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
                  age_seeds=1:9,immunity_seed=1,return_peak=TRUE)
t_peak <- y[[2]]
y <- y[[1]]

dat <- data.frame(age_group=china_cases$age_group, ar=y[,1], thin_propn=thin_propn,r0=actual_r0,
                  peak_time=t_peak)


r0_labels <- unique(all_dats[all_dats$age_group == "80+" & all_dats$r0 > 1.1,c("ar","r0")])
colnames(all_dats)[3] <- "Proportion of adult-adult contacts removed"
p1 <- ggplot(all_dats) + 
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
  theme(legend.position="none",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.title=element_text(size=6),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8)) +
  ylab("Final cumulative incidence") +
  xlab("Age group")
p1
C_use <- setup_C(C, N, beta_scales1)

C <- contact_matrix(polymod,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                    missing.contact.age = "remove",
                    missing.participant.age = "remove")
C <- C$matrix
row.names(C) <- colnames(C)
beta_par <- get_beta_vector(C,age_dat$population,gamma, R0,beta_scales1)
y <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(0,2000,by=1),alphas=alphas, 
                  age_seeds=1:9,immunity_seed=1,return_peak=TRUE)
t_peak <- y[[2]]
y <- y[[1]]

dat1 <- data.frame(age_group=china_cases$age_group, ar=y[,1], thin_propn=thin_propn,r0=actual_r0,
                   peak_time=t_peak)

r0_labels2 <- unique(all_dats2[all_dats2$age_group == "80+" & all_dats2$r0 > 1.1,c("ar","r0")])
p2 <- ggplot(all_dats2) + 
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
(p1 + p2)
#dev.off()
polymod_c_normal <- contact_matrix(polymod,survey.pop=age_dat,age.limits = seq(0,80,by=10),
                                   symmetric=TRUE,
                                   missing.contact.age = "remove",
                                   missing.participant.age = "remove")
C_thin <- polymod_c_normal$matrix
row.names(C_thin) <- colnames(C_thin)
tmp_C <- reshape2::melt(C_thin)
colnames(tmp_C) <- c("age_group","contact_age_group","value")
tmp_C$thin_propn <- "Baseline"
all_C_thinned <- rbind(all_C_thinned, tmp_C)
all_C_thinned$thin_propn <- factor(all_C_thinned$thin_propn, levels=c("Baseline",
                                                                      as.character(seq(0,1,by=0.05))))
p2_supp <- ggplot(all_C_thinned[all_C_thinned$thin_propn %in% c(as.character(seq(0,1,by=0.1)),"Baseline"),]) + 
  geom_tile(aes(x=age_group,y=contact_age_group,fill=value)) + 
  scale_fill_gradient2(low="#0072B2",high="#D55E00",mid="#F0E442",midpoint=6) +
  theme_pubr() +
  ylab("Age group of contact") +
  xlab("Age group") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(legend.position="none",
        axis.text=element_text(size=8),
        axis.text.x=element_text(size=8,angle=90,hjust=0.5,vjust=0.5),
        axis.title=element_text(size=10),
        plot.margin=unit(c(0.2,0.25,0.2,0.2),"cm")) +
  facet_wrap(~thin_propn)

#pdf("fig2_supp.pdf",height=5,width=7)
p2_supp
#dev.off()


all_C_thinned_fig1 <- all_C_thinned %>% filter(thin_propn %in% c("0","0.5","Baseline"))
new_labels <- c("0"="No school contacts","0.5"="50% reduction in adults","Baseline"="POLYMOD")
all_C_thinned_fig1$thin_propn <- new_labels[as.character(all_C_thinned_fig1$thin_propn)]
all_C_thinned_fig1$thin_propn <- factor(all_C_thinned_fig1$thin_propn, levels=c("POLYMOD","No school contacts","50% reduction in adults"))
p_new_fig1_C <- ggplot(all_C_thinned_fig1) + 
  geom_tile(aes(x=age_group,y=contact_age_group,fill=value)) + 
  scale_fill_gradient2(low="#0072B2",high="#D55E00",mid="#F0E442",midpoint=6) +
  theme_pubr() +
  ylab("Age group of contact") +
  xlab("Age group") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(legend.position="none",
        axis.text=element_text(size=7),
        axis.text.x=element_text(size=7,angle=90,hjust=0.5,vjust=0.5),
        axis.title=element_text(size=10),
        strip.text=element_text(size=8),
        plot.margin=unit(c(0.2,0.1,0.2,-0.1),"cm")) +
  facet_wrap(~thin_propn,ncol=1) +
  labs(tag="C")

transmissibility_dist <- transmissibility_dist + labs(tag="D")
lhs <- ((age_dist_p / case_dist_p / transmissibility_dist) | p_new_fig1_C) + plot_layout(widths=c(1,1.2))
lhs  
inset_p1 <- inset_p + theme(legend.position=c(0.55,0.95),legend.justification = 0.5,
                            plot.margin = margin(0.2,1,0.2))+ labs(tag="E")

inset_p1 <- inset_p1 + plot_spacer() + plot_layout(widths=c(20,0.1))
inset_p1
overall_p <- lhs / inset_p1 + plot_layout(heights=c(2,0.8))
overall_p

#rhs <- (inset_p + labs(tag="D") )(transmissibility_dist + labs(tag="E"))# + plot_layout(heights=c(1.2,1))
#png("figS1.png",width=10,height=6,res=300,units="in")
#top_p <- (age_dist_p | case_dist_p | transmissibility_dist) + plot_layout(widths = c(1,1,1.5))
#overall_p <- top_p / p_new_fig1_C / (inset_p + theme(legend.position=c(0.55,0.95),legend.justification = 0.5)+ labs(tag="E") ) + plot_layout(heights=c(1,0.9,1.4))
#overall_p
pdf("fig1.pdf",width=6,height=10)
overall_p
dev.off()

#((age_dist_p / case_dist_p) | p_new_fig1_C | rhs) + plot_layout(widths=c(0.8,1,1.8))
png("fig1.png",width=6,height=10,res=300,units="in")
overall_p
dev.off()



