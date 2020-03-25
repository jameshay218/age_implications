library(ggalluvial)
## POLYMOD contact matrix
data("polymod")

river_ages <- read.csv("data/river_ages.csv")
n1 <- sum(river_ages$n)
propns <- river_ages %>% group_by(group1) %>% summarise(group_n=sum(n)/n1)

N_props1 <- propns$group_n
age_dat_short <- data.frame(lower.age.limit=c(0,20,65),population=N_props1)
polymod_c <- contact_matrix(polymod,survey.pop=age_dat_short,age.limits = c(0,20,65),symmetric=TRUE,
                            missing.contact.age = "sample",
                            missing.participant.age = "remove")
C <- polymod_c$matrix
row.names(C) <- colnames(C)

## Solve base model for comparison
N_props <- age_distribution %>% pull(propn)
N_age_classes <- length(N_props1)
N_immunity_classes <- 1

## Number of people in each age group and age class
N <- matrix(N_props*N_tot,ncol=1,nrow=N_age_classes)

## Relative susceptibility of each age group
alphas <- c(1)

R0 <- 2
gamma <- 6
beta_par <- get_beta(C,age_dat_short$population,gamma, R0)

prop_immune <- rep(0, 3)
beta_scales <- c(0.01,0.2,1)
alphas <- c(1)

N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
N_props_long <- N_props
N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)

C_use <- setup_C(C, N, rep(1,3))
wow <- reshape2::melt(C)
par(mfrow=c(1,2))
colnames(wow) <- c("Infectee","Infector","value")

wow$Infector <- factor(wow$Infector,levels=rev(c("[0,20)","[20,65)","65+")))
wow$Infectee <- factor(wow$Infectee,levels=rev(c("[0,20)","[20,65)","65+")))
cols <- c("#E69F00","#56B4E9", "#009E73")

C1 <- C*0.15
C1[,1] <- C1[,1] * 0.01
C1[,2] <- C1[,2] * 0.2
C1[,1] <- C1[,1]*2
wow1 <- reshape2::melt(C1)
colnames(wow1) <- c("Infectee","Infector","value")

         
         

wow$Infector <- factor(wow$Infector,levels=c("[0,20)","[20,65)","65+"))
wow$Infectee <- factor(wow$Infectee,levels=c("[0,20)","[20,65)","65+"))

p5 <- ggplot(wow, aes(axis1=Infector,axis2=Infectee,y=value)) + 
  geom_alluvium(aes(fill=Infector,col=Infector)) + 
  geom_stratum(fill="grey90") + 
  geom_text(stat="stratum",infer.label=TRUE,size=3) +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols)+
  scale_x_continuous(breaks=c(1,2), labels=c("Infector","Infectee"),position="top") +
  theme_void() + 
  theme(legend.position="none",
        text=element_text(size=8),
        axis.text.x=element_text(size=12,vjust=-2))


wow1$Infector <- factor(wow1$Infector,levels=(c("[0,20)","[20,65)","65+")))
wow1$Infectee <- factor(wow1$Infectee,levels=(c("[0,20)","[20,65)","65+")))
p6 <- ggplot(wow1, aes(axis1=Infector,axis2=Infectee,y=value)) + 
  geom_alluvium(aes(fill=Infector,col=Infector)) + 
  geom_stratum(fill="grey90") + 
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols)+
  geom_text(stat="stratum",infer.label=TRUE,size=3) +
  scale_x_continuous(breaks=c(1,2), labels=c("Infector","Infectee"),position="top") +
  theme_void() + 
  theme(legend.position="none",
        text=element_text(size=8),
        axis.text.x=element_text(size=12,vjust=-2))
p6

png("fig2.png",width=8,height=5,res=300, units="in")
(p5 + p6)
dev.off()



pdf("fig2.pdf",width=8,height=5)
(p5 + p6)
dev.off()


pdf("fig2.pdf",width=8,height=9)
((p3 + labs(tag="A") + xlab("") + theme(plot.margin=unit(c(0,0,0,0),"cm"),plot.tag = element_text(size=14))| 
    p1 + labs(tag="B") + xlab("")+ theme(plot.margin=unit(c(0,0,0,0),"cm"),plot.tag = element_text(size=14)))/
    (p4+ labs(tag="C") + theme(plot.margin=unit(c(0,0,0,0),"cm"),plot.tag = element_text(size=14))| 
       p2+ labs(tag="D")+ theme(plot.margin=unit(c(0,0,0,0),"cm"),plot.tag = element_text(size=14)))/
    (p5 + labs(tag="E")+ theme(plot.margin=unit(c(0,0,0,0),"cm"),plot.tag = element_text(size=14)) | 
       p6 + labs(tag="F")+ theme(plot.margin=unit(c(0,0,0,0),"cm"),plot.tag = element_text(size=14)))) + 
  plot_layout(heights=c(1,1,2.5))
dev.off()
png("fig2.png",width=8,height=9,units="in",res=300)
((p3 + labs(tag="A") + xlab("") + theme(plot.margin=unit(c(0,0,0,0),"cm"),plot.tag = element_text(size=14))| 
    p1 + labs(tag="B") + xlab("")+ theme(plot.margin=unit(c(0,0,0,0),"cm"),plot.tag = element_text(size=14)))/
    (p4+ labs(tag="C") + theme(plot.margin=unit(c(0,0,0,0),"cm"),plot.tag = element_text(size=14))| 
       p2+ labs(tag="D")+ theme(plot.margin=unit(c(0,0,0,0),"cm"),plot.tag = element_text(size=14)))/
    (p5 + labs(tag="E")+ theme(plot.margin=unit(c(0,0,0,0),"cm"),plot.tag = element_text(size=14)) | 
       p6 + labs(tag="F")+ theme(plot.margin=unit(c(0,0,0,0),"cm"),plot.tag = element_text(size=14)))) + 
  plot_layout(heights=c(1,1,2.5))
dev.off()
