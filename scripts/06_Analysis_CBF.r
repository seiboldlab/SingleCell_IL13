################ Here I'm plotting ciliary beat frequency (CBF) data from DDM

library(plotrix)
library(beeswarm)
library(openxlsx)
library(reshape2)
library(scales)
library(lme4)
library(lmerTest)
library(openxlsx)


#########################################################################################
#########################################################################################
####### TOP DOWN IMAGING 
#########################################################################################
#########################################################################################


#Bring in raw data from 6 donors
CBF_df<-read.table("../data/CBF_topdown.txt",header=T,stringsAsFactors=F)

#Make dataset that compares washed and no-washed effects
fourWayDF<-data.frame("CBF"=c(CBF_df$CBF[which(CBF_df$Treatment=="BSA" & CBF_df$Wash=="noWash")],
	CBF_df$CBF[which(CBF_df$Treatment=="IL13" & CBF_df$Wash=="noWash")],
	CBF_df$CBF[which(CBF_df$Treatment=="BSA" & CBF_df$Wash=="yesWash")],
	CBF_df$CBF[which(CBF_df$Treatment=="IL13" & CBF_df$Wash=="yesWash")]),
	"comparison"=c(rep("BSA_noWash",length(CBF_df$CBF[which(CBF_df$Treatment=="BSA" & CBF_df$Wash=="noWash")])),
	rep("IL13_noWash",length(CBF_df$CBF[which(CBF_df$Treatment=="IL13" & CBF_df$Wash=="noWash")])),
	rep("BSA_yesWash",length(CBF_df$CBF[which(CBF_df$Treatment=="BSA" & CBF_df$Wash=="yesWash")])),
	rep("IL13_yesWash",length(CBF_df$CBF[which(CBF_df$Treatment=="IL13" & CBF_df$Wash=="yesWash")]))))
fourWayDF$comparison<-factor(fourWayDF$comparison,levels=unique(fourWayDF$comparison))





###############Do density plots
#####Just the non-washed
pdf("Density_CBF_nonWashed.pdf",width=4,height=4)
par(bty="n")
density1<-density(fourWayDF$CBF[which(fourWayDF$comparison=="IL13_noWash")])
plot(density1,las=1,ylim=c(0,0.3),main="CBF shifts due to IL-13",xlab="Ciliary Beat Frequency",type="n")
polygon(density1, col = alpha("red", 0.5), border="black")
density2<-density(fourWayDF$CBF[which(fourWayDF$comparison=="BSA_noWash")])
polygon(density2, col = alpha("midnightBlue", 0.5), border="black")
dev.off()

#####Just the washed
pdf("Density_CBF_washed.pdf",width=4,height=4)
par(bty="n")
density1<-density(fourWayDF$CBF[which(fourWayDF$comparison=="IL13_yesWash")])
plot(density1,las=1,ylim=c(0,0.16),main="CBF shifts due to IL-13",xlab="Ciliary Beat Frequency",type="n")
polygon(density1, col = alpha("orange", 0.5), border="black")
density2<-density(fourWayDF$CBF[which(fourWayDF$comparison=="BSA_yesWash")])
polygon(density2, col = alpha("lightblue", 0.5), border="black")
dev.off()












############################## Now do box plots for these same comparisons

#No washed
pdf("Boxplots_CBF_noWash.pdf",width=1.6,height=3.5)
color=c("blue","red")
#dev.new(width=1.6,height=3.5)
par(bty="n")
boxplot(fourWayDF$CBF[which(fourWayDF$comparison == "BSA_noWash" | 
	fourWayDF$comparison == "IL13_noWash")]~as.character(fourWayDF$comparison[which(fourWayDF$comparison == "BSA_noWash" | 
	fourWayDF$comparison == "IL13_noWash")]),las=1,ylab="Ciliary beat frequency (Hz)",col=color,outpch=16,outcex=0.5)
dev.off()

pdf("Boxplots_CBF_noWash_noOutliers.pdf",width=1.6,height=3.5)
color=c("blue","red")
#dev.new(width=1.6,height=3.5)
par(bty="n")
boxplot(fourWayDF$CBF[which(fourWayDF$comparison == "BSA_noWash" | 
	fourWayDF$comparison == "IL13_noWash")]~as.character(fourWayDF$comparison[which(fourWayDF$comparison == "BSA_noWash" | 
	fourWayDF$comparison == "IL13_noWash")]),las=1,ylab="Ciliary beat frequency (Hz)",col=color,outline=F)
dev.off()


#Washed
pdf("Boxplots_CBF_yesWash.pdf",width=1.6,height=3.5)
color=c("blue","red")
#dev.new(width=1.6,height=3.5)
par(bty="n")
boxplot(fourWayDF$CBF[which(fourWayDF$comparison == "BSA_yesWash" | 
	fourWayDF$comparison == "IL13_yesWash")]~as.character(fourWayDF$comparison[which(fourWayDF$comparison == "BSA_yesWash" | 
	fourWayDF$comparison == "IL13_yesWash")]),las=1,ylab="Ciliary beat frequency (Hz)",col=color,outpch=16,outcex=0.5)
dev.off()

pdf("Boxplots_CBF_yesWash_noOutliers.pdf",width=1.6,height=3.5)
color=c("blue","red")
#dev.new(width=1.6,height=3.5)
par(bty="n")
boxplot(fourWayDF$CBF[which(fourWayDF$comparison == "BSA_yesWash" | 
	fourWayDF$comparison == "IL13_yesWash")]~as.character(fourWayDF$comparison[which(fourWayDF$comparison == "BSA_yesWash" | 
	fourWayDF$comparison == "IL13_yesWash")]),las=1,ylab="Ciliary beat frequency (Hz)",col=color,outline=F)
dev.off()



############################## Mann-Whitney p-values for all comparisons

#### No wash BSA versus No wash IL-13
wilcox.test(x=fourWayDF$CBF[which(fourWayDF$comparison == "BSA_noWash")],
	y=fourWayDF$CBF[which(fourWayDF$comparison == "IL13_noWash")],paired=F,alternative="greater")

#### Wash BSA versus Wash IL-13 (inserts 1 and 2 only)
wilcox.test(x=fourWayDF$CBF[which(fourWayDF$comparison == "BSA_yesWash")],
	y=fourWayDF$CBF[which(fourWayDF$comparison == "IL13_yesWash")],paired=F,alternative="greater")
















############################## Test for interaction effect for IL-13 contigent on washing
#First, fit a linear mixed model with treatment, wash status, and an interaction term as fixed effects, and donor as a random effect
#Then calculate p-values and degrees of freedom using Satterthwaite approximations
CBF_df_1And2<-CBF_df

#response ~ treatment + wash status + donor (random effect) + treatment:wash status
lmer_CBF<-lmer(CBF_df_1And2$CBF~CBF_df_1And2$Treatment + CBF_df_1And2$Wash + (1 | CBF_df_1And2$Donor) + (CBF_df_1And2$Treatment * CBF_df_1And2$Wash))
summary(lmer_CBF)

#Now use the Satterthwaite approximiation to get p-values for coefficents
#Extract coefficients
coefs <- data.frame(coef(summary(lmer_CBF)))
#Get Satterthwaite-approximated degrees of freedom
coefs$df.Satt <- coef(summary(lmer_CBF))[, 3]
#Get Satterthwaite-approximated p-values
coefs$p.Satt <- coef(summary(lmer_CBF))[, 5]


#Now get marginal effect of IL-13 for washed and no washed separately
##Including donor as random effect or fixed effect returns same coefficients
#Pre-wash
##Isolate just inserts 1 and 2
CBF_df_1And2<-CBF_df[which(CBF_df$Wash == "noWash"),]

#Now do regression model
#response ~ treatment + wash status + donor (random effect) + treatment:wash status
lmer_CBF<-lmer(CBF_df_1And2$CBF~CBF_df_1And2$Treatment + (1 | CBF_df_1And2$Donor))
summary(lmer_CBF)
coefs_noWash <- data.frame(coef(summary(lmer_CBF)))
coefs_noWash$df.Satt <- coef(summary(lmer_CBF))[, 3]
coefs_noWash$p.Satt <- coef(summary(lmer_CBF))[, 5]

#Post-wash
##Isolate just inserts 1 and 2
CBF_df_1And2<-CBF_df[which(CBF_df$Wash == "yesWash"),]

#Now do regression model
#response ~ treatment + wash status + donor (random effect) + treatment:wash status
lmer_CBF<-lmer(CBF_df_1And2$CBF~CBF_df_1And2$Treatment + (1 | CBF_df_1And2$Donor))
summary(lmer_CBF)
coefs_yesWash <- data.frame(coef(summary(lmer_CBF)))
coefs_yesWash$df.Satt <- coef(summary(lmer_CBF))[, 3]
coefs_yesWash$p.Satt <- coef(summary(lmer_CBF))[, 5]






















#########################################################################################
#########################################################################################
####### PROFILE IMAGING 
#########################################################################################
#########################################################################################

#Bring in raw data from 2 donors
CBF_profile<-read.table("../data/CBF_profile.txt",header=T,stringsAsFactors=F)
CBF_profile$Treat_Wash<-paste(CBF_profile$Treatment,CBF_profile$Wash,sep="_")


#Make density plots
pdf("Density_CBF_profile_noWash_FIGURE.pdf",width=4,height=4)
#dev.new(width=4,height=4)
par(bty="n")
density1<-density(CBF_profile$CBF[which(CBF_profile$Treat_Wash=="IL13_noWash")])
density2<-density(CBF_profile$CBF[which(CBF_profile$Treat_Wash=="BSA_noWash")])
plot(density1,las=1,xlim=c(0,max(c(max(density1$x),max(density2$x)))),ylim=c(0,0.3),
	main="CBF shifts due to IL-13 before wash",xlab="Ciliary Beat Frequency",type="n")
polygon(density2, col = alpha("grey", 0.65), border="black")
polygon(density1, col = alpha("black", 0.65), border="black")
dev.off()

pdf("Density_CBF_profile_yesWash_FIGURE.pdf",width=4,height=4)
#dev.new(width=4,height=4)
par(bty="n")
density3<-density(CBF_profile$CBF[which(CBF_profile$Treat_Wash=="IL13_yesWash")])
density4<-density(CBF_profile$CBF[which(CBF_profile$Treat_Wash=="BSA_yesWash")])
plot(density1,las=1,xlim=c(0,max(c(max(density1$x),max(density2$x)))),ylim=c(0,0.25),
	main="CBF shifts due to IL-13 after wash",xlab="Ciliary Beat Frequency",type="n")
polygon(density4, col = alpha("grey", 0.65), border="black")
polygon(density3, col = alpha("black", 0.65), border="black")
dev.off()





#Make box plots
#Plot without outliers
pdf("Boxplots_CBF_profile_noWash_noOutliers_FIGURE.pdf",width=1.65,height=4)
color=c("grey","black")
#dev.new(width=1.65,height=4)
par(bty="n")
boxplot(CBF_profile$CBF[which(CBF_profile$Wash=="noWash")]~CBF_profile[which(CBF_profile$Wash=="noWash"),]$Treatment,
	las=1,ylab="Ciliary beat frequency (Hz)",col=color,outpch=16,outcex=0.5,medcol="white",
	names=F,outline=F,main="Before wash",ylim=c(0,15))
dev.off()

pdf("Boxplots_CBF_profile_yesWash_noOutliers_FIGURE.pdf",width=1.65,height=4)
color=c("grey","black")
#dev.new(width=1.65,height=4)
par(bty="n")
boxplot(CBF_profile$CBF[which(CBF_profile$Wash=="yesWash")]~CBF_profile[which(CBF_profile$Wash=="yesWash"),]$Treatment,
	las=1,ylab="Ciliary beat frequency (Hz)",col=color,outpch=16,outcex=0.5,medcol="white",
	names=F,outline=F,main="After wash",ylim=c(0,15))
dev.off()








####Test for interaction effect for IL-13 contigent on washing
#response ~ treatment + wash status + donor (random effect) + treatment:wash status
lmer_CBF_profile<-lmer(CBF_profile$CBF~CBF_profile$Treatment + CBF_profile$Wash + (1 | CBF_profile$Donor) + (CBF_profile$Treatment * CBF_profile$Wash))
summary(lmer_CBF_profile)

#Now use the Satterthwaite approximiation to get p-values for coefficents
#Extract coefficients
coefs <- data.frame(coef(summary(lmer_CBF_profile)))
#Get Satterthwaite-approximated degrees of freedom
coefs$df.Satt <- coef(summary(lmer_CBF_profile))[, 3]
#Get Satterthwaite-approximated p-values
coefs$p.Satt <- coef(summary(lmer_CBF_profile))[, 5]



#Now get marginal effect of IL-13 for washed and no washed separately
##Including donor as random effect or fixed effect returns same coefficients

#Pre-wash
CBF_profile_noWash<-CBF_profile[which(CBF_profile$Wash == "noWash"),]
#Now do regression model
lmer_CBF_noWash<-lmer(CBF_profile_noWash$CBF~CBF_profile_noWash$Treatment + (1 | CBF_profile_noWash$Donor))
summary(lmer_CBF_noWash)
coefs_noWash <- data.frame(coef(summary(lmer_CBF_noWash)))
coefs_noWash$df.Satt <- coef(summary(lmer_CBF_noWash))[, 3]
coefs_noWash$p.Satt <- coef(summary(lmer_CBF_noWash))[, 5]


#Post-wash
CBF_profile_yesWash<-CBF_profile[which(CBF_profile$Wash == "yesWash"),]
#Now do regression model
lmer_CBF_yesWash<-lmer(CBF_profile_yesWash$CBF~CBF_profile_yesWash$Treatment + (1 | CBF_profile_yesWash$Donor))
summary(lmer_CBF_yesWash)
coefs_yesWash <- data.frame(coef(summary(lmer_CBF_yesWash)))
coefs_yesWash$df.Satt <- coef(summary(lmer_CBF_yesWash))[, 3]
coefs_yesWash$p.Satt <- coef(summary(lmer_CBF_yesWash))[, 5]






