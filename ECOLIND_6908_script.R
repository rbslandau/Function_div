# Script.R: R-Code to reproduce the results of
# Voß, K.and Schäfer, R.B.
# Taxonomic and functional diversity of stream invertebrates along an environmental stress gradient
# Ecological Indicators
# Katharina Voß and Ralf B. Schäfer
# University Koblenz-Landau
# Institute for Environmental Sciences
# Fortstr. 7
# 76829 Landau
# GERMANY
# voss@uni-landau.de
# schaefer-ralf@uni-landau.de

#######################################################################################
#
#     The R code below was tested on R version 3.2.2 on Windows 7 and OS X 10.12
#
########################################################################################

# Code written by K Voß, revised by RB Schäfer
## The code has been written for application to the supplied data sets

# Structure of the code:
# 1 # Gradient building ---------------------------------------------> from line  54 - 135
# 2 # Diversity metrics and gradient --------------------------------> from line 135 - 242
# 3 # Diversity metrics and Organic matter breakdown-----------------> from line 242 - 290
# 4 # Stressor effect on invertertebrate community structure (RDA)---> from line 290 - 341
# 5 # Stressor effect on community trait structure (RDA) ------------> from line 341 - 384
# Fig.1 and Fig.2 from the article ----------------------------------> from line 384 - 444
# Fig. S1 and S2 from Supplementary Informations---------------------> from line 444 - 466

# See README.txt for details
rm(list = ls())

# Libraries Used
pkg <- c("vegan", "pcaPP", "BiodiversityR", "FD")
# Install the packages if missing
# install.packages(pkg)
# Call packages
lapply(pkg, library, character.only=T)

##################
### Setup project
##################

# Set path. Change this line to where you have saved the data files.
# change slashs from \ into /, for example on Windows OS
prj <- "C:/Users/KP/Dropbox/3.paper_revisions/data"
setwd(prj)
# load data
parameter <- read.csv(file.path(prj, "variables_final.csv"), sep=";", dec=".", header=T, check.names=F) 

####################
# 1. Build gradient
####################

rownames(parameter) <- parameter[ ,1]
para <- parameter[ ,-c(1:3)]

# check for normal distribution of continuous variables and visual overview
# spm(para)
# skewness - known from Voss et al. 2015
# double sqrt and log transformation improves distribution

para$"Sqrt(Velocity)" <- sqrt(para$Velocity)
para$"log10(EC)" <- log10(para$EC)
para$"log10(NO2)" <- log10(para$NO2)
para$"log10(NO3)" <- log10(para$NO3)
para$"log10(PO4)" <- log10(para$PO4)
para$"log10(Width riparian zone)" <- log10(para$"Width of riparian zone")

#check for collinearity (omit in case of r > 0.7)
round(cor(para),2)

#   omit variables:
#   A) Skewed variables: Velocity, EC, NO2, NO3, PO4, Width.of.riparian.zone (keep the transformed variables!)
#   B) log_EC (because it strongly correlated with pH (r = 0.74), O2 (r = 0.73), and NO3 (r = 0.79))
#   C) FPOM (all values = 0)

new_par <- para[ ,-c(3,7:10,26,29,31)]

#### Preparation: Values of broken stick model
np_pca <- rda(new_par, scale=TRUE)
library(BiodiversityR)
bs_temp <- PCAsignificance(np_pca)
# in classical PCA 6 axes hold more variance than BS model

#### Conduct sparse PCA
# we estimate the optimal penalty term first for 6 axes
k.max <- 6
##  k.max = max number of considered sparse PCs
oTPO <- opt.TPO(scale(new_par), k.max = k.max, method = "sd")

oTPO$pc
# the model selected by opt. TPO
oTPO$pc$load     
# and related sparse loadings

##  Tradeoff Curves: Explained Variance vs. sparseness
par (mfrow = c (1, k.max))
for (i in 1:k.max)        plot (oTPO, k = i)
# L0 gives number of variables with zero loadings on PC
# lambda opt is the estimated optimal penalty
# ECV1 is the empirical cumulated variance of the PCs

##  Tradeoff Curves: Explained Variance vs. lambda
par (mfrow = c (1, k.max))
for (i in 1:k.max)       plot (oTPO, k = i, f.x = "lambda")
# lambda represents penalty term

# use optimized lambdas to compute sparsePCA
spc <- sPCAgrid(scale(new_par), k = k.max, lambda = oTPO$pc.noord$lambda, method = "sd")
summary(spc)
# first axis captures 40% of variance

# loadings for Table S4
unclass (spc$load)
# first axis captures major part of sparse PCA variance
# and represents most variables

##############################
# Extraction of spc scores    #
##############################
load_spca <- scores(spc, choices = 1, display = "species", scaling = 0)
# write.csv(load_spca, file="Tab_S4.csv")
# uncomment to save loadings (given in Table S4)

# scores will be used as environmental stress gradient in further analysis
pca_axes <- vegan::scores(spc, disp = "sites", choices=c(1:ncol(new_par)), scaling = "sites")
# we extract the first axis
grad <- unlist(pca_axes[ ,1])
# write.csv(pca_axes[,1], file ="stressor gradient")

#####################################################
# 2.    Diversity metrics and gradient
####################################################

# load data on taxonomic diversity metrics and relative abundance
# of each trait modality within the communities
data <- read.csv(file.path(prj, "diversities.csv"), sep="\t", dec=".", header=T, check.names=F)
str(data)

# relationship between taxonomic diversities (Total taxonomic richness and Simpsons Diversity) and gradient

########################
# Taxonomic diversity
########################
# linear regression for taxonomic diversity
m1 <- cor.test(data$TTR, grad)
m2 <- cor.test(data$SD, grad)

p.adjust(c(m1$p.value, m2$p.value), method= "fdr")

# for plot refer to section 5

###########################################################################################################
# Calculation of FD: functional richness, evenness, divergence over all trait modalities for each community
# and for OMB-relevant traits seperately for each community
# to relate it to gradient with simple regression
###########################################################################################################

# load data
trait<- read.csv(file.path(prj, "ID_traits.csv"), sep= "", header=T) 
abun <- read.csv(file.path(prj, "ID_abundance.csv"), sep= "", header=T) 

# calculating proportions of single trait modalities per trait for the fuzzy coded data
col.blocks <- c(7,2,3,4,8,4,5,5,8,9,8,7,8,3,9,4,3,2,3,5,6)
prop_trait <- prep.fuzzy.var(trait, col.blocks, row.w = rep(1, nrow(trait)))
# write.csv(prep.fuzzy.var(trait, col.blocks, row.w = rep(1, nrow(trait))))#, file = "prop_traits.csv")

# calculating distance-based functional diversity indices
fds <- dbFD(x = prop_trait, a = abun, stand.FRic = TRUE)
str(fds)
data$overall_FRic <- fds$FRic
data$overall_FEve <- fds$FEve
data$overall_FDiv <- fds$FDiv

##############################################################################
## Functional diversity of OMB-related traits (n=2, food and feeding habit)

fds_Food <- dbFD(x = prop_trait[, c(47:55)], a = abun, stand.FRic = TRUE)
str(fds_Food)

data$Food_FRic <- fds_Food$FRic
data$Food_FEve <- fds_Food$FEve
data$Food_FDiv <- fds_Food$FDiv

fds_Feeding_habits <- dbFD(x = prop_trait[, c(56:63)], a = abun, stand.FRic = TRUE)

data$FH_FRic <- fds_Feeding_habits$FRic
data$FH_FEve <- fds_Feeding_habits$FEve
data$FH_FDiv <- fds_Feeding_habits$FDiv

########################
#Functional Diversities
########################

#extract FD indices
FD <- data[ ,120:128]

str(FD)
res <- sapply(FD, function(y) {
  test <- cor.test(y, grad, method="pearson")
})

# extract coef. of correlation
cor <- res[4, ]
value <- res[3, ]
p.adjusted <- p.adjust(value, method= "fdr")
fin_FD <- rbind(unlist(cor), unlist(p.adjusted))
fin_FD2 <- as.data.frame(t(fin_FD))
fin_FD2$"FD_metric" <- row.names(fin_FD2)
names(fin_FD2)[1:2] <- c("Correl_coeff", "p_adjusted")
# save (cf. Tab.S5)
# write.csv(fin_FD2, file="cor_FD_gradient.csv", row.names=F)

######################################################################
# Relative abundance of trait modalities (n=113) and gradient
######################################################################

# extract trait modalities
modalities <- data[ ,c(7:119)]

str(modalities)

res1 <- sapply(modalities, function(y) {
  test1 <- cor.test(y, grad, method="pearson")
})

# extract coef. of correlation
cor1 <- res1[4,]
value1 <- res1[3,]
p.adjusted_2 <- p.adjust(value1,method= "fdr")
fin_modal <- rbind(unlist(cor1), unlist(p.adjusted_2))
fin_modal2 <- as.data.frame(t(fin_modal))
fin_modal2$"Modality" <- row.names(fin_modal2)
names(fin_modal2)[1:2] <- c("Correl_coeff", "p")
# save (cf. Tab. S6)
# write.csv(fin_modal2, file="cor_modalities_gradient.csv", row.names = F)

############################################################
## 3.  Diversity metrics and Organic matter breakdown (OMB)
############################################################

# Analyse whether changes in taxonomic and functional diversity propagate to OMB
# relationship between taxonomic diversities and Organic matter breakdown
# relationship between functional diversities (two OMB-relevant FDs and their modalities) and OMB
# simple regression with pearson correlation

# Taxonomic diversity with OMB
m3 <- cor.test(data$TTR, data$OMB)
m4 <- cor.test(data$SD, data$OMB)
p.adjust(c(m3$p.value, m4$p.value), method= "fdr")

#for plot refer to section 5

############################################
#Functional Diversities and OMB
############################################

res_fd <- sapply(FD[ ,-c(1:3)], function(y) {
  test <- cor.test(y, data$OMB, method="pearson")
})

# extract coef. of correlation
cor_fd <- res_fd[4, ]
value_fd <- res_fd[3, ]
p.adjusted2 <- p.adjust(value_fd, method= "fdr")

#food and feeding habit modalities (n=17) and OMB
OMB_mod<- modalities[,47:63]

res2 <- sapply(OMB_mod, function(y) {
  test <- cor.test(y, data$OMB, method="pearson")
})

# extract coef. of correlation
cor2 <- res2[4, ]
value2 <- res2[3, ]
p.ad2 <- p.adjust(value2, method= "fdr")
fin_OMB <- rbind(unlist(cor2), unlist(p.ad2))
fin_OMB2 <- as.data.frame(t(fin_OMB))
fin_OMB2$"OMBrelevant_mod" <- row.names(fin_OMB2)
names(fin_OMB2)[1:2] <- c("Correl_coeff", "p")

# save (cf. Tab.S7)
# write.csv(fin_OMB2, file="cor_OMBmod_OMB.csv", row.names=F)

#############################################################
# 4.  Redundancy Discriminant Analysis for invertebrate taxa#
#############################################################
taxa <- read.csv(file.path(prj, "inv_taxa_rda.csv"), sep="\t", header=T, check.names=F)
rownames(taxa) <- taxa[ ,1]
aggr_taxa <- taxa[-1, ]
#omit ID's

# aggregate some species with same trait information and where species are rare 
# (Cordulegaster sp., Philopotamus sp., Ptychoptera sp and Rhyacophila sp.)
# through summing their abundances

str(aggr_taxa)
aggr_taxa[1:67] <- lapply(aggr_taxa[1:67], as.numeric)

aggr_taxa$"Cordulegaster spec." <- rowSums(aggr_taxa[ ,c(10,11)])
aggr_taxa$"Philopotamus spec." <- rowSums(aggr_taxa[ ,c(38:40)])
aggr_taxa$"Ptychoptera spec." <- rowSums(aggr_taxa[ ,c(48:49)])
aggr_taxa$"Rhyacophila spec." <- rowSums(aggr_taxa[ ,c(51:56)])

# remove aggregated taxa and site information
varetaxa <- aggr_taxa[ ,-c(1,10,11,38:40,48:49,51:56)]

# We employ the hellinger transformation to use RDA with ecological data 
# see Legendre, P., Gallagher, E.D., 2001. Ecologically meaningful transformations for ordination of species data. Oecologia 129, 271–280.

varetaxa_hel <- decostand(varetaxa, "hellinger")

#############
#conduct RDA#
#############

(VAR_rda <- rda(varetaxa_hel ~ grad, scale=FALSE))
# we obtain 1 RDA (constrained) and 27 PCA (unconstrained) axes
# 7% explained by the stressor gradient

# Extraction of canonical coefficients from rda object
coef(VAR_rda)

# Global test of the RDA result
set.seed(111)
anova.cca(VAR_rda, step=1000, by='term')
# grad statistically significant 0.038

# order of species along the gradient (RDA 1 axis - constrained variance = 8%)
sc  <- vegan::scores(VAR_rda, display = 'species', choices = 1, scaling = "species")
hist(sc) #most species in the centre, will not be plotted
par(mfrow=c(1,1))
linestack(sc[abs(sc) > 0.07], rownames(sc)[abs(sc) > 0.07], axis = TRUE, cex = 1.5)
#rownames includes species names

#################################################################
# 5.   Redundancy Discriminant Analysis for invertebrate traits #
#################################################################

#RDA with trait modality data
modalities_hel <- decostand(modalities, "hellinger") # transformation

#############
#conduct RDA#
#############

(mod_rda <- rda(modalities_hel ~ grad))
# 13% explained by stressor gradient (constrained Proportion 0.127)

coef(mod_rda)

# Global test of the RDA result
set.seed(111)
anova.cca(mod_rda, step=1000) #model statistically significant 

#order of modalities along stressor gradient (RDA 1 axis - constrained variance = 13%)
sc1  <- scores(mod_rda, display = 'species', choices = 1)

linestack(sc1[abs(sc1) > .05], rownames(sc1)[abs(sc1) > .05], axis = TRUE, cex = 1.5) 

##############################################
#conduct RDA for OMB-relevant trait modalities

OMBmod_hel <- decostand(OMB_mod, "hellinger") # transformation

(OMBmod_rda <- rda(OMBmod_hel ~ grad))
# 8% explained variance
coef(OMBmod_rda)

set.seed(111)
anova.cca(OMBmod_rda, step=1000)
# p = 0.055

# order of OMB-relevant trait modalities along stressor gradient (RDA 1 axis - constrained variance = 13%)
sc2  <- scores(OMBmod_rda, display = 'species', choices = 1, scaling ="species")
par(mai=c(0.2,0,0.5,2), cex= 1.5)
linestack(sc2[abs(sc2) > .01], rownames(sc2)[abs(sc2) > .01], axis = TRUE) 

######################################
# Code for figures in article
#####################################

#Figure 1
#jpeg('Figure1.jpg', width = 500, height = 150,units = 'mm',res = 600)
par(mfrow=c(1,4), mai=c(0.8, 0.3, 0.5, 3.6), cex= 1.5)
#a)gradient
# account for negative relationship with RDA axes - multiply with -1 to have matching directions
load_spca_n <- load_spca*(-1)
linestack(load_spca_n[abs(load_spca_n) > 0.00003], rownames(load_spca_n)[abs(load_spca_n) > 0.00003], axis = TRUE, cex =1.1)
#b)invertebrates taxa
linestack(sc[abs(sc) > 0.07], rownames(sc)[abs(sc) > 0.07], axis = TRUE, cex =1.1)#rownames includes species names
#c)trait modalities
linestack(sc1[abs(sc1)>.04], rownames(sc1)[abs(sc1)>.04], axis = TRUE, cex =1.1)
#d)OMB-relevant trait modalities
linestack(sc2[abs(sc2) > .01], rownames(sc2)[abs(sc2) > .01], axis = TRUE, cex =1.1)
#dev.off()

#########
#Figure 2
# see above for p adjustments

#jpeg("Figure2a-d.jpg",width = 300, height = 300,units = 'mm',res = 600)
par(mfrow=c(2,2),mar=c(4,5,2,2), cex = 1.5)

# Fig.2a
plot(grad, data$TTR,ylab= "Invertebrate richness", xlab = "Environmental stress gradient", pch=16, yaxt="n")
axis(2,las=2)
# cor.test(data$TTR, grad, method="pearson")
text(-2,15.5, label=expression(paste("r = 0.55; p = 0.004")), cex=0.85)
mtext(text = expression(bold(a)), side = 2, cex=2, las = 1, at = 17, line = 4)
abline(lm(data$TTR~grad), lwd = 1.7)

# Fig 2b
plot(data$SD,data$OMB*1000,yaxt="n", pch = 16, ylab= expression(paste(italic("k"), " invertebrates (",  10^-3, dday^-1,")")), xlab = "Invertebrate Simpson Diversity")
axis(2,las=2)
# cor.test(data$SD,data$OMB,method = c("pearson"))
text(2.5,5, label=expression(paste("r = -0.41, p = 0.06")), cex=0.85)  
mtext(text = expression(bold(b)), side = 2, cex=2, las =1, at = 6.3, line = 4)
abline(lm((data$OMB*1000) ~ data$SD), lwd = 1.7)

# Fig.2c
plot(FD$FH_FEve, data$OMB*1000, ylab= expression(paste(italic("k"), " invertebrates (",  10^-3, dday^-1,")")), xlab = "Functional Evenness - Feeding habit", pch=16, yaxt="n")
axis(2,las=2)
#cor.test(FD$FH_FEve,data$OMB, method="pearson")
text(0.35,5.6, label=expression(paste("r = -0.5; p = 0.03")), cex=0.85)
mtext(text = expression(bold(c)), side = 2, cex=2, las =1, at = 6.3, line = 4)
abline(lm(data$OMB*1000 ~FD$FH_FEve), lwd = 1.7)

#Fig 2d
plot(data$SD, FD$FH_FEve, ylab= "Functional Evenness - Feeding habit", xlab = "Invertebrate Simpson Diversity", pch=16,yaxt="n")
axis(2,las=2)
# cor.test(FD$FH_FEve,data$SD, method="pearson")
text(2.5,0.8, label=expression(paste("r = 0.44; p = 0.02")), cex=0.85)
mtext(text = expression(bold(d)), side = 2, cex=2, las =1, at = 0.9, line = 4)
abline(lm(FD$FH_FEve~data$SD), lwd = 1.7)

dev.off()

##########################################################################
#Supplementary information
############################

#jpeg("Fig. S1.jpg",width = 230, height = 230,units = 'mm',res = 600)

par(mfrow=c(2,2),mar=c(4,5,2,2), cex = 1.5)

#Fig. S1a
plot(data$Gam.Ab, data$SD, yaxt="n", pch = 16, ylab= "Invertebrate Simpson Diversity", xlab ="Number of Gammarids")
axis(2, las=2)
cor.test(data$SD, data$Gam.Ab, method = c("pearson"))
text(200,5.2, label=expression(paste("r = -0.37, p = 0.045")), cex=1)  
mtext(text = expression(bold(a)), side = 2, cex=2, las = 1, at = 6.2, line = 4)
abline(lm(data$SD~data$Gam.Ab), lwd = 1.7)

#Fig. S1b
plot(data$Gam.Ab, data$OMB*1000, yaxt="n", pch = 16, ylab= expression(paste(italic("k"), " invertebrates (",  10^-3, dday^-1,")")), xlab = "Number of Gammarids")
axis(2,las=2)
cor.test(data$Gam.Ab, data$OMB*1000, method = c("pearson"))
text(200, 5.2, label=expression(paste("r = 0.48, p = 0.01")), cex=1)  
mtext(text = expression(bold(b)), side = 2, cex=2, las = 1, at = 6.4, line = 4)
abline(lm(data$OMB*1000~data$Gam.Ab), lwd = 1.7)

#Fig.S1c
plot(data$Gam.Ab, FD$FH_FEve, yaxt="n", pch = 16, ylab= "Funct. Evenness - Feeding habit", xlab = "Number of Gammarids")
axis(2, las=2)
cor.test(data$Gam.Ab, FD$FH_FEve, method = c("pearson"))
text(450, 0.7, label=expression(paste("r = -0.5, p = 0.006")), cex=1)  
mtext(text = expression(bold(c)), side = 2, cex=2, las = 1, at = 0.92, line = 4)
abline(lm(FD$FH_FEve~data$Gam.Ab), lwd = 1.7)

#dev.off()

########
#jpeg("Fig.S2.jpg",width = 115, height = 115,units = 'mm',res = 600)

par(mfrow=c(1,1),mar=c(5,5,4,2), cex = 1.5)

plot(grad,data$rel.Gam.Abu,yaxt="n", pch = 16, ylab= "Relative Gammarid abundance", xlab ="Environmental stress gradient")
axis(2,las=2)

#dev.off()