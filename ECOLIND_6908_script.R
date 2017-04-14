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

# Code written by K Voß, in major parts rewritten by RB Schäfer
## The code has been written for application to the supplied data sets

# Structure of the code:
# 1 # Gradient building ---------------------------------------------
# 2 # Diversity metrics and gradient --------------------------------
# 3 # Stressor effect on invertertebrate community structure (RDA)---
# 4 # Stressor effect on community trait structure (RDA) ------------
# 5 # Diversity metrics and Organic matter breakdown-----------------
# Fig.1 and Fig.2 from the article ----------------------------------
# Fig. S1 and S2 from Supplementary Informations---------------------
# See README.txt for details

# Set path to working directory, where files and history is stored
prj <- "~/Desktop/ECOLIND-6908_Voss_Schaefer/data"
# you need to replace this line with a path on your computer
setwd(prj)

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

# load variables data
parameter <- read.csv("https://raw.githubusercontent.com/rbslandau/Function_div/master/variables_final.csv", sep=";", dec=".", header=T, check.names=F) 

# load trait data for diversity indices (trait - species matrix)
trait <- read.csv("https://raw.githubusercontent.com/rbslandau/Function_div/master/ID_traits.csv", sep= "", header=T, as.is =TRUE) 

# read trait labels for later labelling as function for calculation of FD requires neutral labels
trait_labels <- read.csv("https://raw.githubusercontent.com/rbslandau/Function_div/master/Labels_traits.csv", sep="\t", header=T, check.names=F)

# load taxonomic data for diversity indices (species - site matrix)
abun <- read.csv("https://raw.githubusercontent.com/rbslandau/Function_div/master/ID_abundance.csv", sep= "", header=T) 

# load data on organic matter breakdown
omb_dat <- read.csv("https://raw.githubusercontent.com/rbslandau/Function_div/master/omb.csv", sep= "", header=T, as.is =TRUE)

# load taxonomic data for RDA (species - site matrix, but with real species names)
taxa <- read.csv("https://raw.githubusercontent.com/rbslandau/Function_div/master/inv_taxa_rda.csv", sep="\t", header=T, check.names=F)
# assign site code to rownames
rownames(taxa) <- taxa[ ,1]
# remove ID row and site information
taxa_new <- taxa[-1, -1]

str(taxa_new)
# convert to numeric variable
taxa_new[1:67] <- lapply(taxa_new[1:67], as.numeric)

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
bs_temp <- PCAsignificance(np_pca)
# in classical PCA 6 axes hold more variance than BS model

#### Conduct sparse PCA
# we estimate the optimal penalty term first for 6 axes
k.max <- 6
##  k.max = max number of considered sparse PCs
oTPO <- opt.TPO(scale(new_par), k.max = k.max, method = "sd")

summary(oTPO$pc)
# the model selected by opt. TPO
# first axis captures 41% of variance

oTPO$pc$load*-1    
# and represents most variables
# related sparse loadings for Table S4

##############################
# Extraction of spc scores    #
##############################
# use optimized lambdas to compute sparsePCA
spc <- sPCAgrid(scale(new_par), k = k.max, lambda = oTPO$pc.noord$lambda, method = "sd")

load_spca <- scores(spc, choices = 1, display = "species", scaling = 0)
# write.csv(load_spca, file="Tab_S4.csv")
# uncomment to save loadings (given in Table S4)

# scores will be used as environmental stress gradient in further analysis
pca_axes <- vegan::scores(spc, disp = "sites", choices=c(1:ncol(new_par)), scaling = "sites")
# we extract the first axis, multiplied with -1 to enhance interpretability, see paper
grad <- unlist(pca_axes[ ,1]) * -1
# write.csv(pca_axes[,1], file ="stressor gradient")

#####################################################
# 2.    Diversity metrics and gradient
####################################################

# calculate TTR and SD
TTR <- specnumber(abun, MARGIN = 1)
SD <- diversity(abun, index = "invsimpson", MARGIN = 1)
# relationship between taxonomic diversities (Total taxonomic richness and Simpsons Diversity) and gradient

########################
# Taxonomic diversity
########################
# correlation for taxonomic diversity
m1 <- cor.test(TTR, grad)
m2 <- cor.test(SD, grad)

p.adjust(c(m1$p.value, m2$p.value), method= "fdr")
# strong response of TTR compared to SD suggests that mainly rare species were lost

# recalculate SD downweighing dominant taxa
SD_dw <- diversity(sqrt(abun), index = "invsimpson", MARGIN = 1)

(m2_1 <- cor.test(SD_dw, grad))

# identify dominant species
# use taxa data because contains species labels
tot_abun_spec <- colSums(taxa_new)
sort(tot_abun_spec)
# gammarids dominate abundance

# compute fraction of gammarids per site
tot_gam <- rowSums(taxa_new[ ,names(taxa_new) %in% c("Gammarus roeseli", "Gammarus fossarum", "Gammarus pulex")])
tot_abun <- rowSums(taxa_new)
frac_gam <- tot_gam/tot_abun
# gammarids abundance > 50% in 75% of sites, median: 79%

###########################################################################################################
# Calculation of FD: functional richness, evenness, divergence over all trait modalities
# and over OMB-relevant traits for each community
# to check for correlation
###########################################################################################################

# calculating proportions of single trait modalities per trait for the fuzzy coded data
col.blocks <- c(7,2,3,4,8,4, 5,5,8,9,8,7,8,3,9,4,3,2,3,5,6)
prop_trait <- prep.fuzzy.var(trait, col.blocks, row.w = rep(1, nrow(trait)))
# write.csv(prep.fuzzy.var(trait, col.blocks, row.w = rep(1, nrow(trait))))#, file = "prop_traits.csv")

# calculating distance-based functional diversity indices
fds <- dbFD(x = prop_trait, a = abun, stand.FRic = TRUE, w.abun = TRUE)
str(fds)
FD <- NULL
FD$overall_FRic <- fds$FRic
FD$overall_FEve <- fds$FEve
FD$overall_FDiv <- fds$FDiv
# uncomment for calculating functional redundancy according to Mouillotet al. 2014. Proceedings of the National Academy of Sciences 111, 13757–13762. doi:10.1073/pnas.1317625111
# FD$overall_FRed <- fds$sing.sp/fds$nbsp
  
##############################################################################
## Functional diversity only for OMB-related traits (food and feeding habit)

fds_Food <- dbFD(x = prop_trait[ , c(47:63)], a = abun, stand.FRic = TRUE, w.abun = TRUE)
str(fds_Food)

FD$OMB_FRic <- fds_Food$FRic
FD$OMB_FEve <- fds_Food$FEve
FD$OMB_FDiv <- fds_Food$FDiv
# uncomment for calculating functional redundancy
# FD$OMB_FRed <- fds_Food$sing.sp/fds_Food$nbsp

########################
#Functional Diversities
########################
# correlations
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

fin_FD2
# no significant correlation
# although moderate correlation of functional richness with stress gradient

######################################################################
# Check for correlation with individual modalities
######################################################################
# some modalities have been converted to factors
# because all values were 0
# replace with numeric 0
modalities <- fds$CWM
temp_facmod <- sapply(modalities,is.factor)
modalities[ ,temp_facmod] <- 0

# extract the two modalities
ma <- cor.test(modalities$t10.3, grad)
mb <- cor.test(modalities$t11.3, grad)

p.adjust(c(ma$p.value, mb$p.value), method= "fdr")
# dead plant as food -0.18
# being shredder -0.21

plot(grad, modalities$t11.3)
# no relationship

# May be both modalities are relatively constant?
# summary of variation in both OMB-related traits
summary(modalities$t10.3)
summary(modalities$t11.3)
# mean and standard deviation
mean(fds$CWM$t10.3)
sd(fds$CWM$t10.3) 
# mean and standard deviation
mean(fds$CWM$t11.3)
sd(fds$CWM$t11.3) 

#############################################################
# 3.  Redundancy Discriminant Analysis for invertebrate taxa#
#############################################################

# aggregate some species with same trait information and where species are rare
# (Cordulegaster sp., Philopotamus sp., Ptychoptera sp and Rhyacophila sp.)
# through summing their abundances, though influence on results is negligible
pos_cord <- grep("Cordulegaster", names(taxa_new))
pos_phil <- grep("Philopotamus", names(taxa_new))
pos_ptych <- grep("Ptychoptera", names(taxa_new))
pos_rhya <- grep("Rhyacophila", names(taxa_new))

taxa_new$"Cordulegaster spec." <- rowSums(taxa_new[ ,pos_cord])
taxa_new$"Philopotamus spec." <- rowSums(taxa_new[ ,pos_phil])
taxa_new$"Ptychoptera spec." <- rowSums(taxa_new[ ,pos_ptych])
taxa_new$"Rhyacophila spec." <- rowSums(taxa_new[ ,pos_rhya])

# remove aggregated taxa and site information
varetaxa <- taxa_new[ ,-c(pos_cord, pos_phil, pos_ptych, pos_rhya)]
varetaxa <- aggr_taxa
# We employ the hellinger transformation to use RDA with ecological data 
# see Legendre, P., Gallagher, E.D., 2001. Ecologically meaningful transformations for ordination of species data. Oecologia 129, 271–280.

varetaxa_hel <- decostand(varetaxa, "hellinger")

#############
#conduct RDA#
#############

(VAR_rda <- rda(varetaxa_hel ~ grad, scale=FALSE))
# we obtain 1 RDA (constrained) and 27 PCA (unconstrained) axes
# 7% explained by the stressor gradient

# Global test of the RDA result
set.seed(111)
anova.cca(VAR_rda, step=1000, by='term')
# grad statistically significant 0.039

# order of species along the gradient (RDA 1 axis - constrained variance = 7%)
sc  <- vegan::scores(VAR_rda, display = 'species', choices = 1, scaling = "species")
hist(sc) #most species in the centre, will not be plotted
par(mfrow=c(1,1))
linestack(sc[abs(sc) > 0.07], rownames(sc)[abs(sc) > 0.07], axis = TRUE, cex = 1.5)
#rownames includes species names

#################################################################
# 4.   Redundancy Discriminant Analysis for invertebrate traits #
#################################################################

# RDA with trait modality data
modalities_hel <- decostand(modalities, "hellinger") # transformation

##############
#conduct RDA #
##############

(mod_rda <- rda(modalities_hel ~ grad))
# 8% explained by stressor gradient

# Global test of the RDA result
set.seed(111)
anova.cca(mod_rda, step=1000, by='term') 
# model not statistically significant, p = 0.067 - although slightly more variance explained

# order of modalities along stressor gradient (RDA 1 axis - constrained variance = 8%)
sc1  <- scores(mod_rda, display = 'species', choices = 1)
linestack(sc1[abs(sc1) > .041], trait_labels[abs(sc1) > .041], axis = TRUE, cex = 1.5) 

############################################################
## 5.  Diversity metrics and Organic matter breakdown (OMB)
############################################################

# Analyse relationship between taxonomic and functional diversity and OMB

# Taxonomic diversity with OMB
m3 <- cor.test(TTR, omb_dat$OMB)
m4 <- cor.test(SD, omb_dat$OMB)
p.adjust(c(m3$p.value, m4$p.value), method= "fdr")
# not significant, but SD p = 0.06

# Gammarid abundance might explain breakdown
cor.test(tot_gam, omb_dat$OMB)
# moderate correlation, r = 0.48, p = 0.01

#####################################
#Functional Diversities and OMB
#####################################

res_fd <- sapply(FD, function(y) {
  test <- cor.test(y, omb_dat$OMB, method="pearson")
})

# extract coef. of correlation
cor_fd <- res_fd[4, ]
value_fd <- res_fd[3, ]
p.adjusted2 <- p.adjust(value_fd, method= "fdr")

# correlation with individual traits?
m5 <- cor.test(modalities$t10.3, omb_dat$OMB)
m6 <- cor.test(modalities$t11.3, omb_dat$OMB)
p.adjust(c(m5$p.value, m6$p.value), method= "fdr")
# no correlation

# using absolute abundances rather than relative abundances of shredders
# computation of total trait abundance per site
trait_abun_matx <- data.frame(as.matrix(abun) %*% as.matrix(prop_trait))

# check summed abundance for specific modalities 
m7 <- cor.test(trait_abun_matx[ ,"t10.3"], omb_dat$OMB)
m8 <- cor.test(trait_abun_matx[ ,"t11.3"], omb_dat$OMB)

p.adjust(c(m7$p.value, m8$p.value), method= "fdr")
# both statistically significant

######################################
# Code for figures in article
#####################################

#Figure 1
#jpeg('Figure1.jpg', width = 500, height = 150,units = 'mm',res = 600)
par(mfrow=c(1,3), mai=c(0.8, 0.3, 0.5, 3.6), cex= 1.5)
#a)gradient
# account for negative relationship with RDA axes - multiply with -1 to have matching directions
load_spca_n <- load_spca*(-1)
linestack(load_spca_n[abs(load_spca_n) > 0], rownames(load_spca_n)[abs(load_spca_n) > 0], axis = TRUE, cex =1.1)

#b)invertebrates taxa
linestack(sc[abs(sc) > 0.03], rownames(sc)[abs(sc) > 0.03], axis = TRUE, cex =1.1)#rownames includes species names
#c)trait modalities
linestack(sc1[abs(sc1) > .041], trait_labels[abs(sc1) > .041], axis = TRUE, cex =1.1)
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
text(2.4, 15.8, label=expression(paste("r = -0.55; p = 0.004")), cex=0.85)
mtext(text = expression(bold(a)), side = 2, cex=2, las = 1, at = 17, line = 4)
abline(lm(data$TTR~grad), lwd = 1.7)

# Fig 2b
plot(data$SD,data$OMB * 1000,yaxt="n", pch = 16, ylab= expression(paste(italic("k"), " invertebrates (",  10^-3, dday^-1,")")), xlab = "Invertebrate Simpson Diversity")
axis(2,las=2)
cor.test(data$SD,data$OMB,method = c("pearson"))
text(4.5, 5.5, label=expression(paste("r = -0.41, p = 0.06")), cex=0.85)  
mtext(text = expression(bold(b)), side = 2, cex=2, las =1, at = 6.3, line = 4)
abline(lm((data$OMB*1000) ~ data$SD), lwd = 1.7)

# Fig.2c
plot(grad, data$rel.Gam.Abu, yaxt="n", pch = 16, ylab= "Relative Gammarid abundance", xlab ="Environmental stress gradient")
axis(2,las=2)
# cor.test(grad, data$rel.Gam.Abu, method="pearson")
mtext(text = expression(bold(c)), side = 2, cex=2, las =1, at = 110, line = 4)

#Fig 2d
plot(data$Gam.Ab, data$OMB*1000, yaxt="n", pch = 16, ylab= expression(paste(italic("k"), " invertebrates (",  10^-3, dday^-1,")")), xlab = "Number of Gammarids")
axis(2,las=2)
cor.test(data$Gam.Ab, data$OMB*1000, method = c("pearson"))
text(200, 5.5, label=expression(paste("r = 0.48, p = 0.01")), cex = 0.85)  
mtext(text = expression(bold(d)), side = 2, cex=2, las =1, at = 6.4, line = 4)
abline(lm(data$OMB * 1000 ~ data$Gam.Ab), lwd = 1.7)

dev.off()

############################
#Supplementary information
############################

#jpeg("Fig. S1.jpg",width = 230, height = 230,units = 'mm',res = 600)

par(mar=c(4,5,2,2), cex = 1.5)

#Fig. A1
plot(data$Gam.Ab, data$SD, yaxt="n", pch = 16, ylab= "Invertebrate Simpson Diversity", xlab ="Number of Gammarids")
axis(2, las=2)
cor.test(data$SD, data$Gam.Ab, method = c("pearson"))
text(200,5.2, label=expression(paste("r = -0.37, p = 0.045")), cex=1)  
abline(lm(data$SD ~ data$Gam.Ab), lwd = 1.7)

#dev.off()
