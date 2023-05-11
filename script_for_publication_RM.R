#########################################################################
# Height-diameter allometry for a dominant palm in Puerto Rico to improve understanding carbon and forest dynamics

#########################################################################
# Set working directory on the file the scrit is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#########################################################################
# Organize working directory folders
dir.create('figures')
dir.create('tables')

#########################################################################
# Load libraries and read data
library(EDIutils)
library(BIOMASS)
library(viridis)

data <- read.csv("Data_for_analysis.csv")

#########################################################################
# Q1: **Does P. acuminata exhibit a height-diameter relationship and, if so, what is the best model to fit this relationship?**
#########################################################################

# Initialize lists to hold results of cross-validated and full-data models
mod_list <- vector(mode="list", length=7)
names(mod_list) <- c("linear",
                     "loglinear",
                     "loglog",
                     "powerlaw",
                     "loglogquad-1",
                     "loglogquad-2",
                     "weibull")

# Duplicate for full data models
fullmod_list <- mod_list 

########## Model fitting ##########
set.seed(30)
for(i in 1:100) {
  # Partition data for cross validation
  test <- sample(1:nrow(data),
                 size = ceiling (nrow(data) * 0.632),
                 replace = FALSE)
  fit <- data[rownames(data) %in% test,]
  val <- data[!rownames(data) %in% test,]
  
  # Fit models on partitioned data
  mod_list[[1]][[i]] <- lm(height ~ dbh, data = fit)
  mod_list[[2]][[i]] <- lm(height ~ log(dbh), data = fit)
  mod_list[[3]][[i]] <- lm(log(height) ~ log(dbh), data = fit)
  mod_list[[4]][[i]] <- nls(height ~ a * (dbh ^ b), data = fit, start=list(a=1, b=1))
  mod_list[[5]][[i]] <- lm(log(height) ~ I(log(dbh)^2), data=fit)
  mod_list[[6]][[i]] <- lm(log(height) ~ log(dbh) + I(log(dbh)^2), data=fit)
  mod_list[[7]][[i]] <- modelHD(D=fit$dbh,
                                H=fit$height,
                                method="weibull",
                                useWeight=TRUE)$model
  
  # Evaluate models on partitioned data
  for(j in 1:6){
    val_pred <- predict(mod_list[[j]][[i]], newdata=data.frame(dbh=val$dbh))
    mod_list[[j]][[i]]$rmse <- sqrt(sum((val_pred - val$height)^2) / nrow(val))
  }
  val_pred <- predict(mod_list[[7]][[i]], newdata=data.frame(D=val$dbh))
  mod_list[[7]][[i]]$rmse <- sqrt(sum((val_pred - val$height)^2) / nrow(val))
}

# Fit models to full data
fullmod_list[[1]] <- lm(height ~ dbh, data = data)
fullmod_list[[2]] <- lm(height ~ log(dbh), data = data)
fullmod_list[[3]] <- lm(log(height) ~ log(dbh), data = data)
fullmod_list[[4]] <- nls(height ~ a * (dbh ^ b), data = data, start=list(a=1, b=1))
fullmod_list[[5]] <- lm(log(height) ~ I(log(dbh)^2), data=data)
fullmod_list[[6]] <- lm(log(height) ~ log(dbh) + I(log(dbh)^2), data=data)
fullmod_list[[7]] <- modelHD(D=data$dbh, H=data$height, 
                             method="weibull", useWeight=TRUE)$model

########## Compute AIC (full data) ##########
AIC_table <- data.frame(AIC=round(do.call(rbind, lapply(fullmod_list, AIC)), 2))
# Correction for AIC on models with log-transformed y variable
AIC_table$AIC[c(3,5,6)] <- AIC_table$AIC[c(3,5,6)] + 2*sum(log(data$height))
AIC_table$deltaAIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table <- round(AIC_table, 1)

########## Coefficients (full data) ##########
coeffs <- lapply(fullmod_list, function(x) as.vector(coefficients(x))[seq(1:3)])
coeffs <- as.data.frame(round(do.call(rbind, coeffs), 3))
names(coeffs) <- letters[1:3]

########## R-squared (full data) ##########
r2_list <- lapply(fullmod_list, function(x) summary(x)$r.squared)
r2_list[which(unlist(lapply(r2_list, is.null)))] <- NA
r2 <- round(data.frame(r2=do.call(rbind, r2_list)), 2)

########## RSE (full data) ##########
rse <- do.call(rbind, lapply(fullmod_list, function(x) summary(x)$sigma))
rse <- round(data.frame(rse=rse), 2)

########## RMSE (partitioned data) ##########
rmse_list <- do.call(cbind, lapply(mod_list, function(x) do.call(c, lapply(x, function(y) y$rmse))))
rmse <- round(data.frame(rmse_median = apply(rmse_list, 2, median),
                   rmse_sd = apply(rmse_list, 2, sd)), 2)

########## P-values (full data) ##########
pvals <- lapply(fullmod_list, function(x) summary(x)$coefficients[,ncol(summary(x)$coefficients)][seq(1:3)])
pvals <- as.data.frame(round(do.call(rbind, pvals), 2))
names(pvals) <- paste0("p-val_", letters[1:3])

########## Compile to Goodness of fit table ##########
goodness_fit <- cbind(AIC_table, r2, rse, rmse, coeffs, pvals)
goodness_fit

write.csv(goodness_fit, file="tables/Table2.csv")

################################################
## Figure 1: Compare fits of different H:D models

pdf("figures/Figure 1.pdf")

cols <- c('#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F')

plot(data$dbh, data$height,
     main="", xlab="DBH (cm)", 
     ylab="Stem height (m)", pch=16, 
     col=rgb(0,0,0,0.15), ylim=c(0, 20), bty="L")

nd <- data.frame(dbh=seq(0,50,length.out=1000))

lines(seq(0,50,length.out=1000),
      predict(fullmod_list[[1]], newdata=nd),
      col=cols[1], lwd=2, type="l")

lines(seq(0,50,length.out=1000),
      predict(fullmod_list[[2]], newdata=nd),
      col=cols[2], lwd=2, type="l")

lines(seq(1,50,length.out=1000),
       exp(predict(fullmod_list[[3]], newdata=nd)),
       col=cols[3], lwd=2, type="l")

lines(seq(1,50,length.out=1000),
      predict(fullmod_list[[4]], newdata=nd),
      col=cols[4], lwd=2, type="l")

lines(seq(1,50,length.out=1000),
      exp(predict(fullmod_list[[5]], newdata=nd)),
      col=cols[6], lwd=2, type="l")

lines(seq(1,50,length.out=1000),
      exp(predict(fullmod_list[[6]], newdata=nd)),
      col=cols[7], lwd=2, type="l")

lines(seq(1,50,length.out=1000), 
      predict(fullmod_list[[7]], newdata=data.frame(D=seq(1,50,length.out=1000))),
      lwd=2, col=cols[5])

legend('topleft', legend=c("(1) Linear", 
                           "(2) Log-Linear", 
                           "(3) Log-Log", 
                           "(4) Power Law",
                           "(5) Weibull",
                           "(6) Log-Log Quadratic two-term",
                           "(7) Log-Log Quadratic single-term"), 
       bty='n', lwd=2, cex = 0.9,
       col=cols)

dev.off()

################################################################
# Q2. How does environmental heterogeneity mediate allometry of *P. acuminata*?
################################################################

# Fit ordinary linear regressions of environmental covariates to slenderness ratio
m1_env <- lm(log10(SR) ~ Elev, data = data)
m2_env <- lm(log10(SR) ~ Slope, data = data)
m3_env <- lm(log10(SR) ~ log10(nci), data = data)

# RSE
rse_env <- c(summary(m1_env)$sigma, 
             summary(m2_env)$sigma, 
             summary(m3_env)$sigma)

# R2
r2_env <- c(summary(m1_env)$r.squared, 
            summary(m2_env)$r.squared, 
            summary(m3_env)$r.squared)

# p-values
p_value_intercept <- c(summary(m1_env)$coefficient[1,4],
                       summary(m2_env)$coefficient[1,4],
                       summary(m3_env)$coefficient[1,4])

p_value_slope <- c(summary(m1_env)$coefficient[2,4],
                   summary(m2_env)$coefficient[2,4],
                   summary(m3_env)$coefficient[2,4])

# Final goodness of fit table.
good_env_fit <- data.frame(Model = c("Elevation", "Slope", "NCI"),
                                 rse = rse_env, 
                                 r2 = r2_env, 
                                 int_pval = p_value_intercept, 
                                 slope_pval = p_value_slope)
good_env_fit[,-1] <- round(good_env_fit[,-1], 3)

good_env_fit


################################################
### Figure 2. Environmental effects on slenderness ratio

pdf("figures/Figure 2.pdf", width = 7, height = 2.5)

par(mfrow=c(1,3), mar=c(4.5,4.5,1,1))
plot(data$Elev, data$SR, log='y', pch=16, col=rgb(0,0,0,0.35),
     xlab="Elevation (m)", ylab="Slenderness ratio (H/D)", cex.lab=1.25)
abline(m1_env, lwd=3)
mtext("A", 3, line=-2, adj=0.95, font=2)
plot(data$Slope, data$SR, log='y', pch=16, col=rgb(0,0,0,0.35),
     xlab="Slope (degrees)", ylab="Slenderness ratio (H/D)", cex.lab=1.25)
abline(m2_env, lwd=3)
mtext("B", 3, line=-2, adj=0.95, font=2)
plot(data$nci, data$SR, log='xy', pch=16, col=rgb(0,0,0,0.35),
     xlab="NCI", ylab="Slenderness ratio (H/D)", cex.lab=1.25)
abline(m3_env, lwd=3)
mtext("C", 3, line=-2, adj=0.95, font=2)

dev.off()


################################################################
# Q3: **How do estimates of palm AGB depend on H:D model selection?**
################################################################

# Calculate AGB for our measured palms and compare estimations when using inferred and measured height.
# [Lugo's model](http://www.jstor.org/stable/1942582) is AGB = 4.5 + 7.7 \* Height (Kg).
# [Goodman's model](https://www.sciencedirect.com/science/article/abs/pii/S0378112713006592?via%3Dihub) is applicable for Hstem \> 3 m and 6 \<= D \<= 40 cm.
# Here I apply this model to trees SHORTER of 3m.

# Measured height with the Frangi & Lugo (1985) AGB model
data$agb_FL_measured <- 4.5 + 7.7 * data$height

# Height inferred by model (this study) with the Frangi & Lugo (1985) AGB model
data$agb_FL_inferred <- 4.5 + 7.7 * exp(predict(fullmod_list$`loglogquad-1`))

# Measured DBH with the Goodman et al. (2013) AGB model
data$agb_Goodman <- exp(-3.3488 + 2.7483 * log(data$dbh))

# Divide by 1000 to express in Mg and then by 16 to express it Mg/ha
agb <- data.frame(agb_FL_measured=sum((data$agb_FL_measured/1000)/16, na.rm=T),
                  agb_FL_inferred=sum((data$agb_FL_inferred/1000)/16, na.rm=T),
                  agb_G=sum((data$agb_Goodman/1000)/16, na.rm=T))

# Calculate % difference from base model (agb measured). 
biom_diff <-100+round(rbind(((agb[,1] - agb[,1])/agb[,1])*100,
                        ((agb[,2] - agb[,1])/agb[,1])*100,
                        ((agb[,3] - agb[,1])/agb[,1])*100), digits=1)

################################################
pdf("figures/Figure 3.pdf")
# par(mar = c(4,4,3,7))
plot_agb <- barplot(unlist(agb[1,]), beside=TRUE,
                    names.arg=c("Measured", "Inferred", "Goodman"),
                    col=viridis(3),
                    ylab="Palm Above Ground Biomass (Mg / ha)", 
                    ylim=c(0,3.5), xlab="Alternative AGB models", main="")

axis(1, at=plot_agb, labels=F)
pct1 <- biom_diff
text(plot_agb, agb, labels=paste0(pct1,"%"), pos=3)

dev.off()


################################################################
# **Q4: How does AGB of P. acuminata change during a 30 year period following a major hurricane period?**
################################################################

# Load census data with the [EDUutils](https://docs.ropensci.org/EDIutils/) library

raw1 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", 
                         entityId = "041867f47c9c037a510c082cfa577c78")
census1 <- readr::read_csv(file = raw1)

raw2 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", 
                         entityId = "30e32ee3908061ee549ccc797cd4595f")
census2 <- readr::read_csv(file = raw2)

raw3 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", 
                         entityId = "c4cc8b395ca45deeda6572b8cb621881")
census3 <- readr::read_csv(file = raw3)

raw4 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", 
                         entityId = "2d259e63e4de0c7a9b6cf6d57d32d3b8")
census4 <- readr::read_csv(file = raw4)

raw5 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", 
                         entityId = "41cde6946e1f87efb10582dd273c6a30")
census5 <- readr::read_csv(file = raw5)

raw6 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", 
                         entityId = "325c43057e0dd4e1cd6a13fa5125a76d")
census6 <- readr::read_csv(file = raw6)

lfdp <- rbind(census1, census2, census3, census4, census5, census6)

# Convert mm to cm
lfdp$DBH <- lfdp$DBH * 0.1 

# Calculate basal area
lfdp$basal_area <- (pi*((lfdp$DBH/2)^2)/10000) #m^2

####
# Drop all palms with height of measurement below 1m and >1.6m .
premon_total <- subset(lfdp, Latin=="Prestoea acuminata" & 
                         !is.na(lfdp$DBH) & 
                         lfdp$HOM >= 1 &
                         lfdp$HOM <= 1.6)

# Estimate height using new palm allometric model (this study)
premon_total$height_est <- Hlogquad_single(premon_total$DBH)

# Estimate biomass using the Frangi & Lugo (1985) model
premon_total$agb_lugo <- 4.5 + 7.7 * premon_total$height_est

### Estimate biomass using Family level equation Goodman et al. (2013)
# Family-level models are only valid for individuals with Hstem > 3 m and 6<=DBH<=40 cm all of their models used as response variable Hstem NOT Htotal

premon_total$agb_goodman <- exp(-3.3488 + 2.7483 * log(premon_total$DBH))

lugo <- tapply((premon_total$agb_lugo/1000)/16, 
               premon_total$Census, sum, na.rm=TRUE)

goodman <- tapply((premon_total$agb_goodman/1000)/16, 
                  premon_total$Census, sum, na.rm.=TRUE)

bio_dif <- rbind(NA, ((goodman - lugo) / (lugo)) * 100)

################################################
pdf("figures/Figure 4.pdf")
Y <- barplot(rbind(lugo, goodman), beside = TRUE, 
             names.arg = c("1990", "1994", "2000", "2005", "2011", "2016"),
             col = viridis(2),
             ylab = "Above Ground Biomass (Mg / ha)", 
             ylim = c(0,35), xlab = "Census", main = "")

legend('topleft', legend = c("Frangi and Lugo (1985)", 
                             "Goodman et. al. (2013) "), 
       fill = viridis(2), bty = "n", pt.cex = 2)
axis(1, at = Y, labels = F)
pct <- round(bio_dif, 1)
text(Y[2,], goodman, labels = paste0(pct[2,],"%"), pos = 3)
dev.off()

