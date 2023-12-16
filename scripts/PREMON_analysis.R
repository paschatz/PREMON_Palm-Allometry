#### Height-diameter allometry for a dominant palm in Puerto Rico to improve understanding carbon and forest dynamics####

# Organize working directory folders

if(!dir.exists("figures")){dir.create("figures")}
if(!dir.exists("tables")){dir.create("tables")}

#### Load libraries and read data ####
library(renv)
library(EDIutils)
library(BIOMASS)
library(viridis)
library(FSA)
library(tidyverse)
library(ggpubr)
library(ggExtra)

data <- read.csv('clean_data/Data_for_analysis.csv')

#### Reproducible environment ####

#' This script has been developed within a reproducible environment.
#' We used function renv::snapshot() which took a snapshot of the current state 
#' of the project and saved it in the lockfile. This lockfile can be used to 
#' restore the project to the state it was in when the snapshot was taken.
#' In case you try to reproduce this script and you run into troubles with the 
#' version of R or the packages, you can use the lockfile to restore the
#' project to the state it was in when the snapshot was taken.
#' To do so run the following command:

# renv::restore()

#### Q1: **Does P. acuminata exhibit a height-diameter relationship and, if so, what is the best model to fit this relationship?** ####

# Initialize lists to hold results of cross-validated and full-data models
mod_list <- vector(mode = "list", length = 7)
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
for (i in 1:100) {
  # Partition data for cross validation
  test <- sample(1:nrow(data),
                 size = ceiling (nrow(data) * 0.632),
                 replace = FALSE)
  fit <- data[rownames(data) %in% test, ]
  val <- data[!rownames(data) %in% test, ]

  # Fit models on partitioned data
  mod_list[[1]][[i]] <- lm(height ~ dbh, data = fit)
  mod_list[[2]][[i]] <- lm(height ~ log(dbh), data = fit)
  mod_list[[3]][[i]] <- lm(log(height) ~ log(dbh), data = fit)
  mod_list[[4]][[i]] <- nls(height ~ a * (dbh ^ b), data = fit, start = list(a = 1, b = 1))
  mod_list[[5]][[i]] <- lm(log(height) ~ I(log(dbh)^2), data = fit)
  mod_list[[6]][[i]] <- lm(log(height) ~ log(dbh) + I(log(dbh)^2), data = fit)
  mod_list[[7]][[i]] <- modelHD(D = fit$dbh,
                                H = fit$height,
                                method = "weibull",
                                useWeight = TRUE)$model

  # Evaluate models on partitioned data
  # Compute correction factor for log-transformed prediction
  for (j in 1:6){
    if (j %in% c(3, 5, 6)) {
      cf <- logbtcf(mod_list[[j]][[i]])
    } else {
        cf <- 1
    }
    val_pred <- cf * predict(mod_list[[j]][[i]], newdata = data.frame(dbh = val$dbh))
    mod_list[[j]][[i]]$rmse <- sqrt(sum((val_pred - val$height)^2) / nrow(val))
  }
  val_pred <- predict(mod_list[[7]][[i]], newdata = data.frame(D = val$dbh))
  mod_list[[7]][[i]]$rmse <- sqrt(sum((val_pred - val$height)^2) / nrow(val))
}

# Fit models to full data
fullmod_list[[1]] <- lm(height ~ dbh, data = data)
fullmod_list[[2]] <- lm(height ~ log(dbh), data = data)
fullmod_list[[3]] <- lm(log(height) ~ log(dbh), data = data)
fullmod_list[[4]] <- nls(height ~ a * (dbh ^ b), data = data, start = list(a = 1, b = 1))
fullmod_list[[5]] <- lm(log(height) ~ I(log(dbh)^2), data = data)
fullmod_list[[6]] <- lm(log(height) ~ log(dbh) + I(log(dbh)^2), data = data)
fullmod_list[[7]] <- modelHD(D = data$dbh, H = data$height, method = "weibull", useWeight = TRUE)$model

########## Compute AIC (full data) ##########
AIC_table <- data.frame(AIC=round(do.call(rbind, lapply(fullmod_list, AIC)), 2))
# Correction for AIC on models with log-transformed y variable
AIC_table$AIC[c(3, 5, 6)] <- AIC_table$AIC[c(3, 5, 6)] + 2 * sum(log(data$height))
AIC_table$deltaAIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table <- round(AIC_table, 1)

########## Coefficients (full data) ##########
coeffs <- lapply(fullmod_list, function(x) as.vector(coefficients(x))[seq(1:3)])
coeffs <- as.data.frame(round(do.call(rbind, coeffs), 3))
names(coeffs) <- letters[1:3]

########## R-squared (full data) ##########
r2_list <- lapply(fullmod_list, function(x) summary(x)$r.squared)
r2_list[which(unlist(lapply(r2_list, is.null)))] <- NA
r2 <- round(data.frame(r2 = do.call(rbind, r2_list)), 2)

########## RSE (full data) ##########
rse <- do.call(rbind, lapply(fullmod_list, function(x) summary(x)$sigma))
rse <- round(data.frame(rse = rse), 2)

########## RMSE (partitioned data) ##########
rmse_list <- do.call(cbind, lapply(mod_list, function(x) do.call(c, lapply(x, function(y) y$rmse))))
rmse <- round(data.frame(rmse_median = apply(rmse_list, 2, median),
                   rmse_sd = apply(rmse_list, 2, sd)), 2)

########## P-values (full data) ##########
pvals <- lapply(fullmod_list, function(x) summary(x)$coefficients[, ncol(summary(x)$coefficients)][seq(1:3)])
pvals <- as.data.frame(round(do.call(rbind, pvals), 2))
names(pvals) <- paste0("p-val_", letters[1:3])

########## Compile to Goodness of fit table ##########
goodness_fit <- cbind(AIC_table, r2, rse, rmse, coeffs, pvals)
goodness_fit

write.csv(goodness_fit, file = "tables/Table2.csv")

#### Q1b: **Does P. acuminata exhibit a height-diameter relationship and, if so, what is the best model to fit this relationship?** (basal diameter version) ####


# Initialize lists to hold results of cross-validated and full-data models
bd_mod_list <- vector(mode = "list", length = 7)
names(bd_mod_list) <- c("linear",
                        "loglinear",
                        "loglog",
                        "powerlaw",
                        "loglogquad-1",
                        "loglogquad-2",
                        "weibull")

# Duplicate for full data models
bd_fullmod_list <- bd_mod_list

########## Model fitting ##########
for(i in 1:100) {
  # Partition data for cross validation
  test <- sample(1:nrow(data),
                 size = ceiling(nrow(data) * 0.632),
                 replace = FALSE)
  fit <- data[rownames(data) %in% test, ]
  val <- data[!rownames(data) %in% test, ]
  
  fit <- fit[!is.na(fit$basal_d), ]
  val <- val[!is.na(val$basal_d), ]
  
  # Fit models on partitioned data
  bd_mod_list[[1]][[i]] <- lm(height ~ basal_d, data = fit)
  bd_mod_list[[2]][[i]] <- lm(height ~ log(basal_d), data = fit)
  bd_mod_list[[3]][[i]] <- lm(log(height) ~ log(basal_d), data = fit)
  bd_mod_list[[4]][[i]] <- nls(height ~ a * (basal_d ^ b), data = fit, start = list(a=1, b=2))
  bd_mod_list[[5]][[i]] <- lm(log(height) ~ I(log(basal_d)^2), data = fit)
  bd_mod_list[[6]][[i]] <- lm(log(height) ~ log(basal_d) + I(log(basal_d)^2), data=fit)
  bd_mod_list[[7]][[i]] <- modelHD(D = fit$basal_d,
                                   H = fit$height,
                                   method = "weibull",
                                   useWeight = TRUE)$model
  
# Evaluate models on partitioned data
# Compute correction factor for log-transformed prediction
  for(j in 1:6){
    if(j %in% c(3, 5, 6)){
      cf <- logbtcf(bd_mod_list[[j]][[i]])
    } else {
      cf <- 1
    }
    val_pred <- cf * predict(bd_mod_list[[j]][[i]], newdata = data.frame(basal_d = val$basal_d))
    bd_mod_list[[j]][[i]]$rmse <- sqrt(sum((val_pred - val$height)^2) / nrow(val))
  }
  val_pred <- predict(bd_mod_list[[7]][[i]], newdata = data.frame(D = val$basal_d))
  bd_mod_list[[7]][[i]]$rmse <- sqrt(sum((val_pred - val$height)^2) / nrow(val))
}

# Fit models to full data
bd_fullmod_list[[1]] <- lm(height ~ basal_d, data = data)
bd_fullmod_list[[2]] <- lm(height ~ log(basal_d), data = data)
bd_fullmod_list[[3]] <- lm(log(height) ~ log(basal_d), data = data)
bd_fullmod_list[[4]] <- nls(height ~ a * (basal_d ^ b), data = data, start = list(a=1, b=2))
bd_fullmod_list[[5]] <- lm(log(height) ~ I(log(basal_d)^2), data=data)
bd_fullmod_list[[6]] <- lm(log(height) ~ log(basal_d) + I(log(basal_d)^2), data = data)
bd_fullmod_list[[7]] <- modelHD(D = data$basal_d, H = data$height, 
                                method = "weibull", useWeight = TRUE)$model


########## Compute AIC (full data) ##########
AIC_table <- data.frame(AIC = round(do.call(rbind, lapply(bd_fullmod_list, AIC)), 2))
# Correction for AIC on models with log-transformed y variable
AIC_table$AIC[c(3,5,6)] <- AIC_table$AIC[c(3,5,6)] + 2 * sum(log(data$height))
AIC_table$deltaAIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table <- round(AIC_table, 1)

########## Coefficients (full data) ##########
coeffs <- lapply(bd_fullmod_list, function(x) as.vector(coefficients(x))[seq(1:3)])
coeffs <- as.data.frame(round(do.call(rbind, coeffs), 3))
names(coeffs) <- letters[1:3]

########## R-squared (full data) ##########
r2_list <- lapply(bd_fullmod_list, function(x) summary(x)$r.squared)
r2_list[which(unlist(lapply(r2_list, is.null)))] <- NA
r2 <- round(data.frame(r2 = do.call(rbind, r2_list)), 2)

########## RSE (full data) ##########
rse <- do.call(rbind, lapply(bd_fullmod_list, function(x) summary(x)$sigma))
rse <- round(data.frame(rse = rse), 2)

########## RMSE (partitioned data) ##########
rmse_list <- do.call(cbind, lapply(bd_mod_list, function(x) do.call(c, lapply(x, function(y) y$rmse))))
rmse <- round(data.frame(rmse_median = apply(rmse_list, 2, median),
                         rmse_sd = apply(rmse_list, 2, sd)), 2)

########## P-values (full data) ##########
pvals <- lapply(bd_fullmod_list, function(x) summary(x)$coefficients[,ncol(summary(x)$coefficients)][seq(1:3)])
pvals <- as.data.frame(round(do.call(rbind, pvals), 2))
names(pvals) <- paste0("p-val_", letters[1:3])

########## Compile to Goodness of fit table ##########
bd_goodness_fit <- cbind(AIC_table, r2, rse, rmse, coeffs, pvals)
bd_goodness_fit

write.csv(bd_goodness_fit, file = "tables/Table2-basaldiam.csv")

#### Figure 1: Compare fits of different H:D models ####

pdf("figures/Figure_1.pdf", width = 12, height = 6)

par(mfrow = c(1,2))
cols <- c('#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#A2142F', '#4DBEEE', '#77AC30')

# Panel A (DBH)
plot(data$dbh, data$height,
     main = "", xlab = "DBH (cm)",
     ylab = "Stem height (m)", pch = 16,
     col = rgb(0, 0, 0, 0.15), ylim = c(0, 20), bty = "L")

mtext("A", line = -1, adj = 0.95, font = 2)
nd <- data.frame(dbh = seq(0, 50, length.out = 1000))

cf <- logbtcf(fullmod_list[[5]])

lines(seq(1, 50, length.out = 1000),
      cf * exp(predict(fullmod_list[[5]], newdata = nd)),
      col = cols[5], lwd = 5, type = "l")

lines(seq(0, 50, length.out = 1000),
      predict(fullmod_list[[1]], newdata = nd),
      col = cols[1], lwd = 1.5, type = "l")

lines(seq(0,50,length.out = 1000),
      predict(fullmod_list[[2]], newdata = nd),
      col = cols[2], lwd = 1.5, type = "l")

cf <- logbtcf(fullmod_list[[3]])

lines(seq(1, 50, length.out = 1000),
      cf * exp(predict(fullmod_list[[3]], newdata = nd)),
       col = cols[3], lwd = 1.5, type = "l")

lines(seq(1, 50, length.out = 1000),
      predict(fullmod_list[[4]], newdata = nd),
      col = cols[4], lwd = 1.5, type = "l")

cf <- logbtcf(fullmod_list[[6]])

lines(seq(1, 50, length.out = 1000),
      cf * exp(predict(fullmod_list[[6]], newdata = nd)),
      col = cols[6], lwd = 1.5, type = "l")

lines(seq(1, 50, length.out = 1000),
      predict(fullmod_list[[7]], newdata = data.frame(D = seq(1, 50, length.out = 1000))),
      lwd = 1.5, col = cols[7])

legend('topleft', legend = c("(1) Linear",
                             "(2) Log-Linear",
                             "(3) Log-Log",
                             "(4) Power Law",
                             "(5) Log-Log Quadratic I (single-term)",
                             "(6) Log-Log Quadratic II (double-term)",
                             "(7) Weibull"),
       bty = 'n', lwd = c(rep(1.5, 4), 4, rep(1.5, 2)), cex = 0.9,
       col = cols, text.font = c(rep(1, 4), 2, 1, 1))

# Panel B (basal diameter)
plot(data$basal_d, data$height,
     main = "", xlab = "Basal diameter (cm)", 
     ylab = "Stem height (m)", pch = 16, 
     col = rgb(0, 0, 0, 0.15), ylim = c(0, 20), bty = "L")
mtext("B", line = -1, adj = 0.95, font = 2)

nd <- data.frame(basal_d = seq(0, 50, length.out = 1000))

cf <- logbtcf(bd_fullmod_list[[6]])
lines(seq(1, 50, length.out = 1000),
      cf * exp(predict(bd_fullmod_list[[6]], newdata = nd)),
      col = cols[6], lwd = 5, type = "l")

lines(seq(0, 50, length.out = 1000),
      predict(bd_fullmod_list[[1]], newdata = nd),
      col = cols[1], lwd = 1.5, type = "l")

lines(seq(0, 50, length.out = 1000),
      predict(bd_fullmod_list[[2]], newdata = nd),
      col = cols[2], lwd = 1.5, type = "l")

cf <- logbtcf(bd_fullmod_list[[3]])
lines(seq(1, 50, length.out = 1000),
      cf * exp(predict(bd_fullmod_list[[3]], newdata = nd)),
      col = cols[3], lwd = 1.5, type = "l")

lines(seq(1, 50, length.out = 1000),
      predict(bd_fullmod_list[[4]], newdata = nd),
      col = cols[4], lwd = 1.5, type = "l")

cf <- logbtcf(bd_fullmod_list[[5]])
lines(seq(1, 50, length.out = 1000),
      cf * exp(predict(bd_fullmod_list[[5]], newdata = nd)),
      col = cols[5], lwd = 1.5, type = "l")

lines(seq(1, 50, length.out = 1000), 
      predict(bd_fullmod_list[[7]], newdata = data.frame(D=seq(1, 50, length.out = 1000))),
      lwd = 1.5, col = cols[7])

legend('topleft', legend = c("(1b) Linear",
                             "(2b) Log-Linear",
                             "(3b) Log-Log",
                             "(4b) Power Law",
                             "(5b) Log-Log Quadratic I (single-term)",
                             "(6b) Log-Log Quadratic II (double-term)",
                             "(7b) Weibull"),
       bty = 'n', lwd = c(rep(1.5, 5), 4, rep(1.5, 1)), cex = 0.9,
       col = cols, text.font = c(rep(1, 5), 2, 1))

dev.off()

#### Q2. How does environmental heterogeneity mediate allometry of *P. acuminata*? ####

data$lognci <- log10(data$nci)

# Fit multiple linear regression of environmental covariates to slenderness ratio
m0 <- lm(SR ~ Elev + Slope + lognci, data = data)
summary(m0)

data$Elev.z <- scale(data$Elev)
data$Slope.z <- scale(data$Slope)
data$lognci.z <- scale(data$lognci)

m0.z <- lm(SR ~ Elev.z + Slope.z + lognci.z, data = data)
summary(m0.z)


### STATS
# RSE
rse_env <- round(summary(m0.z)$sigma, 3)

# R2
r2_env <- round(summary(m0.z)$r.squared, 3)

# p-values
p_value_intercept <- summary(m0.z)$coefficient[1, 4]

p_value_slope <- c(summary(m0.z)$coefficient["Elev.z", 4],
                   summary(m0.z)$coefficient["Slope.z", 4],
                   summary(m0.z)$coefficient["lognci.z", 4])

# Final goodness of fit table.
good_env_fit <- data.frame(rse = rse_env,
                           r2 = r2_env,
                           int_pval = round(p_value_intercept, 3),
                           elev_pval = round(p_value_slope[1], 3),
                           slope_pval = round(p_value_slope[2], 3),
                           nci_pval = round(p_value_slope[3], 3))

good_env_fit

#### Figure 2. Environmental effects on slenderness ratio ####

pdf("figures/Figure_2.pdf", width = 7, height = 2.5)

par(mfrow = c(1, 3), mar = c(4.5, 4.5, 1, 1))

# PLOT
par(mfrow=c(1,3), mar=c(4.5,4.5,1,1))
plot(data$Elev, data$SR, log='y', pch=16, col=rgb(0,0,0,0.35),
     xlab="Elevation (m)", ylab="Slenderness ratio (H/D)", cex.lab=1.25)
mtext("A", 3, line=-2, adj=0.95, font=2)

# Predict for Elevation
pred <- data.frame(Elev = seq(min(data$Elev, na.rm=T), max(data$Elev, na.rm=T), length.out=101),
                   Slope = rep(mean(data$Slope, na.rm=T)),
                   lognci = rep(mean(data$lognci, na.rm=T)))
conf_interval <- predict(m0, newdata = pred, interval="confidence", level = 0.95)

# Plot Elevation effect
polygon(x=c(pred$Elev, rev(pred$Elev)),
        y=c(conf_interval[,2], rev(conf_interval[,3])), 
        lty=0, col=rgb(0,0,1,0.5))
lines(pred$Elev, conf_interval[,1], lwd=2, col=1, lty=2)

plot(data$Slope, data$SR, log='y', pch=16, col=rgb(0,0,0,0.35),
     xlab="Slope (degrees)", ylab="Slenderness ratio (H/D)", cex.lab=1.25)
mtext("B", 3, line=-2, adj=0.95, font=2)

# Predict for Slope
pred <- data.frame(Elev = rep(mean(data$Elev, na.rm=T)),
                   Slope = seq(min(data$Slope, na.rm=T), max(data$Slope, na.rm=T), length.out=101),
                   lognci = rep(mean(data$lognci, na.rm=T)))
conf_interval <- predict(m0, newdata = pred, interval="confidence", level = 0.95)

# Plot Slope effect
polygon(x=c(pred$Slope, rev(pred$Slope)),
        y=c(conf_interval[,2], rev(conf_interval[,3])), 
        lty=0, col=rgb(0,0,1,0.5))
lines(pred$Slope, conf_interval[,1], lwd=2, col=1)

plot(data$nci, data$SR, log='xy', pch=16, col=rgb(0,0,0,0.35),
     xlab="NCI", ylab="Slenderness ratio (H/D)", cex.lab=1.25)
mtext("C", 3, line=-2, adj=0.95, font=2)

# Predict for NCI
pred <- data.frame(Elev = rep(mean(data$Elev, na.rm=T)),
                   Slope = rep(mean(data$Slope, na.rm=T)),
                   lognci = log10(seq(min(data$nci, na.rm=T), max(data$nci, na.rm=T), length.out=101)))
conf_interval <- predict(m0, newdata = pred, interval="confidence", level = 0.95)

# Plot NCI effect
polygon(x=10^(c(pred$lognci, rev(pred$lognci))),
        y=c(conf_interval[,2], rev(conf_interval[,3])), 
        lty=0, col=rgb(0,0,1,0.5))
lines(10^pred$lognci, conf_interval[,1], lwd=2, col=1)

dev.off()

#### Q3: **How do estimates of palm AGB depend on H:D model selection?** ####

# Calculate AGB for our measured palms and compare estimations when using inferred and measured height.
# [Frangi and Lugo 1985 model](http://www.jstor.org/stable/1942582)
# [Goodman et al. 2013 model](https://www.sciencedirect.com/science/article/abs/pii/S0378112713006592?via%3Dihub)
# [Avalos et al. 2022 model](https://www.frontiersin.org/articles/10.3389/ffgc.2022.867912/full)

# Frangi & Lugo (1985) AGB model with measured stem height
data$agb_FL_measured <- 4.5 + 7.7 * data$height

# Goodman et al. (2013) family-wise AGB model with measured DBH
data$agb_Goodman <- exp(-3.3488 + 2.7483 * log(data$dbh))

# Avalos et al. (2022) family-wise AGB model with measured DBH
data$agb_Avalos <- 2 * (1.4 * exp(-4.77 + 2.82 * log(data$dbh)))

# Divide by 1000 to express in Mg and then by 16 to express it Mg
agb <- data.frame(agb_FL_measured = sum((data$agb_FL_measured / 1000), na.rm = T),
                  agb_G = sum((data$agb_Goodman / 1000), na.rm = T),
                  agb_A = sum((data$agb_Avalos / 1000), na.rm = T))

biom_diff <- round(rbind(((agb[,1] - agb[,1])/agb[,1])*100,
                         ((agb[,2] - agb[,1])/agb[,1])*100,
                         ((agb[,3] - agb[,1])/agb[,1])*100), digits=1)

#### Figure 3 ####
pdf("figures/Figure_3.pdf")

### Density plots of individual AGB estimates ####
d1 <- density(data$agb_FL_measured, na.rm=T)
plot(d1, col=viridis(3)[1], 
     xlab="Estimated AGB (kg)", lwd=3, main=NA,
     ylim=c(0,0.0325), xlim=c(0,150))
polygon(d1, col=scales::alpha(viridis(3)[1], 0.2), border="NA")

d2 <- density(data$agb_Goodman, na.rm=T)
lines(d2, col=viridis(3)[2], lwd=3)
polygon(d2, col=scales::alpha(viridis(3)[2], 0.2), border="NA")

d3 <- density(data$agb_Avalos, na.rm=T)
lines(d3, col=viridis(3)[3], lwd=3)
polygon(d3, col=scales::alpha(viridis(3)[3], 0.2), border="NA")

legend('right', legend=c("Frangi & Lugo (1985)", 
                            "Goodman et al. (2013)", 
                            "Avalos et al. (2022)"), 
       lwd=3, col=viridis(3), bty='n', pch=22, pt.lwd=1, pt.cex=2, 
       pt.bg=scales::alpha(viridis(3), 0.2))

boxplot(data$agb_FL_measured, data$agb_Goodman, data$agb_Avalos, 
      col = viridis(3),
      horizontal = TRUE,
      add = TRUE,
      at = c(0.0325, 0.03, 0.0275),
      boxwex = 0.002,
      axes = F,
      outcol = viridis(3))

      dev.off()

# What is the median and sd of difference values between the different model estimates?
round(median(data$agb_Goodman - data$agb_FL_measured), 2)
round(sd(data$agb_Goodman - data$agb_FL_measured), 2)

round(median(data$agb_Avalos - data$agb_FL_measured), 2)
round(sd(data$agb_Avalos - data$agb_FL_measured), 2)

# **Q4: How does AGB of P. acuminata change during a 30 year period following a major hurricane period?** ####

# Load census data with the [EDUutils](https://docs.ropensci.org/EDIutils/) library

raw1 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", entityId = "041867f47c9c037a510c082cfa577c78")
census1 <- readr::read_csv(file = raw1)

raw2 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", entityId = "30e32ee3908061ee549ccc797cd4595f")
census2 <- readr::read_csv(file = raw2)

raw3 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", entityId = "c4cc8b395ca45deeda6572b8cb621881")
census3 <- readr::read_csv(file = raw3)

raw4 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", entityId = "2d259e63e4de0c7a9b6cf6d57d32d3b8")
census4 <- readr::read_csv(file = raw4)

raw5 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", entityId = "41cde6946e1f87efb10582dd273c6a30")
census5 <- readr::read_csv(file = raw5)

raw6 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", entityId = "325c43057e0dd4e1cd6a13fa5125a76d")
census6 <- readr::read_csv(file = raw6)

lfdp <- rbind(census1, census2, census3, census4, census5, census6)

# Convert mm to cm
lfdp$DBH <- lfdp$DBH * 0.1 

# Calculate basal area
lfdp$basal_area <- (pi * ((lfdp$DBH / 2) ^ 2) / 10000) #m^2

### Calculate non-palm AGB using the BIOMASS package
# Get wood density values for non-palm trees
wd <- getWoodDensity(genus = unlist(lapply(strsplit(lfdp$Latin, " "), 
                     function(x) x[[1]])),
      species = unlist(lapply(strsplit(lfdp$Latin, " "), 
                     function(x) x[[2]])))

# Compute AGB
# might be necessary to lead library(httr2)
lfdp$AGB <- computeAGB(lfdp$DBH, WD = wd$meanWD, 
                       coord = c(-65.8246, 18.3214))

# Make sure AGB for palms computed in this way is NA
lfdp$AGB[lfdp$Mnemonic == 'PREMON'] <- NA

### Estimate palm biomass
# Drop all palms with height of measurement below 1m and >1.6m .
premon_total <- subset(lfdp, Latin == "Prestoea acuminata" & 
                             !is.na(lfdp$DBH) & 
                             lfdp$HOM >= 1 &
                             lfdp$HOM <= 1.6)

# Estimate height using the best palm H:D model from this study
cf <- logbtcf(fullmod_list$`loglogquad-1`)
premon_total$height_est <- cf * exp(predict(fullmod_list$`loglogquad-1`, 
                                            data.frame(dbh = premon_total$DBH)))

# Estimate the lower and upper bounds of height to quantify a range of uncertainty
premon_total$height_est_min <- cf * exp(predict(fullmod_list$`loglogquad-1`, data.frame(dbh = premon_total$DBH), interval = 'confidence'))[, 2]
premon_total$height_est_max <- cf * exp(predict(fullmod_list$`loglogquad-1`, data.frame(dbh = premon_total$DBH), interval = 'confidence'))[, 3]

# Estimate palm biomass using the Frangi & Lugo (1985) model
premon_total$agb_lugo <- 4.5 + 7.7 * premon_total$height_est
premon_total$agb_lugo_min <- 4.5 + 7.7 * premon_total$height_est_min
premon_total$agb_lugo_max <- 4.5 + 7.7 * premon_total$height_est_max

### Estimate biomass using Family level equation Goodman et al. (2013)
# Family-level models are only valid for individuals with Hstem > 3 m and 6<=DBH<=40 cm all of their models used as response variable Hstem NOT Htotal
premon_total$agb_goodman <- exp(-3.3488 + 2.7483 * log(premon_total$DBH))

# Measured DBH with the Avalos et al. (2022) family-wise AGB model
premon_total$agb_avalos <- 2 * (1.4 * exp(-4.77 + 2.82 * log(premon_total$DBH)))

lugo <- tapply((premon_total$agb_lugo / 1000) / 16, 
                premon_total$Census, sum, na.rm = TRUE)

lugo_min <- tapply((premon_total$agb_lugo_min / 1000) / 16, 
                    premon_total$Census, sum, na.rm = TRUE)

lugo_max <- tapply((premon_total$agb_lugo_max / 1000) / 16, 
                    premon_total$Census, sum, na.rm = TRUE)

goodman <- tapply((premon_total$agb_goodman / 1000) / 16,
                   premon_total$Census, sum, na.rm = TRUE)

avalos <- tapply((premon_total$agb_avalos / 1000) / 16, 
                  premon_total$Census, sum, na.rm = TRUE)

bio_dif <- rbind(NA, ((goodman - lugo) / (lugo)) * 100,
                     ((avalos - lugo) / (lugo)) * 100)

round(100 * (lugo_max - lugo_min) / lugo, 2)

# Figure 4 ####
pdf("figures/Figure_4.pdf")

Y <- barplot(rbind(lugo, goodman, avalos), beside = TRUE, 
             names.arg = c("1990", "1994", "2000", "2005", "2011", "2016"),
             col = viridis(3),
             ylab = "Palm Above Ground Biomass (Mg / ha) in the LFDP", 
             ylim = c(0, 35), xlab = "Census", main = "")

arrows(Y[1, ], lugo_min, Y[1, ], lugo_max, angle = 90, code = 3, len = 0.05)
arrows(Y[1, ], lugo, Y[1, ], lugo_min, angle = 90, code = 2, len = 0.05, col = 'white')

legend('topleft', legend = c("Frangi and Lugo (1985)", 
                             "Goodman et al. (2013)",
                             "Avalos et al. (2022)"), 
                  fill = viridis(3), 
                  bty = "n", pt.cex = 2)
                  axis(1, at = Y[2, ], labels = F)

pct <- round(bio_dif, 1)

text(Y[2, ], goodman - 0.5, labels = paste0(pct[2, ], "%"), 
     pos = 3, cex = 0.65)

text(Y[3, ] + 0.25, avalos - 0.5, labels = paste0(pct[3, ], "%"), 
     pos = 3, cex = 0.65)

dev.off()

# Descriptive Figure 5 with (A) stem density, (B) basal area, and (C) AGB ####

pdf("figures/Figure_5.pdf", width = 7, height = 5)

Year <- c(1990, 1994, 2000, 2005, 2011, 2016)

par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))

### Absolute value plots
# Stem density
plot(Year, as.vector(table(lfdp$Census[!is.na(lfdp$DBH)])) / 16, type = 'b', log = '', 
                   ylim = c(100, 8000), axes = F, pch = 16,
                   ylab = "Stems per ha", )

points(Year, as.vector(table(premon_total$Census[!is.na(lfdp$DBH)])) / 16,
             type = 'b', col = 'blue', pch = 16)

legend('topright', legend = c("All LFDP stems", "P. acuminata"), 
             text.font = c(1, 3),
                         col = c(1, 'blue'), pch = 1, lty = 1, bty = 'n')
axis(1); axis(2, las = 2)
mtext("A", line = -1, adj = 0.05, font = 2)

# Basal area
plot(Year, tapply(lfdp$basal_area[!is.na(lfdp$DBH)], lfdp$Census[!is.na(lfdp$DBH)], 
           sum, na.rm = T) / 16, type = 'b', 
    ylim = c(0, 50),
    axes = F, pch = 16,
    ylab = "Basal area (m^2) per ha", )

points(Year, tapply(premon_total$basal_area[!is.na(premon_total$DBH)], 
                    premon_total$Census[!is.na(premon_total$DBH)], 
                    sum, na.rm = T) / 16, type = 'b', col = 'blue', pch = 16)

axis(1); axis(2, las = 2)
mtext("B", line = -1, adj = 0.05, font = 2)

# AGB
plot(Year, tapply(lfdp$AGB[!is.na(lfdp$DBH)], 
                  lfdp$Census[!is.na(lfdp$DBH)], 
                  sum, na.rm = T) / 16, type = 'b', 
           ylim = c(0, 350),
           axes = F,
           ylab = "Estimated AGB (Mg per ha)", pch = 16)

points(Year, (tapply(premon_total$agb_lugo[!is.na(premon_total$DBH)], 
                     premon_total$Census[!is.na(premon_total$DBH)], 
      sum, na.rm = T) / 16) / 1000, type = 'b', col = 'blue', pch = 16)
      axis(1); axis(2, las = 2)

mtext("C", line = -1, adj = 0.05, font = 2)

### Proportional plots
# Stem density

plot(Year, (as.vector(table(premon_total$Census[!is.na(lfdp$DBH)]))) /
           (as.vector(table(lfdp$Census[!is.na(lfdp$DBH)]))), type = 'b', 
           ylim = c(0, 1), 
           axes = F,
           ylab = "Palm proportion of stems", col = 'red', pch = 16)

axis(1); axis(2, las = 2)
mtext("D", line = -1, adj = 0.05, font = 2)

# Basal area
plot(Year, (tapply(premon_total$basal_area[!is.na(premon_total$DBH)], 
                   premon_total$Census[!is.na(premon_total$DBH)], 
            sum, na.rm = T) / 16) / (tapply(lfdp$basal_area[!is.na(lfdp$DBH)], 
            lfdp$Census[!is.na(lfdp$DBH)], 
            sum, na.rm = T) / 16), type = 'b', 
            ylim = c(0, 1), axes = F,
            ylab = "Palm proportion of basal area", col = 'red', pch = 16)

axis(1); axis(2, las = 2)
mtext("E", line = -1, adj = 0.05, font = 2)

# AGB
plot(Year, ((tapply(premon_total$agb_lugo[!is.na(premon_total$DBH)], 
                    premon_total$Census[!is.na(premon_total$DBH)], 
            sum, na.rm = T) / 16) / 1000) / (tapply(lfdp$AGB[!is.na(lfdp$DBH)], 
            lfdp$Census[!is.na(lfdp$DBH)], 
            sum, na.rm = T) / 16), type = 'b', 
            ylim = c(0, 1), axes = F,
      ylab = "Palm proportion of total estimated AGB", col = 'red', pch = 16)

axis(1); axis(2, las = 2)
mtext("F", line = -1, adj = 0.05, font = 2)

dev.off()

# Supplemental Figures ####
# Figure S1 ####

pdf("figures/Figure_S1.pdf")

# Raw height:DBH plot
raw_dbh <- ggplot(data, aes(x = dbh, y = height)) +
            geom_point(size = 2, shape = 20, color = "gray45") +
            xlab("") +
            ylab("Stem Height (m)") + 
            theme_minimal()

# DBH density plot
den_dbh <- ggplot(data, aes(x = dbh)) + 
            geom_histogram(aes(y = ..density..),
            binwidth = 0.5, #breaks
            colour = "black", fill = "white") +
            geom_density(alpha = 0.5, #density transparency
            color = "black", fill = "blue") +
            xlab("DBH (cm)") + 
            ylab("Density") + 
            theme_minimal()

# Raw height:basal_d plot
raw_basal <- ggplot(data, aes(x = basal_d, y = height)) +
            geom_point(size = 2, shape = 20, color = "gray45") + 
            xlab("") + 
            ylab("") + 
            theme_minimal()

# Basal_diameter density plot
den_basal <- ggplot(data, aes(x = basal_d, na.rm = TRUE)) + 
            geom_histogram(aes(y = ..density..), binwidth = 1, #breaks
            colour = "black", fill = "white") +
            geom_density(alpha = 0.5, #density transparency
            color = "black", fill = "red1") +
            xlab("Basal diameter (cm)") + 
            ylab("") + 
            theme_minimal()

t1 <- ggMarginal(raw_basal + theme_gray() + xlab("Basal Diameter (cm)") +
                 ylab("Stem height (m)") +
                 theme_minimal(), fill = "lightblue")

t2 <- ggMarginal(raw_dbh + theme_gray() + xlab("DBH (cm)") + 
                 ylab("Stem Height") + 
                 theme_minimal(), fill = "tomato")

ggarrange(t1, t2, labels = c("A", "B"), ncol = 1, nrow = 2)

dev.off()


# Figure S2 ####

pdf("figures/Figure_S2.pdf")

hist(census6$DBH[census6$Mnemonic == 'PREMON'] / 10, main = NA, 
                   xlab = "Diameter at breast height (D130) (cm)",
                   ylab = "Number of palms")

hist(data$dbh, add = TRUE, col = 2, breaks = 20)

legend('topleft', legend = c("LFDP 2016 Census", "January 2020 sampling"), 
       pch = 22, pt.bg = c('grey', 2), pt.cex = 2, inset = 0.05)

dev.off()

#### Graphical abstract ####

# Panel A (DBH)
pdf("figures/Graphical_abstract.pdf", width = 8, height = 6)  

par(mar = c(6, 4, 4, 2) + 0.1)  

plot(data$dbh, data$height,
     main = "", xlab = "Diameter at Breast Height (DBH) - cm",
     ylab = "Stem Height - m", pch = 16,
     col = rgb(0, 0, 0, 0.15), ylim = c(0, 20), bty = "L")

nd <- data.frame(dbh = seq(0, 50, length.out = 1000))

title(main = expression('Advancing understanding of tropical forest carbon dynamics through improved H:D allometric models for '*italic(P.~accuminata)*''),
      cex.main = 1.5, line = 2) 

cf <- logbtcf(fullmod_list[[5]])
lines(seq(1, 50, length.out = 1000),
      cf * exp(predict(fullmod_list[[5]], newdata = nd)),
      col = cols[5], lwd = 5, type = "l")

lines(seq(0, 50, length.out = 1000),
      predict(fullmod_list[[1]], newdata = nd),
      col = cols[1], lwd = 1.5, type = "l")

lines(seq(0,50,length.out = 1000),
      predict(fullmod_list[[2]], newdata = nd),
      col = cols[2], lwd = 1.5, type = "l")

cf <- logbtcf(fullmod_list[[3]])
lines(seq(1, 50, length.out = 1000),
      cf * exp(predict(fullmod_list[[3]], newdata = nd)),
      col = cols[3], lwd = 1.5, type = "l")

lines(seq(1, 50, length.out = 1000),
      predict(fullmod_list[[4]], newdata = nd),
      col = cols[4], lwd = 1.5, type = "l")

cf <- logbtcf(fullmod_list[[6]])
lines(seq(1, 50, length.out = 1000),
      cf * exp(predict(fullmod_list[[6]], newdata = nd)),
      col = cols[6], lwd = 1.5, type = "l")

lines(seq(1, 50, length.out = 1000),
      predict(fullmod_list[[7]], newdata = data.frame(D = seq(1, 50, length.out = 1000))),
      lwd = 1.5, col = cols[7])
dev.off()
