library(tidyverse)
library(haven)
library(mirt)
source("mDR_IRT.R")
Rcpp::sourceCpp("drm.cpp")

################################################################################
# KIMS
################################################################################

# read data
data <- read.csv("KIMS/data.csv") %>%
  select(starts_with("Q"))

# data cleaning
unique(c(as.matrix(data)))
data[data==0] <- NA
data <- data[apply(data, 1, function(x)mean(is.na(x)))!=1,]
data <- data[apply(data, 1, function(x)sd(x, na.rm = T))!=0,]

# frequency table
apply(data, 2, table)

# GPCM
fit.gpcm <- mirt::mirt(data, 1 ,quadpts=41,itemtype = "gpcm")

# ------------------------------------------------------------------------------
# complex model

# formula string to run CFA
fs <- paste("f1", paste(colnames(data), collapse = " + "), sep = " ~ ")

# standard deviation of the prior
t_sd <- t_prior_sd(
  N = nrow(data),
  category = 5,
  n = 1
)

# fit a CFA model
fit <- cfa_sim(fs, data, t_prior = t_sd)

# parameter estimates
fit$par_est

# plotting
thresholds <- apply(fit$par_est[[2]][,-1], 1, cut_trans, simplify = FALSE)
plts <- list()

for(i in 1:length(thresholds)){
  plts[[i]] <- plot_DRM(thresholds[[i]], paste0("item",sprintf("%02d", i)))
}
independent_intervals <- do.call(
  gridExtra::arrangeGrob,
  append(
    plts,  # list of grobs
    list(
      ncol = 1,
      top = grid::textGrob(
        "(a) Complex Model",
        x = 0.01,
        hjust = 0,
        gp = grid::gpar(fontsize = 13, fontface = "bold")
      ),
      padding = unit(1, "lines")
    )
  )
)
grid::grid.draw(independent_intervals)
# ------------------------------------------------------------------------------
# parallel model

# formula string to run CFA
fs2 <- paste0(
  fs,"\n",
  paste(paste0(colnames(data), ".t1"), collapse = "=="),"\n",
  paste(paste0(colnames(data), ".t2"), collapse = "=="),"\n",
  paste(paste0(colnames(data), ".t3"), collapse = "=="),"\n",
  paste(paste0(colnames(data), ".t4"), collapse = "=="),
  collapse = "")

cat(fs2)

# standard deviation of the prior
t_sd2 <- t_prior_sd(
  N = nrow(data),
  category = 5,
  n = 39
)
t_2 <- rep(t_sd2, 39)

# fit a CFA model
fit2 <- cfa_sim(fs2, data, t_prior = t_2)

# parameter estimates
fit2$par_est

# plotting
thresholds <- apply(fit2$par_est[[2]][,-1], 1, cut_trans, simplify = FALSE)
plts <- list()

for(i in 1:length(thresholds)){
  plts[[i]] <- plot_DRM(thresholds[[i]], paste0("item",sprintf("%02d", i)))
}
parallel_intervals <- do.call(
  gridExtra::arrangeGrob,
  append(
    plts,  # list of grobs
    list(
      ncol = 1,
      top = grid::textGrob(
        "(b) Parallel Model",
        x = 0.01,
        hjust = 0,
        gp = grid::gpar(fontsize = 13, fontface = "bold")
      ),
      padding = unit(1, "lines")
    )
  )
)
gridExtra::grid.arrange(independent_intervals, parallel_intervals, ncol=2)
# ------------------------------------------------------------------------------
# equal-interval model

# fit a CFA model
fit3 <- cfa_sim(fs, data, t_prior = t_2, eq_interval = TRUE)

# parameter estimates
fit3$par_est

# ------------------------------------------------------------------------------
# Model Selection

# log-likelihood
fit$logL
fit2$logL
fit3$logL
mirt::logLik(fit.gpcm)

# eff par
fit$eff_par
fit2$eff_par
fit3$eff_par
par_gpcm <- 5 * 39

# AIC
fit$AIC
fit2$AIC
fit3$AIC
-2 * mirt::logLik(fit.gpcm) + 2 * par_gpcm

# BIC
fit$BIC
fit2$BIC
fit3$BIC
-2 * mirt::logLik(fit.gpcm) + log(mean(colSums(!is.na(data)))) * par_gpcm

# RMSEA
rmsea_drm(fit, fit2)
rmsea_drm(fit2, fit3)

################################################################################
# CFCS
################################################################################

# read data
data <- read.csv("CFCS/data.csv", sep = "\t", header = TRUE) %>%
  select(starts_with("Q"))

# data cleaning
unique(c(as.matrix(data)))
data[data == 0] <- NA
data[data == -1] <- NA

data <- data[apply(data, 1, function(x)mean(is.na(x)))!=1,]
data <- data[apply(data, 1, function(x)sd(x, na.rm = T))!=0,]

# frequency table
apply(data, 2, table)

# GPCM
fit.gpcm <- mirt::mirt(data, 1 ,quadpts=41,itemtype = "gpcm")

# ------------------------------------------------------------------------------
# complex model

# formula string to run CFA
fs <- paste("f1", paste(colnames(data), collapse = " + "), sep = " ~ ")

# standard deviation of the prior
t_sd <- t_prior_sd(
  N = nrow(data),
  category = 5,
  n = 1
)

# fit a CFA model
fit <- cfa_sim(fs, data, t_prior = t_sd)

# parameter estimates
fit$par_est

# plotting
thresholds <- apply(fit$par_est[[2]][,-1], 1, cut_trans, simplify = FALSE)
plts <- list()

for(i in 1:length(thresholds)){
  plts[[i]] <- plot_DRM(thresholds[[i]], paste0("item",sprintf("%02d", i)))
}
independent_intervals <- do.call(
  gridExtra::arrangeGrob,
  append(
    plts,  # list of grobs
    list(
      ncol = 1,
      top = grid::textGrob(
        "(a) Complex Model",
        x = 0.01,
        hjust = 0,
        gp = grid::gpar(fontsize = 13, fontface = "bold")
      ),
      padding = unit(1, "lines")
    )
  )
)
# ------------------------------------------------------------------------------
# parallel model

# formula string to run CFA
fs2 <- paste0(
  fs,"\n",
  paste(paste0("Q", 1:12, ".t1"), collapse = "=="),"\n",
  paste(paste0("Q", 1:12, ".t2"), collapse = "=="),"\n",
  paste(paste0("Q", 1:12, ".t3"), collapse = "=="),"\n",
  paste(paste0("Q", 1:12, ".t4"), collapse = "=="),
  collapse = "")

cat(fs2)

# standard deviation of the prior
t_sd2 <- t_prior_sd(
  N = nrow(data),
  category = 5,
  n = 12
)
t_2 <- rep(t_sd2, 12)

# fit a CFA model
fit2 <- cfa_sim(fs2, data, t_prior = t_2)

# parameter estimates
fit2$par_est

# plotting
thresholds <- apply(fit2$par_est[[2]][,-1], 1, cut_trans, simplify = FALSE)
plts <- list()

for(i in 1:length(thresholds)){
  plts[[i]] <- plot_DRM(thresholds[[i]], paste0("item",sprintf("%02d", i)))
}
parallel_intervals <- do.call(
  gridExtra::arrangeGrob,
  append(
    plts,  # list of grobs
    list(
      ncol = 1,
      top = grid::textGrob(
        "(b) Parallel Model",
        x = 0.01,
        hjust = 0,
        gp = grid::gpar(fontsize = 13, fontface = "bold")
      ),
      padding = unit(1, "lines")
    )
  )
)
gridExtra::grid.arrange(independent_intervals, parallel_intervals, ncol=2)

# ------------------------------------------------------------------------------
# Model Selection

# log-likelihood
fit$logL
fit2$logL
mirt::logLik(fit.gpcm)

# eff par
fit$eff_par
fit2$eff_par
par_gpcm <- 5 * 12

# AIC
sprintf("%.2f", fit$AIC)
sprintf("%.2f", fit2$AIC)
-2 * mirt::logLik(fit.gpcm) + 2 * par_gpcm

# BIC
sprintf("%.2f", fit$BIC)
sprintf("%.2f", fit2$BIC)
-2 * mirt::logLik(fit.gpcm) + log(mean(colSums(!is.na(data)))) * par_gpcm

# RMSEA
rmsea_drm(fit, fit2)

################################################################################
# PTGI
################################################################################
data <- read_sav("PTGI/PTGI.sav")
data <- data[,4:45] %>% round(0)
colnames(data) <- paste0("q", 1:42)
dim(data)

# data cleaning
unique(c(as.matrix(data)))
data <- data[apply(data, 1, function(x)mean(is.na(x)))!=1,]

# frequency table
apply(data, 2, table)

# GPCM
model <- "f1 = q1-q21
          f2 = q22-q42
          COV = f1*f2"
fit.gpcm <- mirt::mirt(data,model,quadpts=41,itemtype = "gpcm", TOL = 0.0002)

# ------------------------------------------------------------------------------
# complex model

# formula string to run CFA
formula_string <- paste(
  paste("\n f1", paste(colnames(data)[1:21], collapse = " + "), sep = " ~ "),
  paste("\n f2", paste(colnames(data)[22:42], collapse = " + "), sep = " ~ "),
  collapse = "\n"
)

# standard deviation of the prior
t_sd <- t_prior_sd(
  N = 405,
  category = 6,
  n = 1
)

# fit a CFA model
fit.cfa <- cfa_sim(formula_string, data, t_prior = t_sd, threshold = 0.000002)

# plotting
thresholds1 <- t(apply(fit.cfa$par_est[[2]][1:21,-1], 1, cut_trans))
thresholds2 <- t(apply(fit.cfa$par_est[[2]][22:42,-1], 1, cut_trans))
plts <- list()

for(i in 1:nrow(thresholds1)){
  plts[[i]] <- plot_DRM_compare(thresholds1[i,], thresholds2[i,], label = c("P", "N"),item_name = paste0("item pair: ",sprintf("%02d ", i)))
}
do.call(gridExtra::grid.arrange, c(plts, ncol=1))

# ------------------------------------------------------------------------------
# wording-effect model

# formula string to run CFA
formula_string2 <- paste0(
  formula_string,
  "\n",
  paste(paste0("q", 1:21, ".t1"), collapse = "=="),
  "\n",
  paste(paste0("q", 1:21, ".t2"), collapse = "=="),
  "\n",
  paste(paste0("q", 1:21, ".t3"), collapse = "=="),
  "\n",
  paste(paste0("q", 1:21, ".t4"), collapse = "=="),
  "\n",
  paste(paste0("q", 1:21, ".t5"), collapse = "=="),
  "\n",
  paste(paste0("q", 22:42, ".t1"), collapse = "=="),
  "\n",
  paste(paste0("q", 22:42, ".t2"), collapse = "=="),
  "\n",
  paste(paste0("q", 22:42, ".t3"), collapse = "=="),
  "\n",
  paste(paste0("q", 22:42, ".t4"), collapse = "=="),
  "\n",
  paste(paste0("q", 22:42, ".t5"), collapse = "=="),
  collapse = ""
)

# standard deviation of the prior
t_sd <- t_prior_sd(
  N = 405,
  category = 6,
  n = 21
)

# fit a CFA model
fit.cfa2 <- cfa_sim(formula_string2, data, t_prior = t_sd)

# covariance matrix of the latent variables
fit.cfa2$cov_mat


# ------------------------------------------------------------------------------
# parallel model

# formula string to run CFA
formula_string3 <- paste0(
  formula_string,
  "\n",
  paste(paste0("q", 1:42, ".t1"), collapse = "=="),
  "\n",
  paste(paste0("q", 1:42, ".t2"), collapse = "=="),
  "\n",
  paste(paste0("q", 1:42, ".t3"), collapse = "=="),
  "\n",
  paste(paste0("q", 1:42, ".t4"), collapse = "=="),
  "\n",
  paste(paste0("q", 1:42, ".t5"), collapse = "=="),
  collapse = ""
)

# standard deviation of the prior
t_sd <- t_prior_sd(
  N = 405,
  category = 6,
  n = 42
)

# fit a CFA model
fit.cfa3 <- cfa_sim(formula_string3, data, t_prior = t_sd)

# Plotting the wording-effect and parallel models
plts <- list()
thresholds1 <- t(apply(fit.cfa2$par_est[[2]][1:21,-1], 1, cut_trans))
thresholds2 <- t(apply(fit.cfa2$par_est[[2]][22:42,-1], 1, cut_trans))
plts[[1]] <- plot_DRM_compare(thresholds1[1,], thresholds2[1,], item_name = "Wording-effect Model", label = c("P", "N"))
thresholds1 <- t(apply(fit.cfa3$par_est[[2]][1:21,-1], 1, cut_trans))
plts[[2]] <- plot_DRM_compare(thresholds1[1,], thresholds1[1,], item_name = "    Parallel Model        ", label = c("P", "N"))

do.call(gridExtra::grid.arrange, c(plts, ncol=1))


# ------------------------------------------------------------------------------
# Model Selection

# log-likelihood
fit.cfa$logL
fit.cfa2$logL
fit.cfa3$logL
mirt::logLik(fit.gpcm)

# eff par
fit.cfa$eff_par
fit.cfa2$eff_par
fit.cfa3$eff_par
par_gpcm <- 6 * 42 # the correlation between the two latent variables is not counted

# AIC
fit.cfa$AIC
fit.cfa2$AIC
fit.cfa3$AIC
-2 * mirt::logLik(fit.gpcm) + 2 * par_gpcm

# BIC
fit.cfa$BIC
fit.cfa2$BIC
fit.cfa3$BIC
-2 * mirt::logLik(fit.gpcm) + log(mean(colSums(!is.na(data)))) * par_gpcm

# RMSEA
rmsea_drm(fit.cfa, fit.cfa2)
rmsea_drm(fit.cfa2, fit.cfa3)
