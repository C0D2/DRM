library(tidyverse)
source("mDR_IRT.R")
Rcpp::sourceCpp("drm.cpp")

################################################################################
# Data generation by referring to the CFCS results
################################################################################
cfcs_simulation <- function(seed.number, N, nitem){
  dimension <- 1

  item <- matrix(
    c(
       0.73, -0.61, 5.96,
       0.70, -0.15, 4.00,
       0.96, -0.27, 7.60,
       0.89, -0.26, 6.10,
       0.50,  0.38, 3.68,
       0.64, -0.48, 4.64,
       0.58, -0.72, 4.37,
       0.49, -0.45, 5.32,
       0.73, -0.34, 4.39,
       0.80, -0.35, 6.00,
       0.98, -0.29, 8.91,
       0.52, -0.03, 5.59
    ),
    ncol = 3,
    byrow = TRUE
  )

  if(nitem==24) item <- rbind(item, item)
  if(nitem==48) item <- rbind(item, item, item, item)

  thresholds <- list()
  thresholds[[1]] <- c(0.2, 0.4, 0.6, 0.8)
  thresholds[[2]] <- c(0.21, 0.43, 0.54, 0.80)

  data <- matrix(nrow = N, ncol = nitem)
  set.seed(seed.number)
  mu <- rep(0, dimension)
  sigma <- diag(dimension)
  # sigma[1,2] <- sigma[2,1] <- 0.351

  theta <- mvtnorm::rmvnorm(n = N, mean = mu, sigma = sigma)
  for(i in 1:nitem){
    if(i %in% 1:(nitem/2)){
      temp_th <- thresholds[[1]]
    }else{
      temp_th <- thresholds[[2]]
    }
    ppp <- P_DRM(theta,item[i,1:dimension, drop=TRUE],item[i,dimension+1], item[i,dimension+2],cut_score = temp_th)
    for(j in 1:N){
      data[j,i] <- sample(0:length(temp_th), size = 1, prob = ppp[j,])
    }
  }
  colnames(data) <- paste0("q", 1:nitem)
  return(list(data=data,
              item=item))
}

################################################################################
# MODEL FITTING
################################################################################
n_rep <- 100

outcome_measures <- c()
model_selection_results <- c()
for(N in c(500, 1000)){
  for(n_item in c(12, 24, 48)){
    {
      results <- list()
      bias1 <- bias2 <- rmse1 <- rmse2 <- 0
      
      model_selection_AIC_BIC <- matrix(nrow = n_rep, ncol = 2)
      
      for(i in 1:n_rep){
        # data generation
        dataset <- cfcs_simulation(seed.number = i, N = N, nitem = n_item)
        data <- dataset$data
        
        # setting the prior distribution
        t_sd <- t_prior_sd(
          N = N,
          category = 5,
          n = 1
        )
        
        #  run analysis
        # Model 1: unrestricted model
        fit <- efa_sim(data = data,
                       dimension = 1,
                       max_iter = 200,
                       t_prior = t_sd,
                       ngrid = 4000
        )
        
        # Model 2: equal to equal
        model_formula <- paste("f1", paste(colnames(data), collapse = " + "), sep = " ~ ")
        
        item_names <- paste0("q", 1:(n_item/2))
        
        formula_string <- paste0(
          model_formula,"\n",
          paste(
            unlist(
              lapply(item_names, function(item) {
                sprintf("%s.t%d <- 0", item, 1:4)
              })
            ),
            collapse = "\n"
          ),
          collapse = "")
        
        fit1 <- cfa_sim(formula = formula_string,
                       data = data,
                       max_iter = 200,
                       t_prior = t_sd,
                       ngrid = 4000
        )
        
        # Model 3: unequal to equal
        item_names <- setdiff(colnames(data), item_names)
        
        formula_string <- paste0(
          model_formula,"\n",
          paste(
            unlist(
              lapply(item_names, function(item) {
                sprintf("%s.t%d <- 0", item, 1:4)
              })
            ),
            collapse = "\n"
          ),
          collapse = "")
        
        fit2 <- cfa_sim(formula = formula_string,
                        data = data,
                        max_iter = 200,
                        t_prior = t_sd,
                        ngrid = 4000
        )
        
        
        results[[i]] <- fit$par_est
        cat("\n", i,"\n")
        
        # Model selection
        # 0: the restricted model is preferred
        # 1: disagreement between AIC and BIC
        # 2: the unrestricted model is preferred
        AIC_ <- (fit1$AIC > fit$AIC)
        BIC_ <- (fit1$BIC > fit$BIC)
        
        model_selection_AIC_BIC[i,1] <- AIC_ + BIC_
        
        AIC_ <- (fit2$AIC > fit$AIC)
        BIC_ <- (fit2$BIC > fit$BIC)
        
        model_selection_AIC_BIC[i,2] <- AIC_ + BIC_
        
        # BIAS & RMSE
        item <- dataset$item
        
        thresholds <- list()
        thresholds[[1]] <- c(0.2, 0.4, 0.6, 0.8)
        thresholds[[2]] <- c(0.21, 0.43, 0.54, 0.80)
        
        item_par <- cbind(
          item[,1:2],
          log(item[,3]),
          matrix(c(rep(thresholds[[1]], n_item/2), rep(thresholds[[2]], n_item/2)),
                 ncol=4,
                 byrow = TRUE)
        )
        
        pars <- cbind(results[[i]][[1]],
                      log(results[[i]][[2]][,1]),
                      t(apply(results[[i]][[2]][,-1], 1, cut_trans, simplify = TRUE)))
        
        rmse <- (pars - item_par)^2
        bias <- (pars - item_par)
        bias1 <- bias1 + colMeans(bias[1:(n_item/2),])
        bias2 <- bias2 + colMeans(bias[(n_item/2+1):n_item,])
        rmse1 <- rmse1 + colMeans(rmse[1:(n_item/2),])
        rmse2 <- rmse2 + colMeans(rmse[(n_item/2+1):n_item,])
      }
      
      bias1 <- bias1/n_rep
      bias2 <- bias2/n_rep
      rmse1 <- sqrt(rmse1/n_rep)
      rmse2 <- sqrt(rmse2/n_rep)
      
      outcome_measures <- rbind(
        outcome_measures,
        c(N, n_item, "equal",
          paste(formatC(bias1, format="f", digits=3, flag=" "),
                " (",
                formatC(rmse1, format="f", digits=3, flag=""),
                ")",
                sep = "")
        ),
        c(N, n_item, "unequal",
          paste(formatC(bias2, format="f", digits=3, flag=" "),
                " (",
                formatC(rmse2, format="f", digits=3, flag=""),
                ")",
                sep = "")
        )
      )
      model_selection_results <- rbind(
        model_selection_results,
        c(
          sum(model_selection_AIC_BIC[,1] == 0),
          sum(model_selection_AIC_BIC[,1] == 1),
          sum(model_selection_AIC_BIC[,1] == 2),
          sum(model_selection_AIC_BIC[,2] == 0),
          sum(model_selection_AIC_BIC[,2] == 1),
          sum(model_selection_AIC_BIC[,2] == 2)
        )
      )
    }
  }
}



################################################################################
latex_matrix1 <- apply(outcome_measures[,1:6], 1, function(row) paste(row, collapse = " & "))
latex_matrix1 <- paste0("\\hline\n",
                        paste(latex_matrix1, collapse = " \\\\\n"),
                        "\n")
cat(latex_matrix1)

latex_matrix1 <- apply(outcome_measures[,c(1:3,7:10)], 1, function(row) paste(row, collapse = " & "))
latex_matrix1 <- paste0("\\hline\n",
                        paste(latex_matrix1, collapse = " \\\\\n"),
                        "\n")
cat(latex_matrix1)

