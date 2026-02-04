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
      0.61, -0.48,  8.79,
      0.65, -0.11,  4.85,
      0.86, -0.24,  9.56,
      0.82, -0.24,  7.51,
      0.45,  0.33,  4.96,
      0.57, -0.42,  6.05,
      0.53, -0.63,  5.65,
      0.48, -0.41,  5.61,
      0.68, -0.35,  5.45,
      0.72, -0.34,  7.75,
      0.88, -0.28, 11.35,
      0.50, -0.04,  6.28
    ),
    ncol = 3,
    byrow = TRUE
  )

  if(nitem==24) item <- rbind(item, item)
  if(nitem==48) item <- rbind(item, item, item, item)

  thresholds <- list()
  thresholds[[1]] <- c(0.2, 0.4, 0.6, 0.8)
  thresholds[[2]] <- c(0.22, 0.45, 0.55, 0.78)

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

for(N in c(500, 1000)){
  for(n_item in c(12, 24, 48)){
    {
      results <- list()
      bias1 <- bias2 <- rmse1 <- rmse2 <- 0
      
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
        fit <- efa_sim(data = data,
                       dimension = 1,
                       max_iter = 200,
                       t_prior = t_sd,
                       ngrid = 4000
        )
        
        results[[i]] <- fit$par_est
        cat("\n", i,"\n")
        
        # BIAS & RMSE
        item <- dataset$item
        
        thresholds <- list()
        thresholds[[1]] <- c(0.2, 0.4, 0.6, 0.8)
        thresholds[[2]] <- c(0.22, 0.45, 0.55, 0.78)
        
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

