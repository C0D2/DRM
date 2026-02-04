################################################################################
# HELPER FUNCTIONS
################################################################################
t_prior_sd <- function(N, category, n){
  1/sqrt(N/category*2*n)*sqrt(4*pi^2/3)
}

one_hot_3d <- function(M) {
  N <- nrow(M)            # number of persons
  I <- ncol(M)            # number of items
  K <- max(M, na.rm = TRUE) + 1         # number of categories (0-based)

  arr <- array(0, dim = c(I, N, K))

  for (i in 1:I) {
    mat <- matrix(0, N, K)
    mat[cbind(1:N, M[,i] + 1)] <- 1   # +1 because R indices start at 1
    arr[i,,] <- mat
  }

  return(arr)
}

jac_soft <- function(t) {
  npar <- length(t)
  e <- c(1, exp(t))
  D <- sum(e)

  N <- rep(0, npar)
  tmp <- 0
  for(i in 1:npar){
    tmp <- tmp + e[i]
    N[i] <- tmp
  }

  J <- matrix(0, npar, npar)

  for (i in 1:npar) {
    for (j in 1:npar) {
      if (j < i) {
        J[i,j] <- e[j] * (D - N[i]) / D^2
      } else {
        J[i,j] <- - N[i] * e[j] / D^2
      }
    }
  }
  return(J)
}

count_cat <- function(x) length(unique(x))
extract_cat <- function(x) sort(unique(x))

reorder_vec <- function(x){
  match(x, table = extract_cat(x))-1
}

reorder_mat <- function(x) apply(x, MARGIN = 2, FUN = reorder_vec)

P_2PL <- function(theta, a, b) 1/(1 + exp(-(as.vector(theta %*% a) + b)))

cut_trans <- function(x){
  x <- x[!is.na(x)]
  return(
    cumsum(
      c(1,exp(x[-length(x)]))/sum(c(1,exp(x)))
    )
  )
}

P_DRM <- function(theta, a, b, nu, ncats=NULL, cut_score=NULL, return_mu = FALSE){
  p <- P_2PL(theta = theta, a = a, b = b)
  if(is.null(cut_score) & is.null(ncats)){
    stop("Specify either ncat or cut_score.")
  }else if(is.null(cut_score)){
    cut_score <- (1:(ncats-1))/ncats
  }else{
    ncats <- length(cut_score)+1
  }

  probs <- matrix(nrow = nrow(theta), ncol = ncats-1)

  for(i in 1:length(cut_score)){
    probs[,i] <- pbeta(q = cut_score[i],
                       shape1 = p*nu,
                       shape2 = (1-p)*nu)
  }
  if(return_mu){
    return(
      list(prob = cbind(probs,1)-cbind(0,probs),
           mu = p)
    )
  }else {
    return(cbind(probs,1)-cbind(0,probs))
  }
}


llik <- function(data, theta, item, cut_score){
  d <- ncol(theta)
  p_matrix <- P_DRM(theta, item[1:d], item[d+1], item[d+2], cut_score = cut_score)
  likelihood <- sapply(1:nrow(theta), function(j) p_matrix[j, data+1])
  likelihood[is.na(likelihood)] <- 1
  likelihood <- log(likelihood)
  likelihood[likelihood==-Inf] <- -.Machine$double.xmax
  return(likelihood)
}

################################################################################
# RMSEA
################################################################################
rmsea_drm <- function(fit.ref, fit.comp){
  sample_size <- mean(colSums(!is.na(fit.ref$Options$data)))
  df <- fit.ref$eff_par - fit.comp$eff_par
  
  lrt <- -2*(fit.comp$logL - fit.ref$logL)
  lambda <- lrt - df # non-centrality parameter
  f0 <- max(lambda / (sample_size), 0)
  
  if(lambda > 0){
    ci90 <- (qchisq(c(0.05, 0.95), df = df, ncp = lambda) - df)/ (df * sample_size)
  } else {
    ci90 <- NULL
  }
  
  message("\r chi-square: ", sprintf("%.2f", lrt),",  df: ", sprintf("%.2f", df))
  return(
    list(
      rmsea = sqrt(f0/df),
      ci90 = sqrt(ci90)
    )
  )
}

################################################################################
# COMMUNALITY & UNIQUE DISPERSION
################################################################################
null_dispersion <- function(data, cut_score=NULL){
  data <- reorder_vec(data)
  freq <- as.vector(table(data))

  ncats <- length(freq)

  if(is.null(cut_score)) {
    cut_score <- (1:ncats)/ncats
  } else {
    cut_score <- c(cut_score, 1)
  }

  pars <- c(0.5, 1)
  iter <- 0
  repeat{
    iter <- iter + 1

    mu <- pars[1]
    xi <- pars[2]
    nu <- exp(xi)

    p0 <- pbeta(q = cut_score, shape1 = mu*nu,shape2 = (1-mu)*nu)
    p0 <- p0 - c(0, p0[-length(p0)])

    ngrid <- 1000

    beta_grid <- seq(1/2/ngrid, 1-1/2/ngrid, length=ngrid)

    ind_cat <- as.numeric(cut(beta_grid,breaks = c(0,cut_score),labels = 1:ncats))

    p_ <- dbeta(beta_grid, shape1 = nu*mu, shape2 = nu*(1-mu))
    s1 <- log(beta_grid) - digamma(nu*mu)
    s2 <- log(1- beta_grid) - digamma(nu*(1-mu))
    num_mu <- p_*nu*(s1 - s2)
    num_xi <- nu*p_*(digamma(nu) + mu*s1 + (1-mu)*s2)


    l1m <- l1x <- c()
    for(ct in 1:ncats){
      l1m[ct] <- sum(num_mu[ind_cat==ct])/length(beta_grid)
      l1x[ct] <- sum(num_xi[ind_cat==ct])/length(beta_grid)
    }

    L1 <- c(
      freq %*% (l1m/p0),
      freq %*% (l1x/p0)
    )

    L2 <- sum(freq) * t(cbind(l1m, l1x)) %*% diag(1/p0) %*% cbind(l1m, l1x)

    diff <- L1 %*% solve(L2)

    pars <- pars + diff
    if(sum(abs(diff)) < 0.000001 | iter > 30) break
  }
  pars[2] <- exp(pars[2])

  return(as.vector(pars))
}

communality <- function(fit){
  if(inherits(fit, "sg")){
    thresholds <- apply(fit$par_est[[2]][,-1], 1, cut_trans, simplify = FALSE)

    results <- matrix(nrow = ncol(fit$Options$data), ncol = 5)
    for(i in 1:nrow(results)){
      results[i, 1:2] <- null_dispersion(fit$Options$data[,i], thresholds[[i]])
    }

    results[, 3] <- fit$par_est[[2]][,1]

    results[, 5] <- results[, 2]/results[, 3]
    results[, 4] <- 1 - results[, 5]

    colnames(results) <- c("mu", "nu_null", "nu_model", "communality", "unique dispersion")
  } else if(inherits(fit, "mg")){
    group <- fit$Options$group
    n_group <- length(unique(group))
    n_item <- ncol(fit$Options$data)

    results <- list()

    for(g in 1:n_group){
      index <- (1:n_item) + n_item*(g-1)
      thresholds <- apply(fit$par_est[[2]][index,-1], 1, cut_trans, simplify = FALSE)

      tmp <- matrix(nrow = ncol(fit$Options$data), ncol = 5)
      for(i in 1:nrow(tmp)){
        tmp[i, 1:2] <- null_dispersion(fit$Options$data[group == g,i], thresholds[[i]])
      }

      tmp[, 3] <- fit$par_est[[2]][index,1]

      tmp[, 5] <- tmp[, 2]/tmp[, 3]
      tmp[, 4] <- 1 - tmp[, 5]

      colnames(tmp) <- c("mu", "nu_null", "nu_model", "communality", "unique dispersion")

      results[[g]] <- tmp
    }
  }


  return(results)
}

################################################################################
# CFA HELPER FUNCTIONS
################################################################################
library(igraph)

build_loading_matrix <- function(formula_string, variable_order = NULL) {
  # Split the formula string into lines
  lines <- unlist(strsplit(formula_string, "\n"))
  lines <- trimws(lines)  # remove leading/trailing whitespace
  lines <- lines[lines != ""]  # remove empty lines
  lines <- lines[grepl("^f\\w*\\s*~", lines)]

  # Initialize list to store mappings
  mapping <- list()

  for (line in lines) {
    # Split formula into LHS and RHS
    parts <- strsplit(line, "~")[[1]]
    factor <- trimws(parts[1])
    variables <- unlist(strsplit(parts[2], "\\+"))
    variables <- trimws(variables)

    for (var in variables) {
      mapping[[var]] <- c(mapping[[var]], factor)
    }
  }

  # Determine variable and factor order
  all_vars <- unique(names(mapping))
  all_factors <- sort(unique(unlist(mapping)))

  # Use user-specified order if given
  if (!is.null(variable_order)) {
    if (!all(all_vars %in% variable_order)) {
      stop("Some variables in the formula are not in the data.")
    }
    all_vars <- variable_order
  } else {
    all_vars <- sort(all_vars)
  }

  # Create matrix
  mat <- matrix(0, nrow = length(all_vars), ncol = length(all_factors),
                dimnames = list(all_vars, all_factors))

  # Fill in 1s where variables load onto factors
  for (var in names(mapping)) {
    for (fac in mapping[[var]]) {
      mat[var, fac] <- 1
    }
  }

  return(mat)
}
build_constraint_matrix <- function(formula_string, variable_order, factor_order) {
  # Split and trim lines
  lines <- unlist(strsplit(formula_string, "\n"))
  lines <- trimws(lines)
  constraint_lines <- lines[grepl("==", lines)]

  # Initialize zero matrix
  mat <- matrix(0, nrow = length(variable_order), ncol = length(factor_order),
                dimnames = list(variable_order, factor_order))

  # Assign unique constraint group IDs
  for (i in seq_along(constraint_lines)) {
    terms <- trimws(unlist(strsplit(constraint_lines[i], "==")))
    for (term in terms) {
      split_term <- strsplit(term, "\\.")[[1]]
      var <- split_term[1]
      fac <- split_term[2]
      mat[var, fac] <- i
    }
  }

  return(mat)
}
apply_fixed_values <- function(formula_string, init_matrix, loading_matrix) {
  # Extract and trim assignment lines
  lines <- unlist(strsplit(formula_string, "\n"))
  lines <- trimws(lines)
  assign_lines <- lines[grepl("<-", lines, fixed = TRUE)]

  for (line in assign_lines) {
    parts <- unlist(strsplit(line, "<-", fixed = TRUE))
    term <- trimws(parts[1])  # e.g., "r4.f1"
    value <- as.numeric(trimws(parts[2]))  # e.g., 3

    split_term <- strsplit(term, "\\.")[[1]]
    var <- split_term[1]
    fac <- split_term[2]

    # Set value and mark as fixed
    if (var %in% rownames(init_matrix) && fac %in% colnames(init_matrix)) {
      init_matrix[var, fac] <- value
      loading_matrix[var, fac] <- 0
    } else {
      warning(sprintf("Invalid variable/factor name: %s.%s", var, fac))
    }
  }

  return(list(initial = init_matrix, loading = loading_matrix))
}

apply_equal_constraints_seq <- function(param_id_mat, eq_mat) {
  stopifnot(all(dim(param_id_mat) == dim(eq_mat)))

  rn <- rownames(param_id_mat)
  cn <- colnames(param_id_mat)

  # Convert to long form
  df <- as.data.frame(as.table(param_id_mat))
  colnames(df) <- c("row", "col", "param_id")
  df$eq <- as.vector(eq_mat)

  # Only keep cells with non-NA param_id
  df <- df[!is.na(df$param_id), ]

  # Track groups based on equal constraints
  edge_list <- do.call(rbind, lapply(split(df, df$eq), function(group) {
    ids <- unique(group$param_id)
    if (group$eq[1] != 0 && length(ids) > 1) t(combn(ids, 2)) else NULL
  }))

  # Build graph of equal constraints
  if (!is.null(edge_list)) {
    g <- igraph::graph_from_edgelist(apply(edge_list, 2, as.character), directed = FALSE)
  } else {
    g <- igraph::make_empty_graph()
  }

  all_ids <- unique(as.character(df$param_id))
  g <- igraph::add_vertices(g, nv = length(setdiff(all_ids, igraph::V(g)$name)),
                            name = setdiff(all_ids, igraph::V(g)$name))

  comps <- igraph::components(g)
  grouped_ids <- split(as.integer(igraph::V(g)$name), comps$membership)

  # Assign sequential new IDs to grouped values
  id_map <- list()
  next_id <- 1
  for (group in grouped_ids) {
    for (id in group) {
      id_map[[as.character(id)]] <- next_id
    }
    next_id <- next_id + 1
  }

  # Assign new IDs to ungrouped param_ids (not in constraint graph)
  remaining_ids <- setdiff(as.character(df$param_id), names(id_map))
  for (id in remaining_ids) {
    id_map[[id]] <- next_id
    next_id <- next_id + 1
  }

  # Rebuild the output matrix
  out_mat <- param_id_mat
  for (i in seq_len(nrow(out_mat))) {
    for (j in seq_len(ncol(out_mat))) {
      id <- param_id_mat[i, j]
      if (!is.na(id)) {
        out_mat[i, j] <- id_map[[as.character(id)]]
      }
    }
  }

  rownames(out_mat) <- rn
  colnames(out_mat) <- cn
  return(out_mat)
}

group_parameters <- function(df) {
  # Flatten all values, get unique parameter IDs
  vals <- na.omit(as.numeric(unlist(df)))
  parent <- setNames(as.list(vals), vals)  # union-find parent map

  # Union-find helpers
  find <- function(x) {
    while (parent[[as.character(x)]] != x) {
      parent[[as.character(x)]] <- parent[[as.character(parent[[as.character(x)]])]]
      x <- parent[[as.character(x)]]
    }
    x
  }

  union <- function(x, y) {
    root_x <- find(x)
    root_y <- find(y)
    if (root_x != root_y) {
      parent[[as.character(root_y)]] <<- root_x
    }
  }

  # For each row, union all values
  apply(df, 1, function(row) {
    row_vals <- na.omit(as.numeric(row))
    if (length(row_vals) > 1) {
      for (i in 2:length(row_vals)) {
        union(row_vals[1], row_vals[i])
      }
    }
  })

  # Group by connected component (root)
  groups <- list()
  for (v in vals) {
    root <- find(v)
    root_str <- as.character(root)
    if (!root_str %in% names(groups)) {
      groups[[root_str]] <- c()
    }
    groups[[root_str]] <- c(groups[[root_str]], v)
  }

  # Clean and sort
  groups <- lapply(groups, function(x) sort(unique(x)))
  names(groups) <- paste0("[[", seq_along(groups), "]]")
  return(groups)
}


efa_helper_matrices <- function(d, data, eq_interval = FALSE){
  if(is.null(colnames(data))){
    cns <- paste0("v", 1:ncol(data))
  } else {
    cns <- colnames(data)
  }

  formula_string <- paste(paste(paste0("\n f",1:d), paste(cns, collapse = " + "), sep = " ~ "), collapse = "\n")

  # initial matrix for thresholds
  data <- reorder_mat(data)
  category <- apply(data, 2, max, na.rm = TRUE)
  init_mat <- matrix(nrow = ncol(data), ncol = max(category)+1)
  init_mat[,1] <- 2
  for(i in 1:nrow(init_mat)){
    init_mat[i, 2:(category[i]+1)] <- 0
  }

  init_mat <- cbind(matrix(rep(c((d:1)/d), each=ncol(data)), nrow=ncol(data)),
                    as.numeric(scale(colMeans(data, na.rm = TRUE)/category, center = TRUE, scale = TRUE)/2),
                    init_mat)
  load_mat <- build_loading_matrix(formula_string, variable_order = cns)
  load_mat[,1:d][upper.tri(load_mat[,1:d])] <- 0
  fns <- c(colnames(load_mat), "c", "nu", paste0("t", 1:max(category)))
  load_mat <- cbind(load_mat,1,1,1*(!is.na(init_mat[,-(1:(d+2))])))
  init_mat <- init_mat * load_mat
  par_id <- matrix(1:length(as.vector(init_mat)), ncol = ncol(init_mat))
  dimnames(init_mat) <- list(cns, fns)
  dimnames(load_mat) <- list(cns, fns)
  dimnames(par_id) <- list(cns, fns)

  eq_constraint <- build_constraint_matrix(formula_string, cns, fns)

  fxd <- apply_fixed_values(formula_string,init_matrix = init_mat, load_mat)
  init_mat <- fxd$initial
  load_mat <- fxd$loading

  if(eq_interval == TRUE) load_mat[, (d+3):ncol(load_mat)] <- 0

  par_id <- par_id * load_mat
  par_id[par_id == 0] <- NA

  par_id <- apply_equal_constraints_seq(par_id, eq_constraint)

  grouping_index <- group_parameters(par_id)

  return(list(
    init_mat = init_mat,
    load_mat = load_mat,
    eq_constraint = eq_constraint,
    par_id = par_id,
    grouping_index = grouping_index
  ))
}

cfa_helper_matrices <- function(formula_string, data, eq_interval = FALSE){
  if(is.null(colnames(data))){
    cns <- paste0("v", 1:ncol(data))
  } else {
    cns <- colnames(data)
  }

  load_mat <- build_loading_matrix(formula_string, variable_order = cns)
  d <- ncol(load_mat)
  # initial matrix for thresholds
  data <- reorder_mat(data)
  category <- apply(data, 2, max, na.rm = TRUE)
  fns <- c(colnames(load_mat), "c", "nu", paste0("t", 1:max(category)))

  init_mat2 <- matrix(nrow = ncol(data), ncol = max(category)+1)
  init_mat2[,1] <- 2
  for(i in 1:nrow(init_mat2)){
    init_mat2[i, 2:(category[i]+1)] <- 0
  }


  init_mat <- cbind(matrix(rep(1/d, each=ncol(data)*d), nrow=ncol(data)),
                    as.numeric(scale(colMeans(data, na.rm = TRUE)/category, center = TRUE, scale = TRUE)/2)
  )* cbind(load_mat, 1)
  init_mat <- cbind(init_mat, init_mat2)

  load_mat <- cbind(load_mat,1,1,1*(!is.na(init_mat[,-(1:(d+2))])))
  par_id <- matrix(1:length(as.vector(init_mat)), ncol = ncol(init_mat))
  dimnames(init_mat) <- list(cns, fns)
  dimnames(load_mat) <- list(cns, fns)
  dimnames(par_id) <- list(cns, fns)

  eq_constraint <- build_constraint_matrix(formula_string, cns, fns)

  fxd <- apply_fixed_values(formula_string,init_matrix = init_mat, load_mat)
  init_mat <- fxd$initial
  load_mat <- fxd$loading

  if(eq_interval == TRUE) load_mat[, (d+3):ncol(load_mat)] <- 0

  par_id <- par_id * load_mat
  par_id[par_id == 0] <- NA

  par_id <- apply_equal_constraints_seq(par_id, eq_constraint)

  grouping_index <- group_parameters(par_id)

  return(list(
    init_mat = init_mat,
    load_mat = load_mat,
    eq_constraint = eq_constraint,
    par_id = par_id,
    grouping_index = grouping_index,
    d = d
  ))
}


group_helper_matrices <- function(formula_string, data, eq_interval = FALSE, ngroup){
  if(is.null(colnames(data))){
    cns0 <- paste0("v", 1:ncol(data))
  } else {
    cns0 <- colnames(data)
  }
  cns <- c()
  for(i in 1:ngroup){
    cns <- append(cns, paste0(cns0, paste0("g", i)))
  }

  load_mat <- build_loading_matrix(formula_string, variable_order = cns)
  d <- ncol(load_mat)

  # initial matrix for thresholds
  data <- reorder_mat(data)
  category <- apply(data, 2, max, na.rm = TRUE)
  category <- rep(category, ngroup)
  fns <- c(colnames(load_mat), "c", "nu", paste0("t", 1:max(category)))

  init_mat2 <- matrix(nrow = nrow(load_mat), ncol = max(category)+1)
  init_mat2[,1] <- 2
  for(i in 1:nrow(init_mat2)){
    init_mat2[i, 2:(category[i]+1)] <- 0
  }


  init_mat <- cbind(matrix(rep(1/d, each=nrow(load_mat)*d), nrow=nrow(load_mat)),
                    as.numeric(scale(colMeans(data, na.rm = TRUE)/category, center = TRUE, scale = TRUE)/2)
  )* cbind(load_mat, 1)
  init_mat <- cbind(init_mat, init_mat2)

  load_mat <- cbind(load_mat,1,1,1*(!is.na(init_mat[,-(1:(d+2))])))
  par_id <- matrix(1:length(as.vector(init_mat)), ncol = ncol(init_mat))
  dimnames(init_mat) <- list(cns, fns)
  dimnames(load_mat) <- list(cns, fns)
  dimnames(par_id) <- list(cns, fns)

  eq_constraint <- build_constraint_matrix(formula_string, cns, fns)

  fxd <- apply_fixed_values(formula_string,init_matrix = init_mat, load_mat)
  init_mat <- fxd$initial
  load_mat <- fxd$loading

  if(eq_interval == TRUE) load_mat[, (d+3):ncol(load_mat)] <- 0

  par_id <- par_id * load_mat
  par_id[par_id == 0] <- NA

  par_id <- apply_equal_constraints_seq(par_id, eq_constraint)

  grouping_index <- group_parameters(par_id)

  return(list(
    init_mat = init_mat,
    load_mat = load_mat,
    eq_constraint = eq_constraint,
    par_id = par_id,
    grouping_index = grouping_index,
    d = d
  ))
}
################################################################################
# EM & MH-RM ALGORITHMS
################################################################################
sample_mhrm <- function(data, theta_c, item, cut_score, f_cov=NULL, sd=1){
  d <- ncol(theta_c)
  if(is.null(f_cov)){
    Sigma_inv <- diag(d)
  }else{
    Sigma_inv <- solve(f_cov)
  }

  AR <- 0

  theta_p <- theta_c
  for(i in 1:d){
    theta_p[,i] <- theta_p[,i] + rnorm(nrow(theta_p), 0, sd)
    lr <- pmin(llik_ratio_mhrm(data, theta_c, theta_p, item, cut_score, Sigma_inv), 0)
    not_accepted <- (log(runif(nrow(theta_p))) > lr)

    theta_p[not_accepted,i] <- theta_c[not_accepted,i]
    theta_c <- theta_p

    AR <- AR + (1 - mean(not_accepted)) / d
  }
  
  logL <- llik_mhrm(data, theta_c, item, cut_score, Sigma_inv)
  return(
    list(
      theta=theta_c,
      logL=logL,
      AR=AR
      )
  )
}

Mstep_sim <- function(E, item, par_id=NULL, grouping=NULL, ngrid = 1000,
                  max_iter=7, threshold=1e-7, prior, calculate_m=FALSE){
  grid <- E$grid
  d <- ncol(grid)

  n_item <- nrow(item[[1]])
  npar <-  max(par_id, na.rm = TRUE)
  se <- list()
  f_mat <- list()

  e.response <- E$e.response
  q <- nrow(grid)
  beta_grid <- seq(1/2/ngrid, 1-1/2/ngrid, length=ngrid)

  iter <- 0
  div <- 3

  repeat{
    iter <- iter + 1

    L1 <- rep(0, npar)
    L2 <- matrix(0, nrow = npar, ncol = npar)
    IMm <- matrix(0, nrow = npar, ncol = npar)
    se <- matrix(0, nrow = npar, ncol = npar)
    inv_L2 <- matrix(0, nrow = npar, ncol = npar)
    l_prior_w <- 0
    diff <- L1
    grad_norm <- 0
    eff_par <- 0

    # stacking up the gradients and information
    for(i in 1:n_item){
      L1L2 <- L1L2_sim(c(item[[1]][i,],item[[2]][i,]),
                       e.response[i,,],
                       grid, beta_grid, d, calculate_m)

      ind <- par_id[i,]
      ind2 <- !is.na(ind)
      ind <- ind[ind2]

      L1[ind] <- L1[ind] + L1L2$Grad[ind2]
      L2[ind,ind] <- L2[ind,ind] + L1L2$IM[ind2,ind2]
      if(calculate_m){
        IMm[ind,ind] <- IMm[ind,ind] + L1L2$IMm[ind2,ind2]
      }
    }
    L1[!is.finite(L1)] <- sign(L1[!is.finite(L1)]) * 1
    original_L2 <- L2 # to later calculate the effective number of parameters

    # apply prior
    partial_id <- par_id[,(d+3):ncol(par_id)]
    item_t <- item[[2]][,-1]
    ind <- as.vector(partial_id)
    unique_ind <- unique(na.omit(ind))
    current_t <- sapply(unique_ind, function(i){
      mean(item_t[ind==i], na.rm=TRUE)
    })
    # item_number <- which(partial_id %in% unique_ind, arr.ind = TRUE)
    # item_number <- arrayInd(item_number, .dim = dim(partial_id))[, 1]
    item_number <- sapply(unique_ind, function(v) {
      which(partial_id == v, arr.ind = TRUE)[1, 1]
    })
    if(length(current_t) != 0){
      l_prior_w <- mvtnorm::dmvnorm(current_t, mean = rep(0, length(current_t)), sigma = diag(prior[item_number]^2), log = TRUE) # log-prior-weight
      L1[unique_ind] <- L1[unique_ind] - (current_t - 0)/(prior[item_number]^2)
      L2[unique_ind,unique_ind] <- L2[unique_ind,unique_ind] + diag(1/(prior[item_number]^2))
    }

    # calculate parameter changes
    for(i in 1:length(grouping)){
      inv_L2[grouping[[i]],grouping[[i]]] <- solve(L2[grouping[[i]],grouping[[i]]])
      if(calculate_m) se[grouping[[i]],grouping[[i]]] <- solve(
        L2[grouping[[i]],grouping[[i]]] - IMm[grouping[[i]],grouping[[i]]] + L1[grouping[[i]]] %*% t(L1[grouping[[i]]])
        )

      diff[grouping[[i]]] <- as.vector(inv_L2[grouping[[i]],grouping[[i]]] %*% L1[grouping[[i]]])

      grad_norm <- grad_norm + diff[grouping[[i]]] %*% L1[grouping[[i]]]
      eff_par <- eff_par + sum(diag(original_L2[grouping[[i]],grouping[[i]]] %*% inv_L2[grouping[[i]],grouping[[i]]]))
    }

    # convert diff-vector to matrix
    diff <- diff_to_matrix(par_id, diff)

    # update parameters
    if( max(abs(diff)) > div){
      stepsize <- max(sum(abs(diff)),2)
      item[[1]] <- item[[1]] + diff[,1:(d+1)]/stepsize
      item[[2]] <- item[[2]] + diff[,(d+2):ncol(diff)]/stepsize
    } else {
      item[[1]] <- item[[1]] + diff[,1:(d+1)]/2
      item[[2]] <- item[[2]] + diff[,(d+2):ncol(diff)]/2
      div <- max(abs(diff))
    }

    if( div <= threshold | iter > max_iter) break
  }
  if(calculate_m){
    se_delta <- se
    if(length(current_t) != 0){
      for(i in 1:nrow(partial_id)){
        index <- partial_id[i,]
        J <- jac_soft(item[[2]][i,-1])
        se_delta[index,index] <- J %*% se[index,index] %*% t(J)
      }
    }
    I_o <- list()
    for(i in 1:length(grouping)){
      I_o[[i]] <- se[grouping[[i]],grouping[[i]]]
    }
    suppressWarnings({
      se_delta <- diff_to_matrix(par_id, sqrt(diag(se_delta)))
      se_raw <- diff_to_matrix(par_id, sqrt(diag(se)))
    })
    return(list(
      # inv_H=inv_L2,
      I_o=I_o,
      se=se_delta,
      se_raw=se_raw
      ))
  } else {
    return(list(
      item=item,
      grad_norm=grad_norm,
      eff_par=eff_par,
      l_prior_w=l_prior_w,
      H=original_L2,
      H_prior_added=L2))
  }
}

update_gH <- function(data_array, thetalist, item, par_id=NULL, grouping=NULL, ngrid = 1000,
                      prev_g, prev_H, prior, stepsize=1, f_cov=NULL){
  d <- ncol(thetalist[[1]])

  if(is.null(f_cov)){
    f_cov <- diag(d)
  }

  n_item <- nrow(item[[1]])
  npar <-  max(par_id, na.rm = TRUE)

  beta_grid <- seq(1/2/ngrid, 1-1/2/ngrid, length=ngrid)

  n_sample <- length(thetalist)

  # START the g, H calculation for each item
  L1 <- rep(0, npar)
  L2 <- matrix(0, nrow = npar, ncol = npar)
  diff <- L1
  grad_norm <- 0
  eff_par <- 0
  mu_hat <- score_Sigma <- H_Sigma <- 0

  # stacking up the gradients and information
  for(k in 1:n_sample){
    theta <- thetalist[[k]]
    for(i in 1:n_item){
      L1L2 <- L1L2_sim(c(item[[1]][i,],item[[2]][i,]),
                       data_array[i,,],
                       theta, beta_grid, d)

      ind <- par_id[i,]
      ind2 <- !is.na(ind)
      ind <- ind[ind2]

      L1[ind] <- L1[ind] + L1L2$Grad[ind2] / n_sample
      L2[ind,ind] <- L2[ind,ind] + L1L2$IM[ind2,ind2] / n_sample
    }
    # L1[!is.finite(L1)] <- sign(L1[!is.finite(L1)]) * 1

    n <- nrow(theta)
    mu_hat <- mu_hat + colMeans(theta) / n_sample
    Sigma_hat <- cov(theta)
    Sigma_inv <- solve(f_cov)
    score_Sigma <- score_Sigma + (n/2) * as.vector(Sigma_inv - Sigma_inv %*% Sigma_hat %*% Sigma_inv) / n_sample
    H_Sigma <- H_Sigma -(n/2) * kronecker(Sigma_inv, Sigma_inv) / n_sample
  }

  # update gradient and Hessian
  prev_g[[1]] <- prev_g[[1]] + stepsize * (L1 - prev_g[[1]])
  L2 <- prev_H[[1]] <- prev_H[[1]] + stepsize * (L2 - prev_H[[1]])
  original_L2 <- L2 # to later calculate the effective number of parameters

  prev_g[[2]] <- prev_g[[2]] + stepsize * (score_Sigma - prev_g[[2]])
  prev_H[[2]] <- prev_H[[2]] + stepsize * (H_Sigma - prev_H[[2]])
  f_cov <- f_cov + stepsize * solve(prev_H[[2]], score_Sigma)
  # sds <- sqrt(diag(f_cov))
  # f_cov <- f_cov / (sds %o% sds)



  # apply prior
  partial_id <- par_id[,(d+3):ncol(par_id)]
  item_t <- item[[2]][,-1]
  ind <- as.vector(partial_id)
  unique_ind <- unique(na.omit(ind))
  current_t <- sapply(unique_ind, function(i){
    mean(item_t[ind==i], na.rm=TRUE)
  })

  item_number <- sapply(unique_ind, function(v) {
    which(partial_id == v, arr.ind = TRUE)[1, 1]
  })
  if(length(current_t) != 0){
    L1[unique_ind] <- L1[unique_ind] - (current_t - 0)/(prior[item_number]^2)
    L2[unique_ind,unique_ind] <- L2[unique_ind,unique_ind] + diag(1/(prior[item_number]^2))
  }

  # calculate parameter changes
  for(i in 1:length(grouping)){
    inv_L2 <- solve(L2[grouping[[i]],grouping[[i]]])
    diff[grouping[[i]]] <- as.vector(inv_L2 %*% L1[grouping[[i]]])

    grad_norm <- grad_norm + diff[grouping[[i]]] %*% L1[grouping[[i]]]
    eff_par <- eff_par + sum(diag(original_L2[grouping[[i]],grouping[[i]]] %*% inv_L2))
  }

  # convert diff-vector to matrix
  diff <- diff_to_matrix(par_id, diff)

  # update parameters
  item[[1]] <- item[[1]] + stepsize * diff[,1:(d+1)]
  item[[2]] <- item[[2]] + stepsize * diff[,(d+2):ncol(diff)]

  return(list(
    item=item,
    g_approx=prev_g,
    H_approx=prev_H,
    L1=L1,
    H=original_L2,
    H_prior_added=L2,
    grad_norm=grad_norm,
    eff_par=eff_par,
    diff=diff,
    factor_means=mu_hat,
    f_cov=f_cov))
}

L1L2_sim <- function(par, e.response, grid, beta_grid, d, calculate_m = FALSE){
  # par <- par[!is.na(par)]
  npar <- length(par)
  cut_score <- cut_trans(par[(d+3):npar])
  f <- rowSums(e.response)
  ncats <- length(cut_score)+1
  ncut <- ncats-1
  ind_cat <- as.numeric(cut(beta_grid,breaks = c(0,cut_score,1),labels = 1:ncats))

  # 1st derivatives
  l1s <- rep(list(0), npar)

  nu <- exp(par[d+2])

  p0 <- P_DRM_cpp(grid, par[1:d], par[d+1], nu, cut_score = cut_score, return_mu = TRUE)
  mu <- p0$mu
  p0 <- p0$prob

  # pmat <- t(outer(beta_grid, nu*mu-1, FUN = "^")*outer(1-beta_grid, nu*(1-mu)-1, FUN = "^"))/beta(nu*mu,nu*(1-mu)) # probability matrix wo the normalizing factor
  # l1mu <- nu*pmat*t(outer(log(beta_grid/(1-beta_grid)), digamma(nu*mu)-digamma(nu*(1-mu)), FUN = "-"))
  # l1xi <- nu*(digamma(nu)+
  #               mu*t(outer(log(beta_grid), digamma(nu*mu), FUN = "-"))+
  #               (1-mu)*t(outer(log(1-beta_grid), digamma(nu*(1-mu)), FUN = "-"))
  # )*pmat
  # l1m <- matrix(ncol = ncats, nrow = nrow(grid))
  # l1x <- matrix(ncol = ncats, nrow = nrow(grid))
  # for(ct in 1:ncats){
  #   l1m[,ct] <- rowSums(l1mu[,ind_cat==ct])/length(beta_grid)
  #   l1x[,ct] <- rowSums(l1xi[,ind_cat==ct])/length(beta_grid)
  # }
  # l1m <- mu*(1-mu)*l1m
  # for(i in 1:(d+1)){
  #   l1s[[i]] <- sweep(l1m, 1, cbind(grid,1)[,i], FUN = "*")
  # }
  # l1s[[d+2]] <- l1x
  l1s <- compute_l1_cpp(grid, nu, mu, beta_grid, ncats, ind_cat)

  # thresholds
  beta_vals <- matrix(nrow = nrow(grid), ncol=ncut)
  for(b in 1:ncut){
    beta_vals[,b] <- dbeta(cut_score[b], shape1 = nu*mu, shape2 = nu*(1-mu))
  }

  for(b in 1:ncut){
    ind <- c(rep(0,b),rep(1,ncut-b))
    temp <- t((ind-cut_score)*t(beta_vals))*cut_score[1]*exp(par[b+1])
    l1s[[b+d+2]] <- cbind(temp,0)-cbind(0,temp)
  }

  # derivatives
  # Grad <- c()
  # IM <- matrix(0, nrow = npar, ncol = npar)
  # for(r in 1:npar){
  #   Grad <- append(Grad, sum(l1s[[r]]*e.response[,1:ncats]/p0))
  #   for(c in 1:npar){
  #     if(r>=c){
  #       IM[r,c] <- rowSums(l1s[[r]]*l1s[[c]]/p0)%*%f
  #       IM[c,r] <- IM[r,c]
  #     }
  #   }
  # }
  res_cpp <- Grad_IM_cpp(l1s, e.response, p0, f, calculate_m)

  if(calculate_m){
    return(list(
      Grad = res_cpp$Grad,
      IM = res_cpp$IM,
      f = sum(f),
      IMm = res_cpp$IMm
    ))
  }else{
    return(list(
      Grad = res_cpp$Grad,
      IM = res_cpp$IM,
      f = sum(f)
    ))
  }
}

diff_to_matrix <- function(par_id, vec) {
  mat <- as.matrix(par_id)
  replaced <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat),
                     dimnames = dimnames(mat))

  for (i in seq_along(vec)) {
    replaced[mat == i] <- vec[i]
  }

  replaced[is.na(replaced)] <- 0
  return(replaced)
}

################################################################################
# THE MAIN DRIVER FUNCTION
################################################################################
dr_sim <- function(data, dimension=NULL, range=c(-4,4), q=11, t_prior=0.25, initialitem=NULL, ngrid=1000,
               max_iter=200, threshold=0.0000001,contrast_m=NULL, eq_constraint=NULL,
               par_id=NULL, grouping=NULL, est_cov=FALSE, f_cov=NULL, estimation="EM"){
  # data to matrix
  data <- as.matrix(data)
  if(is.null(colnames(data))) colnames(data) <- paste0("v", 1:ncol(data))


  # setting grid on the latent space
  x <- seq(range[1], range[2], length.out=q)
  grid_list <- replicate(dimension, x, simplify = FALSE)
  grid <- as.matrix(do.call(expand.grid, grid_list))

  # saving initial values
  init <- initialitem

  Options = list(initialitem=initialitem,dimension=dimension,data=data,range=range,
                 q=q, max_iter=max_iter, threshold=threshold)

  # preparation for EM
  I <- initialitem
  if(is.null(f_cov)) f_cov <- diag(dimension)
  sds <- sqrt(diag(f_cov))
  prior <- mvtnorm::dmvnorm(grid,
                            mean = rep(0, dimension),
                            sigma = f_cov
  )
  prior <- prior/sum(prior)
  iter <- 0
  diff <- 1

  # preparation for MHRM
  if(estimation == "MHRM"){
    n_sample <- 3 # the number of samples in the stochastic imputation step

    thetalist <- list()
    theta <- matrix(0, nrow = nrow(data), ncol = dimension)
    data_array <- one_hot_3d(data)
    g_approx <- list(0,0)
    H_approx <- list(0,0)

    proposal_sd <- 2.38 # kernel for MH sampling
  }

  # prior for the thresholds
  if(length(t_prior) == 1) t_prior <- rep(t_prior, nrow(initialitem[[1]]))

  # EM
  repeat{
    iter <- iter + 1

    cut_score <- apply(initialitem[[2]][,-1], 1, cut_trans, simplify = FALSE)
    #--------------------------------
    # EM algorithm
    #--------------------------------
    if(estimation == "EM"){
      E <- Estep_cpp(data,
                     cbind(initialitem[[1]],initialitem[[2]][,1]),
                     grid,
                     prior,
                     cut_score)
      dim(E$e.response) <- c(nrow(initialitem[[1]]), nrow(grid), max(data, na.rm = TRUE) + 1)

      M <- Mstep_sim(E,
                     initialitem,
                     par_id,
                     grouping,
                     prior = if(iter<11) t_prior/10 else t_prior,
                     ngrid = if(iter<11) 100 else ngrid)

      initialitem <- M$item
      grad_norm <- M$grad_norm
      eff_par <- M$eff_par
      logL <- E$logL

      # estimating cov matrix
      factor_means <- as.vector(E$Ak%*%E$grid)
      cov_mat <- t(E$grid) %*% sweep(E$grid, 1, E$Ak, FUN = "*") - factor_means %*% t(factor_means)
      sds <- sqrt(diag(cov_mat))
      cov_mat <- cov_mat / (sds %o% sds)

      # update prior
      if(est_cov){
        f_cov <- cov_mat
        prior <- mvtnorm::dmvnorm(grid,
                                  mean = rep(0, dimension),
                                  sigma = f_cov)
        prior <- prior/sum(prior)
      }
      #--------------------------------
      # MHRM algorithm
      #--------------------------------
    } else if (estimation =="MHRM"){
      logL <- 0
      AR <- 0 # acceptance rate
      for(burnin in 1:n_sample){
        thetaUpdate <- sample_mhrm(data, theta, cbind(initialitem[[1]],initialitem[[2]][,1]), cut_score, f_cov, proposal_sd)
        theta <- thetalist[[burnin]] <- thetaUpdate$theta
        logL <- logL + thetaUpdate$logL / n_sample
        AR <- AR + thetaUpdate$AR / n_sample
      }
      if(AR < .3) proposal_sd <- proposal_sd * 0.8
      if(AR > .5) proposal_sd <- proposal_sd * 1.25

      updated_gH <- update_gH(data_array = data_array,
                              thetalist = thetalist,
                              item = initialitem,
                              par_id = par_id,
                              grouping = grouping,
                              ngrid = if(iter<4) 100 else ngrid,
                              prev_g = g_approx,
                              prev_H = H_approx,
                              prior = if(iter<4) t_prior/10 else t_prior,
                              stepsize = if(iter <= 20) 1 else 1/(iter-20),
                              f_cov = f_cov)
      initialitem <- updated_gH$item

      g_approx <- updated_gH$g_approx
      H_approx <- updated_gH$H_approx
      grad_norm <- updated_gH$grad_norm
      if(est_cov){
        f_cov <- updated_gH$f_cov
        sds <- sqrt(diag(f_cov))
        f_cov <- f_cov / (sds %o% sds)
      }
      eff_par <- updated_gH$eff_par
      factor_means <- rep(0, dimension) # updated_gH$factor_means
    }


    # adjust parameters
    if(dimension == 1){
      initialitem[[1]][,1] <- initialitem[[1]][,1] * sds
      initialitem[[1]][,2] <- initialitem[[1]][,2] + initialitem[[1]][,1] * factor_means
    }else{
      initialitem[[1]][,1:dimension] <- sweep(initialitem[[1]][,1:dimension], 2, sds, FUN = "*")
      initialitem[[1]][,dimension+1] <- initialitem[[1]][,dimension+1] + initialitem[[1]][,1:dimension] %*% factor_means
    }
    initialitem[[1]][is.na(par_id[,1:(dimension+1)])] <- init[[1]][is.na(par_id[,1:(dimension+1)])]
    if(!is.null(eq_constraint) & (max(eq_constraint)>0)){
      for(i in 1:max(eq_constraint)){
        initialitem[[1]][which(eq_constraint[,1:(dimension+1)]==i)] <- mean(initialitem[[1]][which(eq_constraint[,1:(dimension+1)]==i)])
        initialitem[[2]][which(eq_constraint[,(dimension+2):ncol(eq_constraint)]==i)] <- mean(initialitem[[2]][which(eq_constraint[,(dimension+2):ncol(eq_constraint)]==i)])
      }
    }

    diff <- c(
      max(abs(I[[1]]-initialitem[[1]]), na.rm = TRUE),
      max(abs(I[[2]]-initialitem[[2]]), na.rm = TRUE)
      )
    max_par <- par_id[which.max(cbind(abs(I[[1]]-initialitem[[1]]),abs(I[[2]]-initialitem[[2]])))]
    diff <- max(diff)
    grad_conv <- grad_norm/(abs(logL)+1e-06)
    if(is.na(grad_conv)) grad_conv <- 1
    I <- initialitem
    if(estimation == "EM"){
      message("\r","\r","EM cycle = ",iter,", logL = ", sprintf("%.2f", logL), ", GradConv = ", sprintf("%.6f",grad_conv),", Max-Change = ", sprintf("%.5f", diff), " at (", paste(max_par, collapse = ","), ")  ", sep="",appendLF=FALSE)
      flush.console()
      if((iter > 11) & (iter >= max_iter | grad_conv < threshold) & (diff < 0.01)) break
    } else if(estimation == "MHRM"){
      message("\r","\r","Phase", if(iter <= 20) 1 else 2, ", MHRM cycle = ",iter, ", (AR|",sprintf("%.1f", AR),"), logL = ", sprintf("%.2f", logL),", Max-Change = ", sprintf("%.5f", diff), " at (", paste(max_par, collapse = ","), ")  ", sep="",appendLF=FALSE)
      flush.console()
      if((iter > 3) & (iter >= max_iter | diff < threshold)) break
    }
  }

  # preparation for outputs
  dimnames(initialitem[[1]]) <- list(colnames(data),c(paste0("f",1:dimension), "c"))
  dimnames(initialitem[[2]]) <- list(colnames(data),c("nu", paste("t", 1:(ncol(initialitem[[2]])-1), sep="")))
  initialitem[[2]][,1] <- exp(initialitem[[2]][,1])


  if (estimation =="EM"){
    theta <- E$posterior%*%E$grid
    theta_se <- sqrt(E$posterior%*%(E$grid^2)-theta^2)

    M <- Mstep_sim(E,
                   I,
                   par_id,
                   grouping,
                   prior = t_prior,
                   ngrid = ngrid,
                   calculate_m = TRUE)

    return(structure(
      list(par_est=initialitem,
           se=M,
           fk=E$freq,
           iter=iter,
           quad=grid,
           diff=diff,
           prior=E$prior,
           posterior=E$posterior,
           Ak=E$Ak,
           theta = theta,
           theta_se = theta_se,
           logL= logL,
           AIC= -2*logL+2*eff_par,
           BIC= -2*logL+log(mean(colSums(!is.na(data))))*eff_par,
           l_prior_w=M$l_prior_w,
           f_means=factor_means,
           cov_mat = f_cov,
           eff_par = eff_par,
           Options = Options # specified argument values
      ),
      class = c("sg", "dr", "list")
    )
    )
  } else if (estimation =="MHRM"){
    return(structure(
      list(par_est=initialitem,
           iter=iter,
           logL= logL,
           f_means=factor_means,
           cov_mat = f_cov,
           eff_par = eff_par,
           Options = Options # specified argument values
      ),
      class = c("sg", "dr", "list")
    )
    )
  }
}

dr_group <- function(
    data, group, dimension=NULL, range=c(-4,4), q=11, t_prior=0.25, initialitem=NULL, ngrid=1000,
    max_iter=200, threshold=0.001, contrast_m=NULL, eq_constraint=NULL,
    par_id=NULL, grouping=NULL, est_cov=FALSE, f_cov=NULL){

  # data to matrix
  data <- as.matrix(data)
  if(is.null(colnames(data))) colnames(data) <- paste0("v", 1:ncol(data))

  # the number of groups
  ngroup <- length(unique(group))
  nitem_g <- nrow(initialitem[[1]])/ngroup

  # setting grid on the latent space
  x <- seq(range[1], range[2], length.out=q)
  grid_list <- replicate(dimension, x, simplify = FALSE)
  grid <- as.matrix(do.call(expand.grid, grid_list))

  # saving initial values
  init <- initialitem

  Options = list(initialitem=initialitem,dimension=dimension,data=data,range=range,
                 q=q, max_iter=max_iter, threshold=threshold, par_id=par_id, group=group)

  # preparation for EM
  I <- initialitem
  prior0 <- mvtnorm::dmvnorm(grid,
                             mean = rep(0, dimension),
                             sigma = if(is.null(f_cov)) diag(dimension) else f_cov
  )
  prior0 <- prior0/sum(prior0)
  prior <- list()
  for(g in 1:ngroup){
    prior[[g]] <- prior0
  }

  iter <- 0
  diff <- 1

  if(length(t_prior) == 1) t_prior <- rep(t_prior, nrow(initialitem[[1]]))

  # EM
  repeat{
    iter <- iter + 1

    # E-step
    Ak <- list()
    freq <- list()
    posterior <- list()
    logL <- 0
    e.response <- c()
    for(g in 1:ngroup){
      item_index <- ((g-1)*nitem_g + 1):(g*nitem_g)
      cut_score <- apply(initialitem[[2]][item_index,-1], 1, cut_trans, simplify = FALSE)

      E <- Estep_cpp(data[group==g,],
                      cbind(initialitem[[1]][item_index,],initialitem[[2]][item_index,1]),
                      grid,
                      prior[[g]],
                      cut_score)

      dim(E$e.response) <- c(nitem_g, nrow(grid), max(data, na.rm = TRUE) + 1)

      e.response <- abind::abind(e.response, E$e.response, along = 1)

      Ak[[g]] <- E$Ak
      freq[[g]] <- E$freq
      posterior[[g]] <- E$posterior
      logL <- logL + E$logL
    }

    # M-step
    M <- Mstep_sim(list(grid = grid, e.response = e.response),
                   initialitem,
                   par_id,
                   grouping,
                   prior = if(iter<4) t_prior/10 else t_prior,
                   ngrid = if(iter<4) 100 else ngrid)
    initialitem <- M$item

    # estimating cov matrix
    factor_means <- list()
    cov_mat <- list()
    sds <- list()
    for(g in 1:ngroup){
      factor_means[[g]] <- as.vector(Ak[[g]] %*% grid)
      cov_mat[[g]] <- t(grid) %*% sweep(grid, 1, Ak[[g]], FUN = "*") - factor_means[[g]] %*% t(factor_means[[g]])
      sds[[g]] <- sqrt(diag(cov_mat[[g]]))
    }


    # update prior
    if(est_cov){
      f_cov <- list()
      prior <- list()
      for(g in 1:ngroup){
        f_cov[[g]] <- cov_mat[[g]]
        prior[[g]] <- mvtnorm::dmvnorm(grid,
                                   mean = factor_means[[g]] - factor_means[[1]],
                                   sigma = f_cov[[g]] / (sds[[1]] %o% sds[[1]]))
        prior[[g]] <- prior[[g]]/sum(prior[[g]])
      }
    }

    # adjust parameters
    if(dimension == 1){
      initialitem[[1]][,1] <- initialitem[[1]][,1] * sds[[1]]
      initialitem[[1]][,2] <- initialitem[[1]][,2] + initialitem[[1]][,1] * factor_means[[1]]
    }else{
      initialitem[[1]][,1:dimension] <- sweep(initialitem[[1]][,1:dimension], 2, sds[[1]], FUN = "*")
      initialitem[[1]][,dimension+1] <- initialitem[[1]][,dimension+1] + initialitem[[1]][,1:dimension] %*% factor_means[[1]]
    }
    initialitem[[1]][is.na(par_id[,1:(dimension+1)])] <- init[[1]][is.na(par_id[,1:(dimension+1)])]
    if(!is.null(eq_constraint) & (max(eq_constraint)>0)){
      for(i in 1:max(eq_constraint)){
        initialitem[[1]][which(eq_constraint[,1:(dimension+1)]==i)] <- mean(initialitem[[1]][which(eq_constraint[,1:(dimension+1)]==i)])
        initialitem[[2]][which(eq_constraint[,(dimension+2):ncol(eq_constraint)]==i)] <- mean(initialitem[[2]][which(eq_constraint[,(dimension+2):ncol(eq_constraint)]==i)])
      }
    }

    diff <- c(
      max(abs(I[[1]]-initialitem[[1]]), na.rm = TRUE),
      max(abs(I[[2]]-initialitem[[2]]), na.rm = TRUE)
    )

    max_par <- par_id[which.max(cbind(abs(I[[1]]-initialitem[[1]]),abs(I[[2]]-initialitem[[2]])))]
    diff <- max(diff)
    grad_conv <- M$grad_norm/(abs(logL)+1e-06)
    if(is.na(grad_conv)) grad_conv <- 1
    I <- initialitem
    message("\r","\r","EM cycle = ",iter,", logL = ", sprintf("%.2f", logL), ", GradConv = ", sprintf("%.6f",grad_conv),", Max-Change = ", sprintf("%.5f", diff), " at (", paste(max_par, collapse = ","), ")  ", sep="",appendLF=FALSE)
    flush.console()
    if((iter > 3) & (iter >= max_iter | grad_conv < threshold)) break
  }

  # preparation for outputs
  dimnames(initialitem[[1]]) <- list(rownames(initialitem[[1]]),c(paste0("f",1:dimension), "c"))
  dimnames(initialitem[[2]]) <- list(rownames(initialitem[[1]]),c("nu", paste("t", 1:(ncol(initialitem[[2]])-1), sep="")))
  initialitem[[2]][,1] <- exp(initialitem[[2]][,1])

  # theta <- rbind(E1$posterior%*%E1$grid, E2$posterior%*%E2$grid)
  # theta_se <- sqrt(rbind(E1$posterior%*%(E1$grid^2), E2$posterior%*%(E2$grid^2))-theta^2)

  return(structure(
    list(par_est=initialitem,
         # se=M[[2]],
         fk=freq,
         iter=iter,
         quad=grid,
         diff=diff,
         prior=prior,
         posterior=posterior,
         Ak=Ak,
         # theta = theta,
         # theta_se = theta_se,
         logL= logL,
         cov_mat = f_cov,
         f_mean = factor_means,
         eff_par = M$eff_par,
         Options = Options # specified argument values
    ),
    class = c("mg", "dr", "list")
  )
  )
}

################################################################################
# EFA, CFA, & Multiple Group Analyses
################################################################################
efa_sim <- function(data,dimension,t_prior=0.25,range=c(-4,4),q=41,ngrid=1000,max_iter=NULL,
                    threshold=NULL, eq_interval = FALSE, estimation="EM"){
  if(is.null(threshold)){
    if(estimation == "EM"){
      threshold <- 0.000001
    }else if(estimation == "MHRM"){
      threshold <- 0.001
    }
  }
  if(is.null(max_iter)){
    if(estimation == "EM"){
      max_iter <- 200
    }else if(estimation == "MHRM"){
      max_iter <- 400
    }
  }

  data <- reorder_mat(as.matrix(data))

  helper_matrices <- efa_helper_matrices(dimension, data, eq_interval = eq_interval)

  fit <- dr_sim(data,
                dimension = dimension,
                contrast_m = helper_matrices$load_mat,
                initialitem = list(
                  helper_matrices$init_mat[,(1:(dimension+1))],
                  helper_matrices$init_mat[,-(1:(dimension+1))]
                ),
                eq_constraint = helper_matrices$eq_constraint,
                par_id = helper_matrices$par_id,
                grouping = helper_matrices$grouping_index,
                ngrid = ngrid,
                max_iter = max_iter,
                threshold = threshold,
                range = range,
                q = q,
                est_cov = FALSE,
                t_prior = t_prior,
                estimation=estimation
  )

  return(fit)
}

cfa_sim <- function(formula,data,t_prior=0.25,range=c(-4,4),q=41,ngrid=1000,max_iter=NULL,
                    threshold=NULL, eq_interval = FALSE, estimation="EM"){
  if(is.null(threshold)){
    if(estimation == "EM"){
      threshold <- 0.000001
    }else if(estimation == "MHRM"){
      threshold <- 0.001
    }
  }
  if(is.null(max_iter)){
    if(estimation == "EM"){
      max_iter <- 500
    }else if(estimation == "MHRM"){
      max_iter <- 1000
    }
  }


  data <- reorder_mat(as.matrix(data))

  helper_matrices <- cfa_helper_matrices(formula, data, eq_interval = eq_interval)

  fit <- dr_sim(data,
                dimension = helper_matrices$d,
                contrast_m = helper_matrices$load_mat,
                initialitem = list(
                  helper_matrices$init_mat[,(1:(helper_matrices$d+1))],
                  helper_matrices$init_mat[,-(1:(helper_matrices$d+1))]
                ),
                eq_constraint = helper_matrices$eq_constraint,
                par_id = helper_matrices$par_id,
                grouping = helper_matrices$grouping_index,
                ngrid = ngrid,
                max_iter = max_iter,
                threshold = threshold,
                range = range,
                q = q,
                est_cov = TRUE,
                t_prior = t_prior,
                estimation=estimation
  )

  return(fit)
}

cfa_group <- function(formula,data,group,t_prior=c(0,0.25),range=c(-4,4),q=41,ngrid=1000,max_iter=200,
                      threshold=0.00001, eq_interval = FALSE){

  data <- reorder_mat(as.matrix(data))

  helper_matrices <- group_helper_matrices(formula, data, eq_interval = eq_interval, ngroup = length(unique(group)))

  fit <- dr_group(data,
                  group = group,
                  dimension = helper_matrices$d,
                  contrast_m = helper_matrices$load_mat,
                  initialitem = list(
                    helper_matrices$init_mat[,(1:(helper_matrices$d+1))],
                    helper_matrices$init_mat[,-(1:(helper_matrices$d+1))]
                  ),
                  eq_constraint = helper_matrices$eq_constraint,
                  par_id = helper_matrices$par_id,
                  grouping = helper_matrices$grouping_index,
                  ngrid = ngrid,
                  max_iter = max_iter,
                  threshold = threshold,
                  range = range,
                  q = q,
                  est_cov = TRUE,
                  t_prior = t_prior
  )

  return(fit)
}

################################################################################
# PLOTTING
################################################################################
library(ggplot2)

# function
plot_DRM <- function(thresholds, item_name = NULL, legend=FALSE){
  thresholds <- c(thresholds, 1)
  thresholds <- c(thresholds[1], diff(thresholds))
  data <- data.frame(
    Category = "a",  # Two bars
    Subgroup = paste0("Cat",1:length(thresholds)),         # Stacked components
    Value = rev(thresholds)               # Values to be proportionally scaled
  )
  ppp <- ggplot(data, aes(x = Category, y = Value, fill = Subgroup)) +
    geom_bar(stat = "identity", position = "fill", show.legend = legend) +  # Normalize stacking to sum to 1
    scale_fill_brewer(palette = "Set2") +
    # geom_point(data = thresholds, aes(x = Category, y = Threshold, color = Color),size = 3, inherit.aes = FALSE) +  # Add points at threshold
    scale_color_identity() +
    scale_x_discrete(expand = c(0, 0)) +
    # coord_cartesian(ylim = c(0, 1)) +
    coord_flip() +  # Flip to horizontal
    scale_y_continuous(
      labels = scales::percent_format(scale = 1),
      breaks = seq(0, 1, by = 1/length(thresholds))  # Set y-axis ticks at 0%, 20%, 40%, etc.
    ) +
    labs(x = item_name, y = NULL, fill = "Category") +
    theme(plot.margin = unit(c(0, 0.1, 0, 0.1),
                             "inches"),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_text(angle = 0),
          panel.background = element_blank(),  # Remove panel background
          plot.background = element_blank(),   # Remove plot background
          panel.grid = element_blank() )

  return(ppp)
}

plot_DRM_compare <- function(..., label = NULL, item_name = NULL, legend = FALSE) {
  # Collect threshold vectors
  t_list <- list(...)
  n_cat <- length(t_list)

  # Default labels if not provided
  if (is.null(label)) label <- paste0("Group", seq_len(n_cat))
  if (length(label) != n_cat) stop("Length of `label` must match number of threshold vectors.")

  # Convert thresholds → interval lengths
  process_threshold <- function(t) {
    t <- c(t, 1)
    t <- c(t[1], diff(t))
    return(rev(t))
  }

  t_processed <- lapply(t_list, process_threshold)

  # Build data for ggplot
  data <- do.call(rbind, lapply(seq_along(t_processed), function(i) {
    data.frame(
      Category = label[i],
      Subgroup = paste0("Cat", seq_along(t_processed[[i]])),
      Value = t_processed[[i]]
    )
  }))

  p <- ggplot(data, aes(x = Category, y = Value, fill = Subgroup)) +
    geom_bar(stat = "identity", position = "fill", show.legend = legend) +
    scale_fill_brewer(palette = "Set2") +
    coord_flip() +
    scale_y_continuous(
      labels = scales::percent_format(),
      breaks = seq(0, 1, length.out = length(t_processed[[1]]) + 1)
    ) +
    labs(x = item_name, y = NULL, fill = "Subgroup") +
    theme(
      plot.margin = unit(c(0, .1, 0, .1), "inches"),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(angle = 0),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid = element_blank()
    )

  return(p)
}
