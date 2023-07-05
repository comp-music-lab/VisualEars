inv_logit <- function(x){
  1/(1+exp(-x))
}

clmm_new_response <- function(participant, levs, theta_eff){
  # print(participant)
  cdf <- inv_logit(theta_eff[, participant["group"] + 1] -
                     participant["random_effect"])
  # print(cdf)
  pmf <- c(cdf, 1) - c(0, cdf)
  # print(pmf)
  return(sample(x = levs,
                size = participant["n_trials"],
                replace = TRUE,
                prob = pmf))
}

clmm_test <- function(random_effects, participants, levs, theta_eff) {
  # generate new responses
  n_participants <- nrow(participants)
  participants <- cbind(participants, random_effects)
  colnames(participants)[4] <- "random_effect"
  new_pas <- apply(participants,
                   1, 
                   clmm_new_response, 
                   levs = levs, 
                   theta_eff = theta_eff)
  
  # preprocess variables
  new_pas <- as.factor(c(new_pas))
  participant_id <- as.factor(rep(participants[, "id"], participants[, "n_trials"]))
  participant_group <- rep(participants[, "group"], participants[, "n_trials"])
  
  # fit CLMM
  sim_fit <- ordinal::clmm(
    formula = new_pas ~ participant_group + (1|participant_id),
    control = ordinal::clmm.control(innerCtrl = "noWarn"))
  
  # converged? if so, reject or fail to reject?
  p_val <- coef(summary(sim_fit))['participant_group', 'Pr(>|z|)']
  if (!is.na(p_val)) {
    if (p_val < .05) {
      return(1)
    }
  }
  return(0)
}

mann_whitney_clmm_test <- function(random_effects, participants, levs, theta_eff) {
  n_participants <- nrow(participants)
  participants <- cbind(participants, random_effects)
  colnames(participants)[4] <- "random_effect"
  new_pas <- apply(participants,
                   1, 
                   clmm_new_response, 
                   levs = levs, 
                   theta_eff = theta_eff)
  max_pas <- as.numeric(as.character(apply(new_pas, 2, max)))
  test <- wilcox.test(x = max_pas[participants[, "group"]==0], 
                      y = max_pas[participants[, "group"]==1], 
                      alternative = "less",
                      mu = 0,
                      paired = F,
                      exact = F)
  if (test$p.value < .05) {
    return(1) 
  }
  return(0)
}

wilcoxon_sr_clmm_test <- function(random_effects, participants, levs, theta_eff) {
  n_participants <- nrow(participants)
  participants <- dplyr::inner_join(
    as.data.frame(participants), 
    data.frame(id = 1:n_participants, random_effects),
    by = "id")
  colnames(participants)[4] <- "random_effect"
  new_pas <- apply(participants,
                   1, 
                   clmm_new_response, 
                   levs = levs, 
                   theta_eff = theta_eff)
  max_pas <- as.numeric(as.character(apply(new_pas, 2, max)))
  test <- wilcox.test(x = max_pas[participants[, "group"]==0], 
                      y = max_pas[participants[, "group"]==1], 
                      alternative = "less",
                      mu = 0,
                      paired = T,
                      exact = F)
  if (test$p.value < .05) {
    return(1)
  }
  return(0)
}

clmm_power <- function(effect, parameters, participants, levs, n_sim = 10^3, test = "clmm", verbose = F) {
  thresholds <- parameters[[1]]
  sigma <- parameters[[2]]
  if (length(levs) != length(thresholds) + 1) {
    stop("The number of levels is not consistent with the thresholds provided.")
  }
  if (verbose) cat(" effect =", effect, ",")
  
  # simulate matrix of random effects; one experiment per row
  n_participants <- nrow(participants)
  random_effects <- matrix(rnorm(n = n_participants * n_sim,
                                 mean = 0, sd = sigma),
                           nrow = n_sim)
  random_effects <- split(random_effects,
                          rep(1:nrow(random_effects), ncol(random_effects)))
  
  # construct interaction between thresholds and effect
  theta_eff <- cbind(thresholds, thresholds - effect)
  
  # run test n_sim times and save results
  if (test == "clmm") f <- "clmm_test"
  if (test == "Mann-Whitney U") f <- "mann_whitney_clmm_test"
  if (test == "Wilcoxon SR") f <- "wilcoxon_sr_clmm_test"
  test_res <- unlist(parallel::mclapply(
    # test_res <- unlist(lapply(
    X = random_effects,
    FUN = f,
    participants,
    levs,
    # theta_eff))
    theta_eff,
    mc.cores = 1))
  if (verbose) cat(" power =", mean(test_res), "\n")
  return(mean(test_res))
}

clmm_power_curve <- function(effects, parameters, participants, levs, n_sim = 10^3, test = "clmm", verbose = F) {
  thresholds <- parameters[[1]]
  if (length(levs) != length(thresholds) + 1) {
    stop("The number of levels is not consistent with the thresholds provided.")
  }
  if (verbose) {
    cat("CLMM thresholds:", thresholds, "\n")
    sigma <- parameters[[2]]
    cat("SD random effect:", sigma, "\n")
  }
  sapply(effects, clmm_power, parameters, participants, levs, n_sim, test, verbose)
}

clmm_power_band <- function(effects, parameters, data, control_group, n_sim = 10^3, test = "clmm", verbose = F) {
  participants <- as.matrix(
    data %>% 
      dplyr::group_by(id, group) %>% 
      dplyr::summarise(group = ifelse(group[1] == control_group, 0, 1),
                       n_trials = n())
  )
  levs <- sort(unique(data$pas))
  thresholds <- parameters$thetas
  if (!is.matrix(thresholds)) thresholds <- as.matrix(t(thresholds))
  K <- ncol(thresholds) + 1
  sigmas <- parameters$sigmas
  n_curves <- length(sigmas)
  power_band <- matrix(NA, n_curves * length(effects))
  for (i in 1:n_curves) {
    if (verbose) tic <- Sys.time()
    parameters_curve <- list(thresholds[i, ], sigmas[i])
    power_band[(i-1) * length(effects) + (1:length(effects))] <-
      clmm_power_curve(effects, parameters_curve, participants, levs, n_sim, test, verbose)
    if (verbose) print(Sys.time() - tic)
  }
  times <- rep(length(effects), n_curves)
  power_band <- data.frame(curve = rep(1:n_curves, times),
                           distr_idx = rep(parameters$distr_idx, times),
                           sigma = rep(sigmas, times),
                           effect = rep(effects, n_curves),
                           power = power_band)
  return(power_band)
}

extract_effects <- function(power_band, power) {
  power_minmax <- power_band %>% 
    dplyr::group_by(effect) %>% 
    dplyr::summarise(min = min(power, na.rm = T), max = max(power, na.rm = T))
  cols <- c("max", "min")
  odds_ratios <- rep(NA, 2)
  for (i in 1:2) {
    i_lower <- max(which(power_minmax[, cols[i]] < power))
    if (is.infinite(i_lower)) stop("Power is too small.")
    if (i_lower == nrow(power_minmax))  stop("Power is too large.")
    m <- (power_minmax[i_lower + 1, cols[i]] - power_minmax[i_lower, cols[i]]) / (exp(power_minmax$effect[i_lower+1]) - exp(power_minmax$effect[i_lower]))
    q <- power_minmax[i_lower, cols[i]] - m * exp(power_minmax$effect[i_lower])
    odds_ratios[i] <- as.numeric((power - q) / m)
  }
  return(odds_ratios)
}

clmm_generate_data <- function(n_participants, 
                               n_trials,
                               control_distribution,
                               effect, 
                               participant_variation,
                               within_subject = F,
                               control_weight = .5) {
  participants <- data.frame(
    id = 1:n_participants,
    n_trials = rep(n_trials, n_participants),
    random_effect = rnorm(n_participants, 0, participant_variation))
  thresholds <- clmm_recover_thresholds(control_distribution, participant_variation)$thetas
  theta_eff <- cbind(t(thresholds), t(thresholds) - effect)
  if (!within_subject) {
    n_control <- round(control_weight * n_participants)
    participants$group <- c(rep(0:1,  c(n_control, n_participants - n_control)))
    data <- data.frame(
      id = rep(participants$id, rep(n_trials, n_participants)), 
      pas = c(apply(participants,
                    1, 
                    clmm_new_response, 
                    levs = 1:length(control_distribution), 
                    theta_eff = theta_eff)))
  } else {
    participants_ws <- rbind(participants, participants)
    participants_ws$group <- c(rep(0:1,  rep(n_participants, 2)))
    data <- data.frame(
      id = rep(participants_ws$id, rep(n_trials, length(participants_ws$id))),
      group = rep(participants_ws$group, rep(n_trials, length(participants_ws$id))),
      pas = c(apply(participants_ws,
                    1, 
                    clmm_new_response, 
                    levs = 1:length(control_distribution), 
                    theta_eff = theta_eff)))
  }
  data <- dplyr::inner_join(data, 
                            participants[, colnames(participants)!="n_trials"],
                            by = "id")
  return(data)
}

clmm_marginal_cdf <- function(theta, sigma, n_sim) {
  K <- length(theta) + 1
  x <- rnorm(n_sim, sd = sigma)
  cdf <- rep(NA, K)
  for (k in 1:(K-1)) {
    cdf[k] <- mean(inv_logit(theta[k] - x))
  }
  cdf[K] <- 1
  return(cdf - c(0, cdf[-K]))
}

clmm_marginal_cdf_wrapper <- function(theta, p, sigma, n_sim) {
  sum(abs(p - clmm_marginal_cdf(theta, sigma, n_sim = n_sim)))
}

clmm_inv_theta<- function(p, sigma, ini, n_sim = 10^5) {
  optim(par = ini, fn = clmm_marginal_cdf_wrapper, p = p, sigma = sigma, n_sim = n_sim)
}

clmm_inv_theta_wrapper = function(target, verbose = F) { # wrapper for parallel::mclapply
  clmm_inv_theta(p = target$distr, 
                 sigma = target$sigma,
                 ini = target$ini[-length(target$distr)])
}

clm_inv_theta <- function(p) {
  cdf <- cumsum(p)
  cdf[length(cdf)] <- 1
  thetas <- -log(1/cdf - 1)
  return(thetas)
}

clm_recover_thresholds <- function(p) {
  cdf <- cumsum(p)
  cdf[length(cdf)] <- 1
  thetas <- -log(1/cdf - 1)
  return(thetas)
}

clmm_recover_thresholds <- function(distrs, sigmas, verbose = F) {
  if (!is.matrix(distrs)) distrs <- as.matrix(t(distrs))
  K <- ncol(distrs)
  targets <- expand.grid(1:nrow(distrs), sigmas)
  targets <- cbind(1:nrow(targets), targets)
  names(targets) <- c("curve", "distr", "sigma")
  
  thresholds <- matrix(NA, nrow = nrow(targets), ncol = K-1)
  
  # initial values for the optimizer
  theta_ini <- t(apply(distrs, 1, clm_recover_thresholds))
  
  # compute thresholds
  targets_list <- vector("list", length = nrow(targets))
  for (i in 1:nrow(targets)){
    targets_list[[i]]$distr <- distrs[targets[i, "distr"],]
    targets_list[[i]]$sigma <- targets[i, "sigma"]
    targets_list[[i]]$ini <- theta_ini[targets[i, "distr"],]
  }
  out <- parallel::mclapply(targets_list,
                            clmm_inv_theta_wrapper,
                            verbose = verbose,
                            mc.cores = 1)
  thresholds <-
    matrix(unlist(lapply(out, function(x) x$par)),
           byrow = T,
           ncol = K-1)
  
  # check that each set of recovered thresholds form an increasing sequence
  sorted_thresholds <- t(apply(thresholds, 1, sort))
  diff_thresholds <- thresholds - sorted_thresholds
  not_sorted <- which(apply(diff_thresholds, 1, function(x) any(x>0)))
  if (length(not_sorted) > 0) {
    for (i in not_sorted) {
      warning(paste0("The recovered thresholds for distr = (", paste(targets_list[[i]]$distr, collapse = ", "),
                     ") and sigma = ", paste(targets_list[[i]]$sigma, collapse = " "), " did not form an increasing sequence. We sorted the thresholds before returning them."))
    }
    thresholds <- sorted_thresholds
  }
  parameters <- list(thetas = thresholds,
                     distr = distrs[targets[, "distr"], ],
                     distr_idx = targets[, "distr"],
                     sigmas = targets[, "sigma"])
  return(parameters)
}