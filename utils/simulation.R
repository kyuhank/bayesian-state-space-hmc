
SimData <- function(draws,
                    SBCtype,
                    ntimes, 
                    Ct,
                    Input,
                    ngen,
                    save_samples_vars = c("K", 
                                          "r", 
                                          "q", 
                                          "sigma2",
                                          "tau2",
                                          "MSY",
                                          "Pt[23]",
                                          "BtoverBmsy[23]",
                                          "HtoverHmsy[23]",
                                          "loglik")) {
  
  # Validate input parameters
  if(ngen > length(draws[[1]])) {
    stop(paste("ngen (", ngen, ") cannot exceed number of available draws (", length(draws[[1]]), ")"))
  }
  
  # Print initial status message
  message(paste("Starting SimData with", ngen, "requested generations and", length(draws[[1]]), "available draws"))
  
  # Initialize storage containers
  variables_list <- list()
  generated_list <- list()
  
  # Initialize counters
  gen_i <- 1          # Current generation index
  draw_idx <- 1       # Current draw index
  rejected_count <- 0 # Count of rejected draws due to biological infeasibility
  
  # Main simulation loop
  while(gen_i <= ngen) {
    
    # Check if all available draws have been exhausted
    if(draw_idx > length(draws[[1]])) {
      warning(paste("Exhausted all available draws before completing", ngen, "generations. Only completed", gen_i-1, "generations out of", ngen, "requested."))
      warning(paste("Total rejected draws:", rejected_count, "out of", length(draws[[1]]), "available draws"))
      break
    }
    
    # Display progress messages at regular intervals
    if(gen_i %% max(1, floor(ngen/10)) == 0 || gen_i == 1) {
      message(paste("Processing generation", gen_i, "of", ngen, "(using draw", draw_idx, ")"))
    }
    
    # Initialize population dynamics vectors
    Pt <- Bt <- ItSim <- Ht <- numeric()
    
    # Extract parameter values from current draw
    draw <- lapply(draws, function(x) x[draw_idx,])
    
    K      <- c(draw$K)
    r      <- c(draw$r)
    q      <- c(draw$q)
    sigma2 <- c(draw$sigma2)
    tau2   <- c(draw$tau2)
    et     <- c(draw$et)
    ut     <- c(draw$ut)
    
    # Calculate standard deviations
    sigma <- sqrt(sigma2)
    tau   <- sqrt(tau2)
    
    # Set initial population state
    Pt[1] <- exp(ut[1] * sigma)
    Bt[1] <- Pt[1] * K
    
    # Calculate biological reference points
    Bmsy <- K / 2
    Hmsy <- r / 2
    msy  <- (r * K) / 4
    
    # Forward simulation of population dynamics
    for(t in 2:(ntimes+1)) {
      Pt[t] <- (Pt[t-1] + r*(1.0 - Pt[t-1])*Pt[t-1] - Ct[t-1]/K) * exp(sigma*ut[t])
      Bt[t] <- Pt[t] * K
      ItSim[t-1] <- q * Bt[t-1] * exp(tau*et[t-1])
      Ht[t-1] <- Ct[t-1] / Bt[t-1]
    }
    
    # Check biological feasibility constraints
    if(all(Bt > 0) & all(Ht >= 0) & all(Ht <= 1)) {
      # Current draw produces biologically feasible trajectory
      
      # Initialize log-likelihood
      loglik <- 0
      ItObs <- Input$It
      
      # Calculate log-likelihood based on SBC type
      for(t in 1:(ntimes)) {
        loglik <- loglik + dnorm(log(ItSim[t]), mean=log(q*Bt[t]), sd=sqrt(tau2), log=TRUE)
        if(SBCtype == "postSBC") {
          loglik <- loglik + dnorm(log(ItObs[t]), mean=log(q*Bt[t]), sd=sqrt(tau2), log=TRUE)
        }
      }
      
      # Initialize variable storage on first successful generation
      if(gen_i == 1) {
        variables_list <- vector("list", length(save_samples_vars))
        names(variables_list) <- save_samples_vars
        for(var_name in save_samples_vars) {
          variables_list[[var_name]] <- numeric(ngen)
        }
      }
      
      # Create mapping of all calculated variables
      all_vars <- list(
        "sigma2" = sigma2,
        "tau2" = tau2,
        "r" = r,
        "K" = K,
        "q" = q,
        "Hmsy" = Hmsy,
        "Bmsy" = Bmsy,
        "MSY" = msy,
        "Pt[20]" = Pt[20],
        "Pt[21]" = Pt[21],
        "Pt[22]" = Pt[22],
        "Pt[23]" = Pt[23],
        "Bt[23]" = Bt[23],
        "Bt[22]" = Bt[22],
        "Bt[21]" = Bt[21],
        "Bt[20]" = Bt[20],
        "BtoverBmsy[23]" = Bt[23]/Bmsy,
        "HtoverHmsy[23]" = Ht[23]/Hmsy,
        "loglik" = loglik
      )
      
      # Store only requested variables
      for(var_name in save_samples_vars) {
        if(var_name %in% names(all_vars)) {
          variables_list[[var_name]][gen_i] <- all_vars[[var_name]]
        } else {
          warning(paste("Variable", var_name, "not found in simulation output"))
        }
      }
      
      # Prepare generated data based on SBC type
      generated_input <- Input
      if(SBCtype == "priorSBC" | SBCtype == "priorSBCdep") {
        generated_input$It <- ItSim
      } else if(SBCtype == "postSBC") {
        generated_input$It2 <- ItSim
      }
      
      # Store generated data for current generation
      generated_list[[gen_i]] <- generated_input
      
      # Move to next generation
      gen_i <- gen_i + 1
      
    } else {
      # Current draw rejected due to biological infeasibility
      rejected_count <- rejected_count + 1
      if(rejected_count %% 50 == 0) {
        message(paste("Warning: Already rejected", rejected_count, "draws due to biological infeasibility"))
      }
    }
    
    # Always advance to next draw regardless of acceptance
    draw_idx <- draw_idx + 1
  }
  
  # Calculate final completion statistics
  actual_gens <- gen_i - 1
  
  if(actual_gens == ngen) {
    # Successfully completed all requested generations
    message(paste("Successfully completed all", ngen, "generations"))
    message(paste("Total draws used:", draw_idx - 1, "out of", length(draws[[1]]), "available"))
    message(paste("Rejection rate:", round(rejected_count/(draw_idx-1)*100, 2), "%"))
  } else {
    # Partially completed due to insufficient draws
    warning(paste("Only completed", actual_gens, "generations out of", ngen, "requested"))
    warning(paste("Consider increasing the number of available draws or relaxing biological constraints"))
  }
  
  # Handle case where no feasible trajectories were found
  if(actual_gens == 0) {
    stop("No biologically feasible trajectories found in the available draws")
  }
  
  # Trim storage containers if fewer generations were completed
  if(actual_gens < ngen) {
    for(var_name in names(variables_list)) {
      variables_list[[var_name]] <- variables_list[[var_name]][1:actual_gens]
    }
    generated_list <- generated_list[1:actual_gens]
  }
  
  # Convert variables list to matrix format
  variables_matrix <- do.call(cbind, variables_list)
  
  # Convert to draws_matrix format required by SBC package
  variables_draws_matrix <- posterior::as_draws_matrix(variables_matrix)
  
  # Create and return SBC_datasets object
  return(SBC:::new_SBC_datasets(variables = variables_draws_matrix,
                                generated = generated_list))
}

SimDataPriorSBC <- function(draws,
                            ntimes,
                            Ct,
                            Input,
                            ngen,
                            save_samples_vars = c("K",
                                                  "r",
                                                  "q",
                                                  "sigma2",
                                                  "tau2",
                                                  "MSY",
                                                  "Pt[23]",
                                                  "BtoverBmsy[23]",
                                                  "HtoverHmsy[23]",
                                                  "loglik")) {
  SimData(
    draws = draws,
    SBCtype = "priorSBC",
    ntimes = ntimes,
    Ct = Ct,
    Input = Input,
    ngen = ngen,
    save_samples_vars = save_samples_vars
  )
}
