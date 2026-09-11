hima2 <- function(formula, 
                  data.pheno, 
                  data.M,  
                  outcome.family = c("gaussian", "binomial", "survival", "quantile"), 
                  mediator.family = c("gaussian", "negbin", "compositional"), 
                  penalty = c("DBlasso", "MCP", "SCAD", "lasso"), 
                  topN = NULL, 
                  scale = TRUE,
                  verbose = FALSE,
                  ...) 
{
  outcome.family <- match.arg(outcome.family)
  mediator.family <- match.arg(mediator.family)
  penalty <- match.arg(penalty)
  
  # Penalty check
  if (penalty == "DBlasso" & outcome.family == "quantile")
  {
    message("Note: Quantile HIMA does not support De-biased Lasso penalty. Switing to 'MCP' (default) ...")
    penalty = "MCP"
  }
  
  if (penalty != "DBlasso" & (outcome.family == "survival" | mediator.family == "compositional"))
  {
    message("Note: Survival HIMA and Compositional HIMA can be only performed using De-biased Lasso. Switing to 'DBlasso' ...")
    penalty = "DBlasso"
  }
  
  # DBlasso
  if (penalty == "DBlasso")
  {
    if (outcome.family %in% c("gaussian", "binomial"))
    {
      if(mediator.family %in% c("gaussian", "negbin"))
      {
        response_var <- as.character(formula[[2]]) 
        ind_vars <- all.vars(formula)[-1]
        
        Y <- data.pheno[,response_var]
        X <- data.pheno[,ind_vars[1]]
        
        if(length(ind_vars) > 1)
          COV <- data.pheno[,ind_vars[-1]] else COV <- NULL
        
        results <- dblassoHIMA(X = X, Y = Y, M = data.M, Z = COV, 
                               Y.family = outcome.family, 
                               topN = topN,
                               scale = scale, verbose = verbose)
        
        attr(results, "variable.labels") <- c("alpha: Effect of exposure on mediator", 
                                              "beta: Effect of mediator on outcome",
                                              "gamma: Total effect of exposure on outcome",
                                              "alpha*beta: Mediation effect",
                                              "% total effect: Percent of mediation effect out of the total effect",
                                              "p.joint: Joint raw p-value of selected significant mediator (based on FDR)")
      } else if (mediator.family == "compositional") {
        response_var <- as.character(formula[[2]]) 
        ind_vars <- all.vars(formula)[-1]
        
        Y <- data.pheno[,response_var]
        X <- data.pheno[,ind_vars[1]]
        
        if(length(ind_vars) > 1)
          COV <- data.pheno[,ind_vars[-1]] else COV <- NULL
        
        res <- microHIMA(X = X, Y = Y, OTU = data.M, COV = COV, FDRcut = 0.3, scale)
        results <- data.frame(alpha = res$alpha, alpha_se = res$alpha_se, 
                              beta = res$beta, beta_se = res$beta_se,
                              FDR = res$FDR, check.names = FALSE)
        rownames(results) <- res$ID
        attr(results, "variable.labels") <- c("alpha: Effect of exposure on mediator", 
                                              "alpha_se: Standard error of the effect of exposure on mediator",
                                              "beta: Effect of mediator on outcome",
                                              "beta_se: Standard error of the effect of mediator on outcome",
                                              "FDR: Hommel's false discovery rate")
      }
    } else if (outcome.family == "survival") {
      response_vars <- as.character(formula[[2]])[c(2,3)]
      ind_vars <- all.vars(formula)[-c(1,2)]
      
      X <- data.pheno[,ind_vars[1]]
      status <- data.pheno[, response_vars[1]]
      OT <- data.pheno[, response_vars[2]]
      
      if(length(ind_vars) > 1)
        COV <- data.pheno[,ind_vars[-1]] else COV <- NULL
      
      res <- survHIMA(X, COV, data.M, OT, status, FDRcut = 0.3, scale, verbose)
      
      results <- data.frame(alpha = res$alpha, alpha_se = res$alpha_se, 
                            beta = res$beta, beta_se = res$beta_se,
                            p.joint = res$p.joint, check.names = FALSE)
      rownames(results) <- res$ID
      attr(results, "variable.labels") <- c("alpha: Effect of exposure on mediator", 
                                            "alpha_se: Standard error of the effect of exposure on mediator",
                                            "beta: Effect of mediator on outcome",
                                            "beta_se: Standard error of the effect of mediator on outcome",
                                            "p.joint: Joint raw p-value of selected significant mediator (based on FDR)")
    }
  } else { # If penalty is not DBlasso
    if (outcome.family %in% c("gaussian", "binomial"))
    {
      response_var <- as.character(formula[[2]]) 
      ind_vars <- all.vars(formula)[-1]
      
      Y <- data.pheno[,response_var]
      X <- data.pheno[,ind_vars[1]]
      
      if(length(ind_vars) > 1)
        COV <- data.pheno[,ind_vars[-1]] else COV <- NULL
      
      results <- hima(X = X, Y = Y, M = data.M, COV.XM = COV, 
                      Y.family = outcome.family, M.family = mediator.family, 
                      penalty = penalty, topN = topN,
                      parallel = FALSE, ncore = 1, scale = scale, verbose = verbose)
      
      attr(results, "variable.labels") <- c("alpha: Effect of exposure on mediator", 
                                            "beta: Effect of mediator on outcome",
                                            "gamma: Total effect of exposure on outcome",
                                            "alpha*beta: Mediation effect",
                                            "% total effect: Percent of mediation effect out of the total effect",
                                            "Bonferroni.p: Bonferroni adjusted p value",
                                            "BH.FDR: Benjamini-Hochberg False Discovery Rate")
    } else if (outcome.family == "quantile") {
      # tau <- readline(prompt = "Enter quantile level(s) (between 0-1, multiple values accepted): ")
      # tau <- eval(parse(text = paste0("c(", tau, ")")))
      
      response_var <- as.character(formula[[2]]) 
      ind_vars <- all.vars(formula)[-1]
      
      Y <- data.pheno[,response_var]
      X <- data.pheno[,ind_vars[1]]
      
      if(length(ind_vars) > 1)
        COV <- data.pheno[,ind_vars[-1]] else COV <- NULL
      
      res <- qHIMA(X = X, M = data.M, Y = Y, Z = COV,
                   Bonfcut = 0.3, penalty = penalty, scale = scale, verbose = verbose, ...)
      
      results <- data.frame(alpha = res$alpha, alpha_se = res$alpha_se, 
                            beta = res$beta, beta_se = res$beta_se,
                            Bonferroni.p = res$Bonferroni.p, tau = res$tau, 
                            check.names = FALSE)
      rownames(results) <- paste0(res$ID, "-q", res$tau*100) 
      
      attr(results, "variable.labels") <- c("alpha: Effect of exposure on mediator", 
                                            "alpha_se: Standard error of the effect of exposure on mediator",
                                            "beta: Effect of mediator on outcome",
                                            "beta_se: Standard error of the effect of mediator on outcome",
                                            "Bonferroni.p: Bonferroni adjusted p value",
                                            "tau: Quantile level of the outcome")
    }
  }
  return(results)
}