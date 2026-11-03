##############################################################################
# Study: 			      KARAT
# Programmer: 		  Lukas Streibl
# Name of program: 	Sample_Size_Calculation_Simulation.R
# Creation date:	  30.10.2026
# Description:      Simulation program for sample size calculation 
##############################################################################

# Set settings ----

# Fixed seed for reproducible simulation result
set.seed(123)

# Load necessary packages
library(survival)
library(ggplot2)
library(gridExtra)

# Assumptions
stage_p <- c(I=.15, II=.35, III=.50)
DFS3 <- c(I=.94, II=.82, III=.72)

mut_p <- .40
DRM_p <- .25

recruitment <- 2.5
study_end <- 5.5

HR_INT_vals <- c(2.844)
HR_M_vals <- c(1)
HR_DRM_vals <- c(1)
dropout_p_vals <- c(.10)

N_vals <- seq(450,465,1)
nsim <- 10000
target_power <- .80

lambda0 <- -log(DFS3)/3

dropout_hazard <- function(p)
  if(p == 0) 0 else -log(1-p)/3

# Study simulation ----

# Simulation of a single study
sim <- function(n, hM, hR, hINT, dropout_p) {
  
  stage <- factor(
    sample(names(stage_p), n, TRUE, stage_p),
    levels = names(stage_p)
  )
  
  M <- rbinom(n, 1, mut_p)
  DRM <- rbinom(n, 1, DRM_p)
  entry <- runif(n, 0, recruitment)
  
  lambda <- lambda0[as.character(stage)] *
    hM^M *
    hR^DRM *
    hINT^(M * DRM)
  
  event <- rexp(n, lambda)
  admin <- study_end - entry
  
  dropout <- if (dropout_p == 0) {
    rep(Inf, n)
  } else {
    rexp(n, dropout_hazard(dropout_p))
  }
  
  time <- pmin(event, admin, dropout)
  status <- as.integer(event <= admin & event <= dropout)
  
  fit <- tryCatch(
    coxph(
      Surv(time, status) ~ M * DRM + strata(stage)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !"M:DRM" %in% names(coef(fit)))
    return(c(p = NA, events = sum(status)))
  
  p <- summary(fit)$coefficients[
    "M:DRM", "Pr(>|z|)"
  ]
  
  c(p = p, events = sum(status))
}

# Sample size calculation ----

# Power for a single sample size
power_N <- function(n, hM, hR, hINT, dropout_p) {
  
  x <- do.call(
    rbind,
    lapply(
      seq_len(nsim),
      function(i)
        sim(n, hM, hR, hINT, dropout_p)
    )
  )
  
  data.frame(
    N = n,
    Power = mean(x[, "p"] < .05, na.rm = TRUE),
    Events = mean(x[, "events"], na.rm = TRUE)
  )
}


# Scenarios
# If more than one set of assumptions is entered above the Power simulation 
# beyond will loop over each scenario separately.
# Sample size is only reproducible for the same set of assumptions
sc <- expand.grid(
  HR_M = HR_M_vals,
  HR_DRM = HR_DRM_vals,
  HR_INT = HR_INT_vals,
  dropout_p = dropout_p_vals
)


# Powersimulation
# Loop over scenarios and test sample sizes 
# Stop per scenario at first sample size with sufficient power 
res_list <- list()
k <- 1

for (i in seq_len(nrow(sc))) {
  
  for (n in N_vals) {
    
    cat(
      "Szenario", i, "von", nrow(sc),
      "- N =", n, "\n"
    )
    
    result <- power_N(
      n = n,
      hM = sc$HR_M[i],
      hR = sc$HR_DRM[i],
      hINT = sc$HR_INT[i],
      dropout_p = sc$dropout_p[i]
    )
    
    result <- cbind(sc[i, ], result)
    
    res_list[[k]] <- result
    k <- k + 1
    
    # Stop if tested sample size has sufficient power (accounting also for Monte-Carlo error)
    if (!is.na(result$Power) &&
        result$Power >= target_power + sqrt(target_power*(1-target_power)/nsim)) {
      break
    }
  }
}

res <- do.call(rbind, res_list)
rownames(res) <- NULL


# Required sample size
req <- do.call(
  rbind,
  lapply(
    split(
      res,
      list(
        res$HR_M,
        res$HR_DRM,
        res$HR_INT,
        res$dropout_p
      ),
      drop = TRUE
    ),
    function(x) {
      
      x <- x[order(x$N), ]
      ok <- x[x$Power >= target_power + sqrt(target_power*(1-target_power)/nsim), ]
      
      data.frame(
        HR_M = x$HR_M[1],
        HR_DRM = x$HR_DRM[1],
        HR_INT = x$HR_INT[1],
        Dropout = 100 * x$dropout_p[1],
        N = if (nrow(ok)) ok$N[1] else NA,
        Power = if (nrow(ok)) ok$Power[1] else NA,
        Events = if (nrow(ok)) ok$Events[1] else NA
      )
    }
  )
)

# Export of results ----

# Create pdf with simulation results
pdf("simulation_results.pdf",6,10)

grid.table(
  round(res[,c(1:7)],3),
  rows=NULL,
  theme=ttheme_minimal(base_size=8)
)

dev.off()


# Create power-plot with all scenarios
res$Scenario <- with(
  res,
  paste0(
    "HR_M=", HR_M,
    ", HR_DRM=", HR_DRM,
    ", HR_INT=", HR_INT
  )
)

res$Dropout <- paste0(
  100 * res$dropout_p,
  "%"
)

ggplot(res, aes(N, Power, color = Scenario)) +
  geom_line() +
  geom_point() +
  geom_hline(
    yintercept = target_power,
    linetype = 2
  ) +
  facet_wrap(~ Dropout) +
  labs(
    x = "Patients recruited",
    y = "Power",
    color = "Scenario"
  ) +
  theme_minimal()

ggsave(
  "Power_Cox.pdf",
  width=10,
  height=7
)

# Save Workspace
save.image("simulation_results.RData")
