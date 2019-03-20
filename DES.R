library(flexsurv) # for Gompertz distribution
library(dplyr)
library(tidyr)

# Source the parameters -------------------------------------------------------
# NOTE - This requires you to have set the working directory to the directory
# containing this file
source(file = "parameters.R")

# Define function to perform discrete event simulation ------------------------
do_DES <- function(params, n_DES = 1000, .reshape = FALSE, .summarise = FALSE) {
  with(params, {
    qaly <- function(a, b) {
      (exp(-r * a) * ((a ^ 2 * u2 + a * u1 + u0) * r ^ 2 + (2 * a * u2 + u1) *
                        r + 2 * u2) - exp(-r * b) * ((b ^ 2 * u2 + b * u1 + u0) * r ^ 2 + (2 * b *
                                                                                             u2 + u1) * r + 2 * u2)) / r ^ 3
    }
    
    DES <- data.frame(
      iter = seq(1, n_DES),
      X1.treatment = rweibull(n = n_DES,
                              scale = lambda_treatment,
                              shape = k1),
      X1.control = rweibull(n = n_DES,
                            scale = lambda_control,
                            shape = k1),
      X2 = rgompertz(n = n_DES,
                     rate = b2,
                     shape = a2),
      X3 = rlnorm(
        n = n_DES,
        meanlog = mu3,
        sdlog = sigma3
      )
    ) %>% mutate(
      progressed.treatment = (X1.treatment < X2),
      progressed.control = (X1.control < X2),
      LY_stable.treatment = pmin(X1.treatment, X2),
      LY_stable.control = pmin(X1.control, X2),
      LY.treatment = if_else(progressed.treatment, X1.treatment + X3, X2),
      LY.control = if_else(progressed.control, X1.control + X3, X2),
      QALY_stable.treatment = v_sd * qaly(0, LY_stable.treatment),
      QALY_stable.control = v_sd * qaly(0, LY_stable.control),
      QALY_progressive.treatment = if_else(
        progressed.treatment,
        v_pd * qaly(X1.treatment, X1.treatment + X3),
        0
      ),
      QALY_progressive.control = if_else(progressed.control,
                                         v_pd * qaly(X1.control, X1.control + X3),
                                         0),
      QALY.treatment = QALY_stable.treatment + QALY_progressive.treatment,
      QALY.control = QALY_stable.control + QALY_progressive.control,
      cost_stable.treatment = c_treatment / r * (1 - exp(-r * LY_stable.treatment)),
      cost_stable.control = c_control / r * (1 - exp(-r * LY_stable.control)),
      cost_progression.treatment = if_else(
        progressed.treatment,
        c_progression * exp(-r * X1.treatment),
        0
      ),
      cost_progression.control = if_else(progressed.control, c_progression *
                                           exp(-r * X1.control), 0),
      cost_death.treatment = c_death * exp(-r * LY.treatment),
      cost_death.control = c_death * exp(-r * LY.control),
      cost_progressive.treatment = if_else(progressed.treatment, c_pd / r * (exp(-r * X1.treatment) - exp(
        -r * (X1.treatment + X3)
      )), 0),
      cost_progressive.control = if_else(progressed.control, c_pd / r * (exp(-r *
                                                                               X1.control) - exp(
                                                                                 -r * (X1.control + X3)
                                                                               )), 0),
      cost.treatment = cost_stable.treatment + cost_progression.treatment + cost_death.treatment + cost_progressive.treatment,
      cost.control = cost_stable.control + cost_progression.control + cost_death.control + cost_progressive.control
    )
    if (.reshape) {
      DES <- DES %>% gather("key", "value",-iter,-X2,-X3) %>%
        separate(key, c("variable", "arm"), sep = "\\.") %>%
        mutate_at(vars(arm), ~ factor(
          .,
          levels = c("treatment", "control", "comparison"),
          labels = c("Treatment", "Control", "Comparison")
        )) %>%
        spread(variable, value) %>%
        select(arm, iter, cost, QALY) %>%
        group_by(arm) %>%
        mutate(NMB = threshold * QALY - cost)
    }
    
    if (.summarise) {
      return(DES %>% summarise_all(list(
        mean = ~mean, se = ~ sd(.) / n()
      )))
    } else {
      return(DES)
    }
  })
}


# Apply the function to the parameters ----------------------------------------
results <- do_DES(params_base, .reshape = TRUE, .summarise = TRUE)
