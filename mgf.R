library(flexsurv) # for Gompertz distribution
library(dplyr)    # for mutate and transmute
library(purrr)    # for map and partial

# Source the parameters -------------------------------------------------------
# NOTE - This requires you to have set the working directory to the directory
# containing this file
source(file = "parameters.R")

# Define function to perform moment-generating function method ----------------
do_MGF <- function(params) {
  with(params, {
    # Convenience functions
    f_X1_treatment <- function(x) dweibull(x, scale = lambda_treatment, shape = k1)
    F_X1_treatment <- function(x) pweibull(x, scale = lambda_treatment, shape = k1)
    f_X1_control   <- function(x) dweibull(x, scale = lambda_control, shape = k1)
    F_X1_control   <- function(x) pweibull(x, scale = lambda_control, shape = k1)
    f_X2           <- function(x) dgompertz(x, rate = b2, shape = a2)
    F_X2           <- function(x) pgompertz(x, rate = b2, shape = a2)
    
    semi_infinite_integral <- partial(integrate, lower = 0, upper = Inf, rel.tol = 1e-8)
    siiv <- function(...) semi_infinite_integral(...)$value
    
    MGF_qaly <- function(mgf_a_0, mgf_a_1, mgf_a_2,
                         mgf_b_0, mgf_b_1, mgf_b_2) {
      1 / r ^ 3 * ((
        (u0 * mgf_a_0 + u1 * mgf_a_1 + u2 * mgf_a_2) * r ^ 2 +
          (u1 * mgf_a_0 + 2 * u2 * mgf_a_1) * r +
          (2 * u2 * mgf_a_0)
      ) - (
        (u0 * mgf_b_0 + u1 * mgf_b_1 + u2 * mgf_b_2) * r ^ 2 +
          (u1 * mgf_b_0 + 2 * u2 * mgf_b_1) * r +
          (2 * u2 * mgf_b_0)
      ))
    }
    
    p_treatment <- siiv(function(x) f_X2(x) * F_X1_treatment(x))
    p_control   <- siiv(function(x) f_X2(x) * F_X1_control(x))
    
    GM_X1_treatment <- map_dbl(
      0:2,
      ~ siiv(
        function(x)
          x ^ (.) * exp(-r * x) * f_X1_treatment(x) * (1 - F_X2(x)) / p_treatment))
    GM_X1_control <- map_dbl(
      0:2,
      ~ siiv(
        function(x) x ^ (.) * exp(-r * x) * f_X1_control(x) * (1 - F_X2(x)) / p_control))

    GM_X2_treatment <- map_dbl(
      0:2,
      ~ siiv(
        function(x)
          x ^ (.) * exp(-r * x) * f_X2(x) * (1 - F_X1_treatment(x)) / (1 - p_treatment)
        ))
    GM_X2_control <- map_dbl(
      0:2,
      ~ siiv(
        function(x)
          x ^ (.) * exp(-r * x) * f_X2(x) * (1 - F_X1_control(x)) / (1 - p_control)
        ))
    
    GM_X3 <- map_dbl(
      0:2,
      ~ siiv(
        function(x)
          x ^ (.) * exp(-r * x) * dlnorm(x, meanlog = mu3, sdlog = sigma3)
        ))
    
    MGF <- data.frame(
      arm = factor(c("Treatment", "Control")),
      cost_stable = c(
        c_treatment / r * (
          1 - p_treatment * GM_X1_treatment[1] - (1 - p_treatment) * GM_X2_treatment[1]
        ),
        c_control / r * (
          1 - p_control * GM_X1_control[1] - (1 - p_control) * GM_X2_control[1]
        )
      ),
      cost_progression = c(p_treatment, p_control) * c_progression *
        c(GM_X1_treatment[1], GM_X1_control[1]),
      cost_death = c_death * c(
        p_treatment * GM_X1_treatment[1] * GM_X3[1] + (1 - p_treatment) * GM_X2_treatment[1],
        p_control * GM_X1_control[1] * GM_X3[1] + (1 - p_control) * GM_X2_control[1]
      ),
      cost_progressive = c(p_treatment, p_control) * c_pd / r *
        c(GM_X1_treatment[1], GM_X1_control[1]) * (1 - GM_X3[1]),
      QALY_stable = v_sd * c(
        p_treatment * MGF_qaly(
          1,
          0,
          0,
          GM_X1_treatment[1],
          GM_X1_treatment[2],
          GM_X1_treatment[3]
        ) +
          (1 - p_treatment) * MGF_qaly(
            1,
            0,
            0,
            GM_X2_treatment[1],
            GM_X2_treatment[2],
            GM_X2_treatment[3]
          ),
        p_control * MGF_qaly(1, 0, 0, GM_X1_control[1], GM_X1_control[2], GM_X1_control[3]) +
          (1 - p_control) * MGF_qaly(1, 0, 0, GM_X2_control[1], GM_X2_control[2], GM_X2_control[3])
      ),
      QALY_progressive = v_pd * c(p_treatment, p_control) * c(
        MGF_qaly(
          GM_X1_treatment[1],
          GM_X1_treatment[2],
          GM_X1_treatment[3],
          GM_X1_treatment[1] * GM_X3[1],
          GM_X1_treatment[2] * GM_X3[1] + GM_X1_treatment[1] * GM_X3[2],
          GM_X1_treatment[3] * GM_X3[1] + 2 * GM_X1_treatment[2] * GM_X3[2] + GM_X1_treatment[1] *
            GM_X3[3]
        ),
        MGF_qaly(
          GM_X1_control[1],
          GM_X1_control[2],
          GM_X1_control[3],
          GM_X1_control[1] * GM_X3[1],
          GM_X1_control[2] * GM_X3[1] + GM_X1_control[1] * GM_X3[2],
          GM_X1_control[3] * GM_X3[1] + 2 * GM_X1_control[2] * GM_X3[2] + GM_X1_control[1] *
            GM_X3[3]
        )
      )
    ) %>% transmute(
      arm = arm,
      cost = cost_stable + cost_progression + cost_death + cost_progressive,
      QALY = QALY_stable + QALY_progressive,
      NMB = QALY * threshold - cost
    )
    return(MGF)
  })
  
}

# Apply the function to the parameters ----------------------------------------
results <- do_MGF(params_base)
