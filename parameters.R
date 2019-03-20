# Declare parameters ----------------------------------------------------------
params_base <- list(
  threshold = 20000,
  dr = 0.035,
  u0 = 0.95,
  u1 = -0.002,
  u2 = -0.0005,
  v_sd = 0.9,
  v_pd = 0.6,
  c_treatment = 480,
  c_control = 200,
  c_progression = 3000,
  c_death = 5000,
  c_pd = 1000,
  lambda_control = 1.5,
  hr_treatment = 0.56,
  k1 = 2,
  a2 = 0.4,
  b2 = 0.1,
  mu3 = 0,
  sigma3 = 1
)
params_base$lambda_treatment <- with(
  params_base,
  lambda_control * hr_treatment ^ (-1 / k1)
)
params_base$r <- log(1 + params_base$dr)