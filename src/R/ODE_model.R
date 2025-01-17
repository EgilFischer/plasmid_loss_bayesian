# Define the coupled differential equations
dp_dn_dt <- function(time, state, parameters) {
  p <- state[1] # Current value of p
  N <- state[2] # Current value of N
  with(as.list(parameters), {
    # ODEs
    dp <- (psiT - psiR) * p * (1 - p) - rho * psiT * p + gamma * N * p * (1 - p)
    dN <- psiR * p * N + psiT * (1 - p) * N
    list(c(dp, dN))
  })
}