# Estimate using Stan ####
library(rstan)
library(rstanarm)
library(ggplot2)
library(deSolve)
library(tidyverse)


# Set options for rstan 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Plasmid loss model ####
#run it for each of the experiments ####
data.conj <- readxl::read_xlsx("./data/Stability_Egil_data.xlsx", sheet = "For_R")
data.growth.rates <- readxl::read_xlsx("./data/Stability_Egil_data.xlsx", sheet = "Growth_rates_for_R")
data.conj[data.conj$Time == 0,] <- data.conj%>%filter(Time == 0)%>%mutate(Nax = Nax / 100,
                                       CFX = CFX / 100)#at time 0 a different solution was used
# Stan model code (saved as a separate file, e.g., "model.stan")
stan_file <- "src/ODE_plasmid_v2.stan"




#select growth rates
strain <- "MH1"
donor <- "1D"
psiR <- c(data.growth.rates%>%
  select(all_of(donor), all_of(paste0("Growth rates ", donor)))%>%
  filter(.data[[donor]] == strain)%>%select(paste0("Growth rates ", donor)) )[1]
psiT <- c(data.growth.rates%>%
            select(all_of(donor), all_of(paste0("Growth rates ", donor)))%>%
            filter(.data[[donor]] == paste0(strain,"-T"))%>%select(paste0("Growth rates ", donor)) )[1]

# Define data
T <- 3  # Number of time points
R <- 3  # Number of replicates
time <- c(0, 1, 2)  # Time points
p_obs <- matrix((data.conj%>%
                   filter(Donor == donor & Strain_ID==strain))$CFX/(data.conj%>%
                                                                      filter(Donor == donor & Strain_ID==strain))$Nax, 
                nrow = R, byrow = FALSE)
cfx_obs <-matrix((data.conj%>%
                    filter(Donor == donor & Strain_ID==strain))$CFX, 
                 nrow = R, byrow = FALSE)


N <- matrix((data.conj%>%filter(Donor == donor & Strain_ID==strain))$Nax, nrow = R, byrow = FALSE) # Input variable N at each time point
N0 <- mean(N[,1])

# Package data for Stan
stan_data <- list(
  T = T,
  R = R,
  ts = time,
  p_obs = p_obs,
  N = N,
  p0 = mean(p_obs[,1]),
  N0 = N0,
  psiT = unlist(psiT),
  psiR = unlist(psiR)
  
)


# Compile the Stan model
stan_model <- stan_model(file = stan_file)
# Fit the model
rm(fit)


fit <- sampling(
  stan_model,
  data = stan_data,
  iter = 1,     # Number of iterations
  chains = 1,      # Number of chains
  warmup = 0,    # Number of warmup iterations
  thin = 1,        # Thinning interval
  init = function() {
    list(rho = runif(1, 0, .1), 
         log10_gamma = runif(1, -15, -5))
  }, #force initial values of parameters
  algorithm = "Fixed_param",
  seed = 123       # Random seed for reproducibility
  
)

fit@sim$samples[[1]]$`y_pred[1]`
fit@sim$samples[[1]]$`y_pred[2]`
fit@sim$samples[[1]]$`y_pred[3]`


fit@sim$samples[[1]]$`N_pred[1]`
fit@sim$samples[[1]]$`N_pred[2]`
fit@sim$samples[[1]]$`N_pred[3]`


stan_diag(fit)

stan_trace(fit)

stan_hist(fit)

stan_ess(fit)

# Plot posterior predictive distribution
y_pred <- rstan::extract(fit,"y_pred")$y_pred
N_pred <- rstan::extract(fit,"N_pred")$N_pred

# Summarize predicted values
predicted_mean <- apply(y_pred, 2, mean)
predicted_sd <- apply(y_pred, 2, sd)

# Plot observed vs predicted
df <- data.frame(
  Time = rep(time, each = R),
  Observed = as.vector(p_obs),
  Predicted = rep(predicted_mean, each = R)
)

df.N <- data.frame(
  Time = rep(time, each = R),
  Observed = as.vector(N),
  Predicted = rep(apply(N_pred, 2, mean), each = R)
)


ggplot(df, aes(x = Time)) +
  geom_point(aes(y = Observed), color = "blue", size = 2) +
  geom_line(aes(y = Predicted), color = "red", size = 1) +
  theme_minimal() +
  labs(
    title = "Observed vs Predicted",
    x = "Time",
    y = "Proportion"
  )


ggplot(df.N, aes(x = Time)) +
  geom_point(aes(y = Observed), color = "blue", size = 2) +
  geom_line(aes(y = Predicted), color = "red", size = 1) +
  theme_minimal() +
  labs(
    title = "Observed vs Predicted",
    x = "Time",
    y = "N"
  )




# Define the coupled differential equations
dp_dn_dt <- function(time, state, parameters) {
  p <- state[1] # Current value of p
  N <- state[2] # Current value of N
  with(as.list(parameters), {
    # ODEs
    dp <- (psiT - psiR) * p * (1 - p) - rho * psiT* p + gamma * N * p * (1 - p)
    dN <- psiR * p * N + psiT * (1 - p) * N
    list(c(dp, dN))
  })
}

# Time points for the simulation
time_points <- seq(0, 2, length.out = 100) # Example time points

# Parameters
parameters <- list(
  psiT = unlist(psiT),  # Example psiT
  psiR = unlist(psiR),  # Example psiR
  rho = summary(fit)$summary["rho",1],   # Example rho
  gamma = exp(summary(fit)$summary["gamma",1]) # Example gamma
)

# Initial conditions
initial_state <- c(
  p = mean(p_obs[,1]),#summary(fit)$summary["y0[1]",1], # Initial value of p
  N = N0  # Initial value of N
)

# Solve the coupled ODEs
ode_solution <- ode(
  y = initial_state,
  times = time_points,
  func = dp_dn_dt,
  parms = parameters
)

# Extract the results
ode_results <- as.data.frame(ode_solution)

ggplot() +
  geom_point(aes(x = df$Time, y = df$Observed), color = "blue", size = 2) +
  geom_line(aes(x = ode_results$time, y = ode_results$p), color = "red", size = 1) +
  theme_minimal() +
  labs(
    title = "Observed vs Predicted",
    x = "Time",
    y = "Proportion"
  )
ggplot() +
  geom_line(aes(x = ode_results$time, y = ode_results$N), color = "red", size = 1) +
  theme_minimal() +
  labs(
    title = "Observed vs Predicted",
    x = "Time",
    y = "Proportion"
  )



