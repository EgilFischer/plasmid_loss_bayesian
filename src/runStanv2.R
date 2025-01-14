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
# run it for each of the experiments ####
data_conjugation <- readxl::read_xlsx("./data/Raw/Stability_Egil_data.xlsx",
  sheet = "For_R"
)
data_growth_rates <- readxl::read_xlsx("./data/Raw/Stability_Egil_data.xlsx",
  sheet = "Growth_rates_for_R"
)

data_conjugation[data_conjugation$Time == 0, ] <- data_conjugation %>%
  filter(Time == 0) %>%
  mutate(
    Nax = Nax / 100,
    CFX = CFX / 100
  ) # at time 0 a different solution was used



# function to create a data input file for stan
create_stan_data <- function(strain, donor, data_growth_rates, data_conjugation) {
  #get the growth rates for this strain and donor
  #without plasmid
  psiR <- c(data_growth_rates %>%
    select(all_of(donor), all_of(paste0("Growth rates ", donor))) %>%
    filter(.data[[donor]] == strain) %>%
    select(paste0("Growth rates ", donor)))[1]
  #with plasmid
  psiT <- c(data_growth_rates %>%
    select(all_of(donor), all_of(paste0("Growth rates ", donor))) %>%
    filter(.data[[donor]] == paste0(strain, "-T")) %>%
    select(paste0("Growth rates ", donor)))[1]

  # Define data from conjugation experiment wiht donor and strain
  use_data <- data_conjugation %>%
    filter(Donor == donor & Strain_ID == strain);
  R <- length(unique(use_data$Replicate)) # Number of replicates
  time <- sort(unique(use_data$Time)) # Time points
  time_points <- length(time) # Number of time points
  
    #proportion of 
  p_obs <- matrix(
    use_data$CFX / use_data$Nax,
    nrow = R, byrow = FALSE
  )

  # population size
  population_size <- matrix(use_data$Nax,
    nrow = R, byrow = FALSE
  ) 

  # Package data for Stan
  stan_data <- list(
    T = time_points,
    R = R,
    ts = time,
    p_obs = p_obs,
    N = population_size,
    p0 = mean(p_obs[, 1]),
    N0 = mean(population_size[, 1]),
    psiT = unlist(psiT),
    psiR = unlist(psiR)
  )
}

#check stan_data
check_stan_data <- function(stan_data){
  #set an error message
  error_message = NULL;
  #check length of ts = T
  if(stan_data$T != length(stan_data$ts)){
    error_message = "\n ts does not have length T; \n ";
  }
  #check dimensions p_obs
  if(min(dim(stan_data$p_obs)==c(stan_data$R,stan_data$T))==0)
  {
    error_message = paste(error_message, "dim(p_obs) is not", stan_data$T, "x", stan_data$R ,"; \n")
  }
  #check dimensions N
  if(min(dim(stan_data$N)==c(stan_data$R,stan_data$T))==0)
  {
    error_message = paste(error_message, "dim(N) is not", stan_data$T, "x", stan_data$R ,";")
  }
  if(!is.null(error_message)) stop(error_message)
}

#create object for Stan 
stan_data <- create_stan_data(strain = "MH1",donor = "1D", 
                             data_growth_rates = data_growth_rates, 
                             data_conjugation = data_conjugation)

#to be changed but psi does not fit wiht N so change
stan_data$psiT<-2.7*stan_data$psiT/stan_data$psiR; 
stan_data$psiT<-2.7;

#check data
check_stan_data(stan_data)

# Stan model code (saved as a separate file, e.g., "model.stan")
stan_file <- "./src/ODE_plasmid_v2.stan"

# Compile the Stan model
stan_model <- stan_model(file = stan_file)

# Fit the model - currently it just simulates data as there is some bug in the estimation procedure.
fit <- sampling(
  stan_model,
  data = stan_data,
  iter = 5, # Number of iterations
  chains = 2, # Number of chains
  warmup = 1, # Number of warmup iterations
  thin = 1, # Thinning interval
  init = function() {
    list(
      rho = runif(1, 0, .1),
      log10_gamma = runif(1, -15, -5)
    )
  }, # force initial values of parameters
 # algorithm = "Fixed_param",
  seed = 123 # Random seed for reproducibility
)

# function to obtain the simulated data
print_sim_data <- function(sim, var_names) {
  for (name in var_names) {
    print(paste(name, ":", sim@sim$samples[[1]][name]))
  }
}


# print output of simulated data for y = proportion and N is total population size
print_sim_data(fit, c("y_pred[1]", "y_pred[2]", "y_pred[3]"))
print_sim_data(fit, c("N_pred[1]", "N_pred[2]", "N_pred[3]"))
# print all output
print_sim_data(fit, fit@sim$fnames_oi)




stan_diag(fit)

stan_trace(fit)

stan_hist(fit)

stan_ess(fit)

# Plot posterior predictive distribution
y_pred <- rstan::extract(fit, "y_pred")$y_pred
N_pred <- rstan::extract(fit, "N_pred")$N_pred

# Summarize predicted values
predicted_mean <- apply(y_pred, 2, mean)
predicted_sd <- apply(y_pred, 2, sd)

# Plot observed vs predicted
df <- with(stan_data, {data.frame(
  Time = rep(ts, each =R),
  Observed = as.vector(p_obs),
  Predicted = rep(predicted_mean, each =R)
)})

df.N<- with(stan_data, { data.frame(
  Time = rep(ts, each = R),
  Observed = as.vector(N),
  Predicted = rep(apply(N_pred, 2, mean), each = R)
)
})

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
    dp <- (psiT - psiR) * p * (1 - p) - rho * psiT * p + gamma * N * p * (1 - p)
    dN <- psiR * p * N + psiT * (1 - p) * N
    list(c(dp, dN))
  })
}

# Time points for the simulation
time_points <- seq(0, 2, length.out = 100) # Example time points

# Parameters
parameters <- list(
  psiT = unlist(psiT), # Example psiT
  psiR = unlist(psiR), # Example psiR
  rho = summary(fit)$summary["rho", 1], # Example rho
  gamma = exp(summary(fit)$summary["gamma", 1]) # Example gamma
)

# Initial conditions
initial_state <- c(
  p = mean(p_obs[, 1]), # summary(fit)$summary["y0[1]",1], # Initial value of p
  N = N0 # Initial value of N
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
