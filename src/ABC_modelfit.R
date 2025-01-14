# Load necessary packages
library(deSolve)  # For solving ODEs
library(EasyABC)  # For Approximate Bayesian Computation
library(ggplot2)
library(tidyverse)
library("GGally")

#run it for each of the experiments ####
data.conj <- readxl::read_xlsx("./data/Stability_Egil_data.xlsx", sheet = "For_R")
data.growth.rates <- readxl::read_xlsx("./data/Stability_Egil_data.xlsx", sheet = "Growth_rates_for_R")
data.conj[data.conj$Time == 0,] <- data.conj%>%filter(Time == 0)%>%mutate(Nax = Nax / 100,
                                                                          CFX = CFX / 100)#at time 0 a different solution was used

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

# Define the coupled differential equations
dc_dr_dt <- function(time, state, parameters) {
  x <- state[1] # Current value of c
  r <- state[2] # Current value of r
  with(as.list(parameters), {
    # ODEs
    dx <- psiT * x - rho * psiT * x + gamma * x * r 
    dr <- psiR * r + rho * psiT * x - gamma * x * r
    list(c(dx, dr))
  })
}


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


# Simulate the model (a wrapper)
simulate_model_size <- function(params, time_points = c(0,1,2)) {
  # Parameters to simulate
  # Parameters
  parameters <- list(
    psiT = 2.7,#unlist(psiT),  # Example psiT
    psiR = 2.7,#unlist(psiR),  # Example psiR
    rho = params[1],   # Example rho
    gamma = 10^params[2] # Example gamma
  )
  
  # Initial conditions
  initial_state <- c(
    x = mean(cfx_obs[,1]),
    r = N0-mean(cfx_obs[,1])  
  )
  
  
  # Solve ODE
  out <- ode(y = initial_state, times = time_points, func = dc_dr_dt, parms = parameters)
  return(as.data.frame(out))
}


simulate_model_prop <- function(params, time_points = c(0,1,2)) {
  # Parameters to simulate
  # Parameters
  parameters <- list(
    psiT = 2.7,#unlist(psiT),  # Example psiT
    psiR = 2.7,#unlist(psiR),  # Example psiR
    rho = params[1],   # Example rho
    gamma = 10^params[2] # Example gamma
  )
  
  # Initial conditions
  initial_state <- c(
    p = mean(cfx_obs[,1])/N0,
    r = N0  
  )
  
  
  # Solve ODE
  out <- ode(y = initial_state, times = time_points, func = dp_dn_dt, parms = parameters)
  return(as.data.frame(out))
}


sim.size <- simulate_model_size(c(0.1, -7))
sim.prop <- simulate_model_prop(c(0.1, -7))


# Observed data
observed_data_size <- data.frame(
  time = rep(c(0, 1, 2), each =3),
  x_mean = c(cfx_obs),
  r_mean = c(N - cfx_obs )
)

observed_data_prop <- data.frame(
  time = rep(c(0, 1, 2), each =3),
  p_mean = c(p_obs),
  n_mean = c(N)
)

# Distance function to compare simulated and observed data
distance_function_size <- function(params) {
  sim_data_size <- simulate_model(params)
  # Compute the sum of squared errors (SSE)
  sse <- sum((sim_data$x - observed_data$x_mean)^2 + (sim_data$r - observed_data$r_mean)^2)
  return(sse)
}

distance_function_prop <- function(params) {
  sim_data_prop <- simulate_model_prop(params)
  # Compute the sum of squared errors (SSE)
  sse <- sum(sapply(c(1:3), function(i){sum((sim_data_prop$p[i] - observed_data_prop$p_mean[c(1:3)+(i-1)*3])^2)})) 
  return(sse)
}

steps <- 100
range_rho  <- c(0.00, 0.2)
range_loggamma <- c(-10,-6)
search_area <- data.frame(rho = rep(seq(range_rho[1],range_rho[2],diff(range_rho)/steps), each =steps+1),
                          loggamma = rep(seq(range_loggamma[1],range_loggamma[2],diff(range_loggamma)/steps), steps+1))

area.high.fit<- data.frame(t(mapply(FUN = function(rho,gamma){c(rho, gamma,distance_function_prop(c(rho,gamma)))},
                       search_area$rho,
                       search_area$loggamma)))


ggplot(area.high.fit, aes(x = X1, y = X2, z = X3)) +
  geom_contour_filled(binwidth = 0.025)+
  theme_minimal() +                     # Clean theme
  labs(title = "2D Density Plot ~Sum of squares",
       x = "rho",
       y = "log10(gamma)")



# 
# ABC parameter inference
abc_results <- ABC_sequential(
  method = "Lenormand",
  model = distance_function_prop,
  prior = list(c("unif", 0.0, 1.0), c("unif", -15, -5)),  # Prior distributions for rho and log10(gamma)
  nb_simul = 500,  # Number of simulations
  summary_stat_target = 0.0,  # Placeholder for custom summary statistics
  p_acc_min = 0.1               # Minimum acceptance rate

)


# View results
summary(abc_results$param)

ggplot(data.frame(abc_results$param), aes(x = X1, y = X2)) +
  geom_point(alpha = 0.5, size = 1) +    # Add scatter plot of accepted parameters
  geom_density_2d(color = "blue") +     # Add 2D density contours
  stat_density_2d(aes(fill = ..level../500), geom = "polygon", alpha = 0.4) +  # Density heatmap
  scale_fill_viridis_c() +              # Color gradient for density
  theme_minimal() +                     # Clean theme
  labs(title = "2D Density Plot of Accepted Parameters",
       x = "rho",
       y = "log10(gamma)")



med_sim <- simulate_model_prop(apply(abc_results$param,2,median))

ggplot()+geom_path(aes(time, p),data = med_sim)+
  geom_point(aes(time, p), data = data.frame(time = rep(c(0,1,2),each = 3),p = c(p_obs)))



abc_MCMC_results <- ABC_mcmc(
  method = "Marjoram",
  model = distance_function_prop,
  prior = list(c("unif", 0.0, 1.0), c("unif", -15, -5)),  # Prior distributions for rho and log10(gamma)
  n_between_sampling= 50,  # Number of simulations between samples
  n_rec = 10000,
  summary_stat_target = 0.0  # Placeholder for custom summary statistics
  #p_acc_min = 0.1               # Minimum acceptance rate
  
)

mcmc.res <- data.frame(cbind(abc_MCMC_results$param, cat = abc_MCMC_results$dist))

ggplot(mcmc.res, aes(x = V1, y = V2, colour = cat, alpha = cat))+geom_point()

#View results
summary(abc_MCMC_results$param)

ggplot(data.frame(cbind(abc_MCMC_results$param,abc_MCMC_results$dist*abc_MCMC_results$stats_normalization)), aes(x = X1, y = X2, colour = X3)) +
  geom_point(alpha = 0.5, size = 1 ) +    # Add scatter plot of accepted parameters
 # geom_density_2d(color = "blue") +     # Add 2D density contours
#  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.4) +  # Density heatmap
 # scale_fill_viridis_c() +              # Color gradient for density
  theme_minimal() +                     # Clean theme
  labs(title = "2D Density Plot of Accepted Parameters",
       x = "rho",
       y = "log10(gamma)")



med_sim <- simulate_model_prop(apply(abc_MCMC_results$param,2,median))

ggplot()+geom_path(aes(time, p),data = med_sim)+
  geom_point(aes(time, p), data = data.frame(time = rep(c(0,1,2),each = 3),p = c(p_obs)))
