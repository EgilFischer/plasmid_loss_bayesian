# Estimate using Stan ####
library(rstan)
library(rstanarm)
library(ggplot2)
library(deSolve)
library(tidyverse)


# Set options for rstan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#user defined functions #### 
# function to create a data input file for stan
create_stan_data <- function(strain, donor, data_growth_rates, data_conjugation) {
  #' @title create data for the plasmid loss estimation procedure in Stan
  #' @description
    #' Extracts data for the specified strain and donor and creates a list as input for Stan model.
  #' @param strain name of the strain to extract data for (string)
  #' @param donor name of plasmid donor to extract data for (string)
  #' @param data_growth_rate data,frame containing bacterial growth rates
  #' @param data_conjugation data.frame containing conjugation data CFX and Nax columns
  #' @return stan_data list with inputs for Stan model:
    #'  
    #' T: number of time points
    #' 
    #' R: number of replicates
    #' 
    #' ts: time points
    #' 
    #' p_obs: observed proportion with plasmid
    #' 
    #' N: observed population size
    #' 
    #' p0: observed mean initial proportion with plasmid
    #' 
    #' N0:  observed mean initial population size
    #' 
    #' psiT: input growth rate with plasmid
    #' 
    #' psiR: input growth rate without plasmid  
    #' 
  
  #get the growth rates for the specified strain and donor
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
  
  #deal with unknown growth rates
  if(length(unlist(psiR))==0) {
    warning(paste(strain, "and" , donor, "combination does not exist in the growth rate data"));
  psiR <- 1; psiT<- 1}
  
  # Define data from conjugation experiment with specified donor and strain
  use_data <- data_conjugation %>%
    filter(Donor == donor & Strain_ID == strain);
  R <- length(unique(use_data$Replicate)) # Number of replicates
  time <- sort(unique(use_data$Time)) # Time points
  time_points <- length(time) # Number of time points
  
  #proportion plasmid carrying 
  p_obs <- matrix(
    use_data$CFX / use_data$Nax,
    nrow = R, byrow = FALSE
  )
  
  #total population size
  population_size <- matrix(use_data$Nax,
                            nrow = R, byrow = FALSE
  ) 
  
  # Package data for Stan
  stan_data <- list(
    T = time_points, #number of time points
    R = R, #number of replicates
    ts = time, #time points
    p_obs = p_obs, #observed proportion with plasmid
    N = population_size, #observed population size
    p0 = mean(p_obs[, 1]), #observed mean initial proportion with plasmid
    N0 = mean(population_size[, 1]), #observed mean initial population size
    psiT = unlist(psiT), #input growth rate with plasmid
    psiR = unlist(psiR) #input growth rate without plasmid
  )
  
  return(stan_data)
}

#check if the data in stan_data has the right dimensions.
check_stan_data <- function(stan_data){
  #' @title Data dimension check stan data
  #' @description
    #' Checks the dimension of observed data. It willt throw an error if the dimensions are incorrect 
  #' @param stan_data is a list produced by create_stan_data function
  #' @seealso create_stan_data
    
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
  
  #the error message is build up during each step. If no errors are found it will still be NULL otherwise it throws an error.
  if(!is.null(error_message)) stop(error_message)
}

# function to obtain the simulated data
print_sim_data <- function(sim, var_names) {
  #' @title Print data simulated by Stan model
  #' @param sim output of Stan fit function
  #' @param var_names names to get output from  
  for (name in var_names) {
    print(paste(name, ":", sim@sim$samples[[1]][name]))
  }
}

# load data ####
data_conjugation <- readxl::read_xlsx("./data/Stability_Egil_data.xlsx",
  sheet = "For_R"
)
data_growth_rates <- readxl::read_xlsx("./data/Stability_Egil_data.xlsx",
  sheet = "Growth_rates_for_R"
)

# at time 0 a different solution was used therefore we need to divide by 100
data_conjugation[data_conjugation$Time == 0, ] <- data_conjugation %>%
  filter(Time == 0) %>%
  mutate(
    Nax = Nax / 100,
    CFX = CFX / 100
  ) 




#create object for Stan #### 
stan_data <- create_stan_data(strain = "MH1",donor = "2D", 
                             data_growth_rates = data_growth_rates, 
                             data_conjugation = data_conjugation)

#TO DO: to be changed after discussion with empiricist but psi does not fit with N so change
stan_data$psiR<-2.7*stan_data$psiR/stan_data$psiT; 
stan_data$psiT<-2.7;

#check data
check_stan_data(stan_data)

# Stan model code (saved as a separate file, e.g., "model.stan")
stan_file <- "./src/Stan/ODE_plasmid_v2.stan"
# 
# # Compile the Stan model ~do not care about message : "hash mismatch so recompiling; make sure Stan code ends with a blank line"
# stan_model <- stan_model(file = stan_file)
# 
# # Fit the model - currently it just simulates data as there is some bug in the estimation procedure.
# rm(fit);fit <- sampling(
#   stan_model,
#   data = stan_data,
#   iter = 20000, # Number of iterations
#   chains = 1, # Number of chains
#   warmup = 7500, # Number of warmup iterations
#   thin = 1, # Thinning interval,
#   algorithm = "NUTS",
#   control = list(max_treedepth = 12),
#   init = function() {
#     list(
#       rho = runif(1, 0, .1),
#       log10_gamma = runif(1, -15, -3)
#     )
#   },
#   
#   # force initial values of parameters
#  # algorithm = "Fixed_param",
#   seed = 123 # Random seed for reproducibility
# )


# print output of simulated data for y = proportion and N is total population size
#print_sim_data(fit, c("y_pred[1]", "y_pred[2]", "y_pred[3]"))
#print_sim_data(fit, c("N_pred[1]", "N_pred[2]", "N_pred[3]"))
# print all output
#print_sim_data(fit, fit@sim$fnames_oi)



# model diagnostics

stan_trace(fit, pars = list("rho","log10_gamma"))

stan_hist(fit, pars = list("rho","log10_gamma"))

stan_ess(fit, pars = list("rho","log10_gamma"))

pairs(fit,pars = c("rho" , 
                   "log10_gamma", 
                   "gamma"))

# extract posterior predictive distribution
y_pred <- rstan::extract(fit, "y_pred")$y_pred
N_pred <- rstan::extract(fit, "N_pred")$N_pred

# Summarize predicted values
predicted_mean <- apply(y_pred, 2, mean)

predicted_sd <- apply(y_pred, 2, sd)

# Plottable data for observed vs predicted
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


#plots ####
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


#run Stan over multipe Strains and donor combinations ####

#function iterating over strain_donor combinations
stan_fit_strain_donor <- function(strain_donor_combinations, #data frame containing combination of strain and donor name
                                  data_conjugation, #data with conjugation experiment data
                                  data_growth_rates, #data with growth rate experiment data
                                  iter = 30000,
                                  warmup = 15000
){
  #' @title Fit the ode-system for strain-donor combinations
  #'
  #' This function processes strain-donor combinations and related experimental data
  #' (e.g., conjugation and growth rates) to fit a model or perform specific analyses.
  #'
  #' @param strain_donor_combinations Data frame. A table containing combinations of strain and donor names.
  #' @param data_conjugation Data frame. Experimental data from conjugation experiments.
  #' @param data_growth_rates Data frame. Experimental data from growth rate experiments.
  #' @param iter Number of iterations Default 30 000
  #' @param warmup Number of warmup Default 15 000 
  #' @return Model or processed data. The output depends on the implementation of the function.
  #' @examples
  #' strain_donor_combinations <- data.frame(strain = c("strain1", "strain2"),
  #'                                          donor = c("donor1", "donor2"))
  #' data_conjugation <- data.frame(time = c(1, 2, 3), conjugation_rate = c(0.1, 0.2, 0.3))
  #' data_growth_rates <- data.frame(strain = c("strain1", "strain2"), growth_rate = c(1.2, 1.5))
  #' stan_fit_strain_donor(strain_donor_combinations, data_conjugation, data_growth_rates)
  
  #iterate over the combinations of strains and donors
  output <- NULL;
  for(i in c(1:dim(strain_donor_combinations)[1])){
  strain_donor <- c(strain_donor_combinations[i,]);
  print(strain_donor);
  #create the stan_data 
  stan_data <- create_stan_data(strain = unlist(strain_donor$Strain_ID),
                                donor = unlist(strain_donor$Donor),
                                data_growth_rates = data_growth_rates,
                                data_conjugation = data_conjugation);
  check_stan_data(stan_data)
  #determine the growth rate in this data set
  overall_psi <- coefficients(lm(logN~t, data = data.frame(logN = c(log(stan_data$N)), 
                                          t = rep(stan_data$ts ,each = stan_data$R))))["t"]
  #adjust to ratio psiT to psiR 
  alpha <- with(stan_data, psiT/psiR)
  stan_data$psiR <-with(stan_data,{ overall_psi/((1-p0) *alpha+p0) })
  stan_data$psiT <- alpha*stan_data$psiR
  
  #run the fitting procedure
  {rm(fit);fit <- sampling(
    stan_model,
    data = stan_data,
    iter = iter, # Number of iterations
    chains = 1, # Number of chains
    warmup = warmup, # Number of warmup iterations
    thin = 1, # Thinning interval,
    algorithm = "NUTS",
    control = list(max_treedepth = 12),
    init = function() {
      list(
        rho = runif(1, 0, .1),
        log10_gamma = runif(1, -15, -3)
      )
    },
    
    # force initial values of parameters
    # algorithm = "Fixed_param",
    seed = 123 # Random seed for reproducibility
  )}
  
  #if fit contains values
  if(is.null(fit)) stop()
  #extract relevant data
  output <- rbind(output,data.frame(strain = strain_donor[1],
                       donor = strain_donor[2],
                       rho = rstan::extract(fit,"rho"),
                       log10_gamma = rstan::extract(fit,"log10_gamma"),
                       prop_predict_t0 = rstan::extract(fit, "y_pred")$y_pred[,1],
                       prop_predict_t1 = rstan::extract(fit, "y_pred")$y_pred[,2],
                       prop_predict_t2 = rstan::extract(fit, "y_pred")$y_pred[,3]))
  
  }
  return(output)
}


# All combinations of donor and strain ------------------------------------
#get the donor strain combinations
strain_donor_combinations <- unique(data_conjugation[,c("Strain_ID","Donor")])

#remove a few with unclear data 
strain_donor_combinations <- strain_donor_combinations%>%filter(!(Strain_ID %in% c("4H","5H","8H","4B3")))



stan_fit_strain_donor_output <- stan_fit_strain_donor(strain_donor_combinations, #data frame containing combination of strain and donor name
                                  data_conjugation, #data with conjugation experiment data
                                  data_growth_rates, #data with growth rate experiment data
                                  iter = 5000,
                                  warmup = 1000
)

saveRDS(stan_fit_strain_donor_output, file = paste0("./results/",gsub("-","",Sys.Date()),"_stan_fit_strain_donor_output.RDS"))

stan_fit_strain_donor_output <- readRDS(file = paste0("./results/20250116_stan_fit_strain_donor_output.RDS"))

prior_distribution <- data.frame(rho = runif(10^6,0,1 ),
                                 log10_gamma =  rnorm(10^6,-7,3 ))

ggplot(stan_fit_strain_donor_output)+
   geom_histogram(aes(x = rho, 
                      y = after_stat(density), 
                      fill = paste(Strain_ID,Donor)))+
  geom_density(data = prior_distribution,
                aes(x = rho,y = after_stat(density)),
                colour = "red", linewidth = 0.5)+
  facet_grid(Strain_ID~Donor)+
  ggtitle("rho")
                
                

ggplot(stan_fit_strain_donor_output)+
  geom_histogram(aes(x = log10_gamma, 
                     y = after_stat(density), 
                     fill = paste(Strain_ID,Donor)))+
  geom_density(data = prior_distribution,
                aes(x = log10_gamma,y = after_stat(density)),
                colour = "red", linewidth = .5)+
  facet_grid(Strain_ID~Donor)+
  ggtitle("log10_gamma")


ggplot(stan_fit_strain_donor_output)+
  geom_bin2d(aes(x = log10_gamma, 
                     y = rho), bins = 25)+
  scale_fill_continuous(type = "viridis")+
  facet_grid(Strain_ID~Donor)+
  ggtitle("log10_gamma vs rho")

#organize data for comparison
comp_data <- data_conjugation%>%mutate(p_obs = CFX/Nax)

#calculate mean square error ####
mse_fit <- function(comp_data, stan_fit_strain_donor_output){
  #number of entries
  n <- length(comp_data$Replicate);
  #initialize output 
  mse <- data.frame(mse = rep(Inf, n))
  #iterate over the comparison data
  for(i in c(1:n))
  {
      strain <- comp_data$Strain_ID[i];
      donor <- comp_data$Donor[i];
      time <- comp_data$Time[i];
      stan_fit<-stan_fit_strain_donor_output %>% filter(Strain_ID == strain & Donor == donor) 
      if(length(stan_fit)<1)
        mse$mse[i]<- Inf
      else
      {
         mse$mse[i] <- sum((stan_fit[,time + 5] -    comp_data$p_obs[i])^2 )
      }
      
  }
  return(cbind(comp_data, mse))
  
}

comp_data_mse = mse_fit(comp_data, stan_fit_strain_donor_output )
#combine for all time points
comp_data_mse_aggregate <- comp_data_mse%>%reframe(.by = c(Strain_ID,Replicate,Donor, Origin),
                                                   sum_mse = sum(mse))

ggplot(comp_data_mse_aggregate)+geom_point(aes(x = Strain_ID, y = sum_mse, colour = as.factor(Replicate)))+facet_grid(Donor~.)

ggplot(stan_fit_strain_donor_output)+
  geom_point(aes(x = 0, y = prop_predict_t0), alpha = 0.01)+
  geom_point(aes(x = 1, y = prop_predict_t1, alpha = 0.01))+
  geom_point(aes(x = 2, y = prop_predict_t2), alpha = 0.01)+
  geom_path(aes(x = Time, y = p_obs, group = Replicate), data =comp_data)+
  facet_grid(Strain_ID~Donor,scales = "free_y")
  


# 
# # In this section solve the ODE's without Stan ####
# source("./src/R/ODE_model.R")
# 
# 
# # Time points for the simulation
# time_points <- seq(0, 2, length.out = 100) # Example time points
# 
# # Parameters
# parameters <- list(
#   psiT = unlist(stan_data$psiT), # Example psiT
#   psiR = unlist(stan_data$psiR), # Example psiR
#   rho = summary(fit)$summary["rho", 1], # Example rho
#   gamma = summary(fit)$summary["gamma", 1] # Example gamma
# )
# 
# # Initial conditions
# initial_state <- c(
#   p = stan_data$p0,  # Initial value of p
#   N = stan_data$N0  # Initial value of N
# )
# 
# # Solve the coupled ODEs
# ode_solution <- ode(
#   y = initial_state,
#   times = time_points,
#   func = dp_dn_dt,
#   parms = parameters
# )
# 
# # Extract the results
# ode_results <- as.data.frame(ode_solution)
# 
# ggplot() +
#   geom_point(aes(x = df$Time, y = df$Observed), color = "blue", size = 2) +
#   geom_line(aes(x = ode_results$time, y = ode_results$p), color = "red", size = 1) +
#   theme_minimal() +
#   labs(
#     title = "Observed vs Predicted",
#     x = "Time",
#     y = "Proportion"
#   )
# ggplot() +
#   geom_line(aes(x = ode_results$time, y = ode_results$N), color = "red", size = 1) +
#   theme_minimal() +
#   labs(
#     title = "Observed vs Predicted",
#     x = "Time",
#     y = "Proportion"
#   )
