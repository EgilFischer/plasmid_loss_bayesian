//differential equations defining proportion with plasmid and total population size
functions {
  vector ode_system(real t,
             vector y,
             real rho, 
             real gamma,
             real psiT, 
             real psiR) {
    vector[2] dydt;
    dydt[1] = (psiT - psiR) * y[1] * (1 - y[1]) - psiT * rho * y[1] + gamma * y[2] * y[1] * (1 - y[1]);
    dydt[2] = (psiT* y[1] + psiR *(1-y[1])) * y[2];
    return dydt;
  }
  
}

//input data
data {
  int<lower=1> T; //number of time points
  int<lower=1> R; //number of replicates
  array[T] vector[R] p_obs;//proportion observed with plasmid
  array[T] vector[R] N; //total population size
  real N0; //initial population size
  real p0; //initial proportion with plasmid
  array[T] real ts; //time points
  real psiT; // growth rate with plasmid
  real psiR; // growth rate without plasmid
}

// parameters to be estimated
parameters {
 // real<lower =0, upper = 1> rho;                       // Parameter rho
  real<upper = 0> log10_gamma;  // Log-transformed gamma
  
}

// transformed parameters
transformed parameters {
  real gamma =pow(10,-7) ;//pow(10, log10_gamma); // Back-transform log10_gamma to gamma
  real rho = 0.1;
  array[T] real y_sim;// =[0.1,.1,.1]';                 // Predicted values for the equation of interest
  matrix[R, T] sigma;               // Standard deviations for each replicate and time
  vector[2] y0;
  y0 = [.9,437037]';
  
 array[T] vector[2] ode_solution = ode_rk45(ode_system, y0, -10^-7, ts,rho, gamma, psiT, psiR);
  y_sim = ode_solution[,1];
  // array[T]vector[2] ode_solution = ode_rk45(ode_system, y0, -10^7, ts,  rho, gamma, psiT, psiR);
  // // Solve the ODE system
  //     for (t in 1:T) {
  //      y_sim[t] = ode_solution[t][1]; // Extract the specific ODE equation of interest (e.g., state[1])
  //    }
   
  // Compute standard deviations for each replicate and time point
     for (r in 1:R) {
       for (t in 1:T) {
         sigma[r, t] = 0.25 /3;//sqrt(y_sim[t] * (1 - y_sim[t]) / R);
       }
     }

}

//model to be fitted
model {
  // Priors
  // rho ~ normal(0, 1);                      // Example prior for rho
   log10_gamma ~ uniform(-15, -5); // Prior for log10_gamma
  
  // 
  // Likelihood
  for (r in 1:R) {
     for (t in 1:T) {
      p_obs[r, t] ~ normal(y_sim[t], sigma[r, t]); // Likelihood for replicate r at time t
     }
   }
}

// generate quantities to compare observed to simulation
generated quantities {
  vector[T] y_pred;
  vector[T] N_pred;
  array[T]vector[2] mu = ode_rk45(ode_system, y0, -10^-7, ts,rho, gamma, psiT, psiR);
  {
    for (t in 1:T) {
      y_pred[t] = mu[t,1];
      N_pred[t] = mu[t,2];
    }

  }

}



