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

data {
  int<lower=1> T;
  int<lower=1> R;
  array[T] vector[R] p_obs;
  array[T] vector[R] N;
  real N0;
  real p0;
  array[T] real ts;
  real psiT;
  real psiR;
}

parameters {
 // real<lower =0, upper = 1> rho;                       // Parameter rho
  real<upper = 0> log10_gamma;  // Log-transformed gamma
  
}

transformed parameters {
  real gamma =pow(10,-7) ;//pow(10, log10_gamma); // Back-transform log10_gamma to gamma
  real rho = 0.1;
  vector[T] y_sim =[0.1,.1,.1]';                 // Predicted values for the equation of interest
  matrix[R, T] sigma;               // Standard deviations for each replicate and time
  vector[2] y0;
  y0 = [.9,437037]';
  
  array[T]vector[2] ode_solution = ode_rk45(ode_system, y0, -10^-7, ts,rho, gamma, psiT, psiR);
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
// 
// parameters {
//   vector[2] y0;
//   real<lower=0> sigma;
//   real<lower=0> tau;
//   real rho;
//   real gamma;
// }
// 
// model {
//   sigma ~ normal(0, 0.5);
//   tau ~ normal(0,1000);
//   rho ~ uniform(0, 1);
//   gamma ~ normal(-15, 2.5);
//   y0[1] ~ uniform(0.95*p0,1.05*p0);
//   y0[2] ~ normal(N0, 0.05 * N0);
//   
//   array[T] vector[2] mu = ode_rk45(sho, y0, -10^-7, ts, rho, exp(gamma), psiT, psiR);
//   for (t in 1:T) {
//     
//     for(r in 1:R) {
//       y[t,r] ~ normal(mu[t,1], sigma);
//       N[t,r] ~ normal(mu[t,2], tau);
//     }
//     
//   }
//   
// }
// 
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



