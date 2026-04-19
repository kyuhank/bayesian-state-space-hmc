// SBC / non-centered / differentiable soft constraints
// Constraint type:
//   1) Soft-min for effective catch on log scale via log_sum_exp
//   2) Softplus floor for next biomass median
// Notes:
//   - Non-centered process noise uses ut ~ N(0, 1) and Pt = PtMedian * exp(sigma * ut).
//   - smoothness controls how closely the approximation follows hard bounds.
//   - tau2Prior[3] toggles whether the index likelihood is included (1 = on, else off).
data {
  int ntimes;
  vector[ntimes] It;
  vector[ntimes] Ct;
  
  vector[2] logrPrior;
  vector[2] logKPrior;
  vector[2] sigma2Prior;
  vector[3] tau2Prior;
  
  vector[2] logqPrior;
}
parameters {
  real logr;
  real logK;
  real<lower=logqPrior[1], upper=logqPrior[2]> logq;
  real<lower=0> sigma2;
  real<lower=0> tau2;
  
  vector[ntimes + 1] ut;
}
transformed parameters {
  // Parameter transforms
  real sigma = sqrt(sigma2);
  real tau = sqrt(tau2);
  
  real r = exp(logr);
  real K = exp(logK);
  real q = exp(logq);
  
  real MSY = (r * K) / 4;
  real Hmsy = r / 2;
  real Bmsy = K / 2;
  
  vector[ntimes] HtoverHmsy;
  vector[ntimes + 1] BtoverBmsy;
  
  vector<lower=0>[ntimes + 1] Pt;
  vector<lower=0>[ntimes + 1] PtMedian;
  vector<lower=0>[ntimes + 1] Bt;
  vector[ntimes] Ht;
  vector[ntimes] PredCt;
  vector[ntimes] PredIt;
  
  // Higher values approach hard min/max behavior while staying differentiable.
  real smoothness = 50;
  
  real cdiff_scale = 0.01;
  
  // Initial condition on relative biomass median
  PtMedian[1] = 1;
  Pt[1] = PtMedian[1] * exp(sigma * ut[1]);
  Bt[1] = Pt[1] * K;
  
  for (t in 2 : (ntimes + 1)) {
    // Surplus production term before harvest removal
    real surplus_production = r * Pt[t - 1] * (1.0 - Pt[t - 1]);
    
    real log_expected_catch = log(Ct[t - 1]);
    real log_max_sustainable = log(K) + log(Pt[t - 1]) + log(0.99);
    
    // Differentiable approximation of min(Ct, 0.99 * K * Pt)
    real log_soft_min = -log_sum_exp(
                                     [smoothness * (-log_expected_catch),
                                      smoothness * (-log_max_sustainable)])
                        / smoothness;
    
    PredCt[t - 1] = exp(log_soft_min);
    
    real biomass_after_production = Pt[t - 1] + surplus_production;
    real harvest_impact = PredCt[t - 1] / K;
    real raw_next_biomass = biomass_after_production - harvest_impact;
    
    // Differentiable lower bound for biomass median
    PtMedian[t] = log1p_exp(smoothness * raw_next_biomass) / smoothness
                  + 1e-4;
    
    Pt[t] = PtMedian[t] * exp(sigma * ut[t]);
    
    Bt[t] = Pt[t] * K;
    Ht[t - 1] = PredCt[t - 1] / Bt[t - 1];
  }
  
  HtoverHmsy = Ht / Hmsy;
  BtoverBmsy = Bt / Bmsy;
  PredIt = q * Bt[1 : ntimes];
  
  // Penalty term that keeps effective catch close to input catch
  vector[ntimes] cdiff = (log(PredCt) - log(Ct)) / cdiff_scale;
}
model {
  // Priors
  logq ~ uniform(logqPrior[1], logqPrior[2]);
  logr ~ normal(logrPrior[1], logrPrior[2]);
  logK ~ normal(logKPrior[1], logKPrior[2]);
  
  sigma2 ~ inv_gamma(sigma2Prior[1], sigma2Prior[2]);
  tau2 ~ inv_gamma(tau2Prior[1], tau2Prior[2]);
  
  ut ~ std_normal();
  
  // Index likelihood can be disabled for prior-only experiments.
  
  if (tau2Prior[3] == 1) {
    It ~ lognormal(log(q * Bt[1 : ntimes]), tau);
  }
  
  cdiff ~ std_normal();
}
generated quantities {
  real loglik = 0;
  
  for (t in 1 : ntimes) {
    loglik += normal_lpdf(log(It[t]) | log(q * Bt[t]), tau);
  }
}
