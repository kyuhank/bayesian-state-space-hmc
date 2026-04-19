// SBC / centered / hard-if catch constraint
// Constraint type:
//   PredCt[t] = Ct[t] unless Ct exceeds feasible harvest, then clipped to 0.99 * K * PtMedTemp
// Notes:
//   - This branch-based constraint is intentionally non-smooth at the boundary.
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
  
  vector<lower=0>[ntimes + 1] Pt;
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
  
  vector[ntimes + 1] PtMedian;
  vector[ntimes + 1] Bt;
  vector[ntimes] Ht;
  vector[ntimes] PredIt;
  vector[ntimes] PredCt;
  
  real<lower=0> PtMedTemp;
  
  real cdiff_scale = 0.01;
  
  // Initial condition on relative biomass median
  PtMedian[1] = 1;
  Bt[1] = Pt[1] * K;
  
  for (t in 2 : (ntimes + 1)) {
    // Deterministic biomass median before catch removal
    PtMedTemp = Pt[t - 1] + r * (1.0 - Pt[t - 1]) * Pt[t - 1];
    
    // Hard branch constraint on effective catch (non-differentiable boundary)
    if (PtMedTemp <= (Ct[t - 1] / K)) {
      PredCt[t - 1] = K * PtMedTemp * 0.99;
    } else {
      PredCt[t - 1] = Ct[t - 1];
    }
    
    PtMedian[t] = Pt[t - 1] + r * (1.0 - Pt[t - 1]) * Pt[t - 1]
                  - PredCt[t - 1] / K;
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
  
  // Process model
  Pt ~ lognormal(log(PtMedian), sigma);
  
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
