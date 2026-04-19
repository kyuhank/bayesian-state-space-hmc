// SBC / centered / soft biomass floor + soft harvest-rate penalty
// Design intent:
//   - Use observed catch Ct directly (no PredCt reconstruction).
//   - Prevent negative/near-zero biomass with a differentiable soft floor.
//   - Encourage exploitation rate Ht <= Hmax via a soft penalty term.
// Notes:
//   - tau2Prior[3] toggles the index likelihood (1 = on, otherwise off).
data {
  // Time series length
  int ntimes;
  // Abundance index observations
  vector[ntimes] It;
  // Input catch series (used directly in the process update)
  vector[ntimes] Ct;

  // Prior hyperparameters: mean/sd on log scale for r and K
  vector[2] logrPrior;
  vector[2] logKPrior;
  // Prior hyperparameters: inv-gamma(shape, scale) for process/obs variances
  vector[2] sigma2Prior;
  vector[3] tau2Prior;

  // Uniform prior bounds for logq
  vector[2] logqPrior;
}

parameters {
  // Log-scale population parameters
  real logr;
  real logK;
  real<lower=logqPrior[1], upper=logqPrior[2]> logq;
  // Process and index variances
  real<lower=0> sigma2;
  real<lower=0> tau2;

  // Relative biomass states (centered parameterization)
  vector<lower=0>[ntimes+1] Pt;
}

transformed parameters {
  // Standard deviations
  real sigma = sqrt(sigma2);
  real tau = sqrt(tau2);

  // Natural-scale parameters
  real r = exp(logr);
  real K = exp(logK);
  real q = exp(logq);

  // Reference points
  real MSY = (r * K) / 4;
  real Hmsy = r / 2;
  real Bmsy = K / 2;

  // Derived monitoring quantities
  vector[ntimes] HtoverHmsy;
  vector[ntimes+1] BtoverBmsy;

  // Latent trajectory summaries
  vector<lower=0>[ntimes+1] PtMedian;
  vector<lower=0>[ntimes+1] Bt;
  vector[ntimes] Ht;
  // Predicted index trajectory
  vector[ntimes] PredIt;

  // Fixed tuning values
  real smoothness = 50;  // biomass soft-floor smoothness
  real Hmax = 0.999;  // harvest-rate cap
  real H_pen_scale = 0.01;  // penalty strength

  // Smoothing for the soft hinge on the harvest-rate constraint
  real<lower=0> H_pen_smoothness = 50;  // e.g., 20-50

  // Initial median relative biomass
  PtMedian[1] = 1;
  Bt[1] = Pt[1] * K;

  for (t in 2:(ntimes+1)) {
    // Logistic surplus production term
    real surplus_production = r * Pt[t-1] * (1.0 - Pt[t-1]);

    // Input catch Ct is used directly (not modified)
    real biomass_after_production = Pt[t-1] + surplus_production;
    real harvest_impact = Ct[t-1] / K;
    real raw_next_biomass = biomass_after_production - harvest_impact;

    // Differentiable lower bound for next-step biomass median
    PtMedian[t] = log1p_exp(smoothness * raw_next_biomass) / smoothness + 1e-4;

    // Biomass and exploitation trajectory
    Bt[t] = Pt[t] * K;

    // Harvest rate based on input catch
    Ht[t-1] = Ct[t-1] / Bt[t-1];
  }

  // Expected abundance index
  PredIt = q * Bt[1:ntimes];

  HtoverHmsy = Ht / Hmsy;
  BtoverBmsy = Bt / Bmsy;
}

model {
  // Priors
  logq ~ uniform(logqPrior[1], logqPrior[2]);
  logr ~ normal(logrPrior[1], logrPrior[2]);
  logK ~ normal(logKPrior[1], logKPrior[2]);

  sigma2 ~ inv_gamma(sigma2Prior[1], sigma2Prior[2]);
  tau2   ~ inv_gamma(tau2Prior[1], tau2Prior[2]);

  // Process likelihood: centered state model
  Pt ~ lognormal(log(PtMedian), sigma);

  // Index likelihood
  if (tau2Prior[3] == 1) {
    It ~ lognormal(log(q * Bt[1:ntimes]), tau);
  }

  // Soft inequality constraint: penalize only positive excess above Hmax
  for (t in 1:ntimes) {
    real z = Ht[t] - Hmax;
    real excess = log1p_exp(H_pen_smoothness * z) / H_pen_smoothness;
    target += normal_lpdf(excess | 0, H_pen_scale);
  }
}

generated quantities {
  // Pointwise log-likelihood sum for diagnostics/comparison
  real loglik = 0;
  for (t in 1:ntimes) {
    loglik += normal_lpdf(log(It[t]) | log(q * Bt[t]), tau);
  }
}
