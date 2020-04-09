data {
    int<lower=1> N_U;
    int<lower=1> N_L;
    int<lower=1> C;
    simplex[C] A_U[N_U];
    simplex[C] A_L[N_L];
    int<lower=1,upper=C> G_L[N_L];
    real<lower=0> alpha;
    real<lower=0> beta;
    real<lower=0> epsilon;
    real<lower=0> delta;
}
parameters{
    simplex[C] M[C];
    vector<lower=0>[C] theta;
    simplex[C] p;
    vector<lower=0>[C] gamma;
}
transformed parameters {
    vector<lower=0>[C] M_prior[C];
    for(i in 1:C) {
        for(j in 1:C) {
            M_prior[i,j] = gamma[i] * epsilon + gamma[i] * (i == j);
        }
    }
}
model {
    // Prior on p
    vector[C] log_p = log(p);
    p ~ dirichlet(rep_vector(delta, C));
    // Prior on theta
    theta ~ normal(0, 5);
    // Prior on M
    for(i in 1:C) {
        M[i] ~ dirichlet(M_prior[i]);
    }
    // Prior on gamma
    for(i in 1:C) {
      gamma[i] ~ gamma(alpha, beta);  
    }
    // Likelihood for data with observed uncertainty
    for(r in 1:N_L) {
        A_L[r,] ~ dirichlet(theta[G_L[r]] * M[G_L[r]]);
    }
    // Dirichlet mixture likelihood
    for(r in 1:N_U) {
        vector[C] lp = log_p;
        for(i in 1:C) {
            lp[i] += dirichlet_lpdf(A_U[r] | theta[i] * M[i]);
        }
        target += log_sum_exp(lp);
    }
}
