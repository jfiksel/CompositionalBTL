functions {
  matrix vector_array_to_matrix(vector[] x) {
    matrix[size(x), rows(x[1])] y;
    for (m in 1:size(x))
      y[m] = x[m]';
    return y;
  }
  
}
data {
    int<lower=1> N;
    int<lower=1> C;
    simplex[C] A[N];
    matrix<lower=0,upper=1> [N,C] G;
    real<lower=0> alpha;
    real<lower=0> beta;
    real<lower=0> epsilon;
}
parameters{
    simplex[C] M[C];
    real<lower=0> theta;
    vector<lower=0>[C] gamma;
}
transformed parameters {
    vector<lower=0>[C] M_prior[C];
    matrix<lower=0,upper=1> [N,C] mu;
    mu = G * vector_array_to_matrix(M);
    for(i in 1:C) {
        for(j in 1:C) {
            M_prior[i,j] = gamma[i] * epsilon + gamma[i] * (i == j);
        }
    }
}
model {
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
    // Likelihood for data with observed classes
    for(r in 1:N) {
        A[r] ~ dirichlet(to_vector(theta * mu[r]));
    }
}
