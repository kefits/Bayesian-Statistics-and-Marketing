data {
  int<lower=0> D;
  int<lower=0> N;
  int<lower=0> M;
  int<lower=0, upper=1> y[N,D+1];
  matrix[N*D, M] X;
  
  matrix[M,M] A;
  real<lower=0> nu;
}
transformed data {
  int<lower=0> N_pos;
  int<lower=1,upper=N> n_pos[N];
  int<lower=1,upper=(D+1)> d_pos[N];
  int<lower=0> N_neg;
  int<lower=1,upper=N> n_neg[N*D];
  int<lower=1,upper=(D+1)> d_neg[N*D];
  
  vector[D] zero;
  matrix[D, D] I;
  zero = rep_vector(0, D);
  I = diag_matrix(rep_vector(1, D));
  
  
  N_pos = N;
  N_neg = N*D;
  {
    int i;
    int j;
    i = 1;
    j = 1;
    for (n in 1:N) {
      for (d in 1:(D+1)) {
        if (y[n,d] == 1) {
          n_pos[i] = n;
          d_pos[i] = d;
          i = i + 1;
        } else {
          n_neg[j] = n;
          d_neg[j] = d;
          j = j + 1;
        }
      }
    }
  }
}
parameters {
  vector[M] beta;
  cov_matrix[D] Sigma;
  vector<lower=0>[N] z_pos;
  vector<upper=0>[N*D] z_neg;
}
transformed parameters {
  vector[D+1] z[N];
  vector[D+1] temp[N];
  vector[D] w[N];
  vector[D] w_ast[N];
  vector[N*D] mu;
  matrix[D, D] L;
  matrix[M, M] R;
  
  for (n in 1:N_pos) {
    z[n_pos[n], d_pos[n]] = z_pos[n];
  }
  for (n in 1:N_neg) {
    z[n_neg[n], d_neg[n]] = z_neg[n];
  }
  for (n in 1:N) {
    temp[n] = z[n] - z[n][D+1];
    w_ast[n] = temp[n][1:D];
  }
  
  R = cholesky_decompose(A);
  L = cholesky_decompose(Sigma);
  mu = X*beta;
  for (n in 1:N){
    w[n] = w_ast[n] + mu[((n-1)*D+1):n*D];
  }
}
model {
  Sigma ~ inv_wishart(nu, nu*I);
  beta ~ multi_normal_cholesky(rep_vector(0,M), R);
  target += multi_normal_cholesky_lpdf(w_ast | zero, L);
}
generated quantities {
  vector[M] beta_tilde;
  cov_matrix[D] Sigma_tilde;
  
  beta_tilde = beta/sqrt(Sigma[1,1]);
  Sigma_tilde = Sigma/Sigma[1,1];
}