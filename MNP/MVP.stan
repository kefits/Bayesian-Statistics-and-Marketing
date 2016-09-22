data {
  int<lower=0> D;
  int<lower=0> N;
  int<lower=0, upper=1> y[N,D+1];
  matrix[N, D] X;
}
transformed data {
  int<lower=0> N_pos;
  int<lower=1,upper=N> n_pos[N];
  int<lower=1,upper=D> d_pos[N];
  int<lower=0> N_neg;
  int<lower=1,upper=N> n_neg[N*(D-1)];
  int<lower=1,upper=D> d_neg[N*(D-1)];
  N_pos = N;
  N_neg = N*D-N;
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
  real beta;
  cov_matrix[D] Sigma;
  vector<lower=0>[N] z_pos;
  vector<upper=0>[N*(D-1)] z_neg;
}
transformed parameters {
  vector[D+1] z[N];
  vector[D+1] temp[N];
  vector[D] w[N];
  
  for (n in 1:N_pos) {
    z[n_pos[n], d_pos[n]] = z_pos[n];
  }
  for (n in 1:N_neg) {
    z[n_neg[n], d_neg[n]] = z_neg[n];
  }
  for (n in 1:N) {
    temp[n] = z[n] - z[n][D+1];
    w[n] = temp[n][1:D];
  } 
}
model {
  int nu;
  matrix[D, D] I;
  nu = D+2;
  I = diag_matrix(rep_vector(1, D));
  
  Sigma ~ inv_wishart(nu, nu*I);
  beta ~ normal(0, 100);
  w ~ multi_normal(X*beta, Sigma);
}
generated quantities {
  real beta_tilde;
  cov_matrix[D] Sigma_tilde;
  
  beta_tilde = beta/Sigma[1,1];
  Sigma_tilde = Sigma/Sigma[1,1];
}