data{
  int N;
  vector[N] t;
  vector[N] y; 
  vector[N] trtm; // 1: Negative ctrl, 2: Positive ctrl, 3: Recovery treatment
}

parameters{
  real<lower=0,upper=1> sigma; // common sd
  real<lower=0,upper=2> k_dec; // time-constant for decrease
  real<lower=0,upper=2> k_inc; // time-constant for recovery
  real<lower=0> y_max; // upper 
  real<lower=0> y_min; // lower
  real<lower=0> amp_inc; // amplitude of the recovery
  real<lower=0> t_cp; // time at the changing point
}

model{
  for(i in 1:N){
    if(trtm[i] == 1){ // for NC 
      target += normal_lpdf(y[i] | y_max, sigma);
    } else if(trtm[i] == 2){ // for PC
      target += normal_lpdf(y[i] | y_min + (y_max - y_min) * exp(-k_dec * t[i]), sigma);
    } else { // for RT
      if(t[i] <= t_cp){ // before the changing point
        target += normal_lpdf(y[i] | y_min + (y_max - y_min) * exp(-k_dec * t[i]), sigma);
      } else { // before the changing point
        target += normal_lpdf(y[i] | y_min + (y_max - y_min) * exp(-k_dec * t[i]) + amp_inc * (1 + exp(k_inc * (t_cp - t[i]))), sigma);
      }
    }
  }
}

generated quantities{
  vector[16] PC;
  vector[16] NC;
  vector[16] RT;
  
  for(i in 1:16){
    NC[i] = normal_rng(y_max, sigma);
    PC[i] = normal_rng(y_min + (y_max - y_min) * exp(-k_dec * (i - 1)), sigma);
    RT[i] = normal_rng(y_min + (y_max - y_min) * exp(-k_dec * (i - 1)), sigma);
    if((i - 1) > t_cp)
      RT[i] += amp_inc * (1 + exp(k_inc * (t_cp - (i - 1))));
  }
}
