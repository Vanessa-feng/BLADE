#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rmath.h>
#include <iostream>

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends("RcppArmadillo")]]

// calculate the density of Poisson distribution
double dpois(int x, double tau, bool log_value = true){
  if(tau <= 0){
    throw std::logic_error("tau should be positive");
  }
  double r;
  if(x == 0){
    r = - tau;
  }else{
    vec s1 = arma::linspace(1, x, x);
    r = -sum(log(s1)) + x * log(tau) - tau;
  }
  if(log_value){return r;}
  else{return exp(r);}
}

// MH update theta
double compute_log_likelihood(const vec& d, double theta, double sigma_d=0.1) {
  double log_likelihood = 0.0;
  for (size_t i = 0; i < d.n_elem; i++) {
    log_likelihood += -0.5 * log(2 * M_PI * sigma_d * sigma_d)
    - (d(i) - theta) * (d(i) - theta) / (2 * sigma_d * sigma_d);
  }
  return log_likelihood;
}

double update_theta(const double theta_old, const vec d, const int z, const vec mu, 
                    const vec Sigma) {
  double theta_prop = arma::randn() * sqrt(Sigma(z)) + mu(z);
  
  double log_prior_old = log_normpdf(theta_old, mu(z), sqrt(Sigma(z)));
  double log_prior_new = log_normpdf(theta_prop, mu(z), sqrt(Sigma(z)));
  
  double log_likelihood_old = compute_log_likelihood(d, theta_old);
  double log_likelihood_new = compute_log_likelihood(d, theta_prop);
  
  double log_ratio = (log_likelihood_new + log_prior_new) - (log_likelihood_old + log_prior_old);
  
  double theta_new;
  if (log_ratio > log(randu())) {
    theta_new = theta_prop;
  } else {
    theta_new = theta_old;
  }
  return theta_new;
}

// update mu
vec update_mu(vec &mu_t, const vec group_t, const vec theta, const vec Sigma, const double mu0, const double p0){
  int K = mu_t.n_elem;
  for(int k=0; k < K; k++){
    uvec index = find(group_t == k);
    int n_k = index.n_elem;
    double mu_star = (sum(theta(index)) + mu0)/ (n_k + p0);
    double Sigma_star = Sigma(k) / (n_k + p0);
    mu_t(k) = randn(distr_param(mu_star, sqrt(Sigma_star)));
  }
  return mu_t;
}


// update Sigma
vec update_sigma(vec &Sigma_t, const vec group_t, const vec theta, const vec mu, 
                 const double tau, const double alpha, const double beta) {
  
  int K = Sigma_t.n_elem;
  for(int k=0; k < K; k++){
    uvec index = find(group_t == k);
    int n_k = index.n_elem;
    double rss = sum(pow(theta(index) - mu(k),2));
    
    double sigma_inv = randg(distr_param(n_k/2.0 + alpha, 1.0/(rss/2.0 + beta)));
    Sigma_t(k) = 1.0 / sigma_inv;
  }
  
  return Sigma_t;
}

//function to calculate factors for a new cluster
vec VN(const int Kmax, const int N, const double gamma, const double tau){
  vec vn = vec(Kmax, fill::zeros);
  int iters = max(N+100, 1000);
  for(int k = 1; k < Kmax +1 ; k++){
    double r = -datum::inf;
    for(int t = k; t < iters + 1; t++){
      double b = 0;
      vec s1 = arma::linspace(t-k+1, t, k);
      b +=  sum(log(s1));
      vec s2 = arma::linspace(t*gamma, t*gamma + N -1, N);
      b += - sum(log(s2));
      double s3 = dpois(t-1, tau);
      b +=  s3;
      double m = max(b, r);
      r = log(exp(r-m) + exp(b-m)) + m;
    }
    vn(k-1) = r;
  }
  return vn;
}


//function to calculate MRF energy function
double MRF_energy(const int &k, const int &l, const vec &group_t, const umat &G,
                  const double f, const bool log_value = true){
  
  uvec G_l = G.row(l).t();
  uvec neighbor = G_l.elem(find(G_l > 0)) - 1;
  
  uvec neighbor_id = find(group_t(neighbor) == k);
  double result = f * neighbor_id.n_elem;
  
  if(log_value){return result;}
  else{ return exp(result); }
}


//function to calculate probability of existing clusters
double prob_existing(const int k, const int l, const vec group_t, const vec mu_t, const vec Sigma_t,
                     const double y, const umat G, const double f, const double GAMMA){
  double result = 0;
  uvec index = find(group_t == k);
  int n_k;
  if(group_t(l) == k){
    n_k = index.n_elem - 1;
  }else{ n_k = index.n_elem;}
  
  result += log(n_k + GAMMA) + MRF_energy(k, l, group_t, G, f);
  
  result += log_normpdf(y, mu_t(k), sqrt(Sigma_t(k)));
  
  return result;
}

// integration
double integral(double x, double mu0, double tau, double alpha, double beta){
  double s = 1 + 1.0 / tau; 
  double term1 = -0.5 * (log(2*3.14159) +  log(s)) + alpha * log(beta) - 
    (alpha + 0.5)*log(beta + (x -mu0)*(x - mu0)/(2*s) );
  double term2 = lgamma(alpha + 0.5) - lgamma(alpha);
  return term1 + term2;
}
// function to calculate probability of a new cluster
double prob_new(const int K, const vec vn, const double y, double mu0, 
                double tau, double alpha, double beta, const double GAMMA){
  
  double h = integral(y, mu0, tau, alpha, beta);
  
  double result = max((h + log(GAMMA) + vn(K) - vn(K-1)), -999999.0);
  return result;
}

void update_group(vec &group_t, vec &mu_t, vec &Sigma_t, const vec &Y_t, 
                  const vec &vn, const umat &G, const double &f, const double &GAMMA, const double mu0, 
                  const double &tau, const double alpha, const double beta, const double &temperature){
  
  int Kmax = vn.n_elem - 1;
  int N = Y_t.n_elem;
  vec clusters_name = unique(group_t);
  int K = clusters_name.n_elem;
  
  for(int i = 0; i < N; i++){
    int label_old = group_t(i);
    int label_new;
    uvec group_i = find(group_t == label_old);
    int c_size = group_i.n_elem;
    
    
    if(c_size > 1){
      vec prob_i = zeros<vec>(K+1);
      for(int k = 0; k < K; k++){
        prob_i(k) = prob_existing(k, i, group_t, mu_t, Sigma_t, Y_t(i), G, f, GAMMA);
      }
      
      prob_i(K) = prob_new(K, vn, Y_t(i), mu0, tau, alpha, beta, GAMMA);
      
      
      prob_i = (prob_i - max(prob_i))/temperature;
      prob_i = exp(prob_i) / sum(exp(prob_i));
      
      label_new = as_scalar(Rcpp::RcppArmadillo::sample(regspace(0, K), 1, false, prob_i));
      
      if(label_new >= K && K < Kmax){
        
        //initial parameters for new cluster
        vec Sigma_t_new = zeros<vec>(K+1);
        Sigma_t_new.head(K) = Sigma_t;
        Sigma_t_new(K) = 1.0 / randg(distr_param(alpha, 1.0/beta));
        Sigma_t.swap(Sigma_t_new);
        
        vec mu_t_new = zeros<vec>(K+1);
        mu_t_new.head(K) = mu_t;
        double mu_star =  Y_t(i)/(1 + tau);
        double sd_star = sqrt(Sigma_t(K)/(1 + tau));
        mu_t_new(K) = randn(distr_param(mu_star, sd_star));
        mu_t.swap(mu_t_new);
        
        
        group_t(i) = label_new;
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
      }
      else if(label_new >= K && K >= Kmax){
        
        group_t(i) = label_old;
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
      }
      else{
        group_t(i) = label_new;
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
      }
      
    }
    
    else{
      
      //The cluster is a singleton, only K choices
      vec prob_i = vec(K+1, fill::zeros);
      for(int k = 0; k < K; k++){
        prob_i(k) = prob_existing(k, i, group_t, mu_t, Sigma_t, Y_t(i), G, f, GAMMA);
      }
      
      prob_i(K) = prob_new(K, vn, Y_t(i), mu0, tau, alpha, beta, GAMMA);
      
      prob_i = (prob_i - max(prob_i))/temperature;
      prob_i = exp(prob_i) / sum(exp(prob_i));
      
      label_new = as_scalar(Rcpp::RcppArmadillo::sample(regspace(0, K), 1, false, prob_i));
      
      if(label_new == label_old || label_new == K){
        group_t(i) = label_old;
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
      }
      else{
        uvec index = find(clusters_name != label_old);
        group_t(i) = label_new;
        for( int id = 0; id < N; id++){
          if(group_t(id) > label_old){
            group_t(id) += -1;
          }
        }
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
        
        if( K >= 1){
          
          mu_t = mu_t(index);
          Sigma_t = Sigma_t(index);
          
        }
      }
    }
  }
}

//add K_prior

// [[Rcpp::export]]
Rcpp::List runMCMC(vec &group_t, vec &theta_t, vec &mu_t, vec &Sigma_t, const Rcpp::List &X, 
                   const umat &G, const double &f, const double &tau, const double &mu0, 
                   const double &alpha, const double &beta, const int &K_prior, const double GAMMA=1,
                   const int max_iters = 100, const int seed = 12569){
  int N = group_t.n_elem;
  int Kmax = 50;  // Maximum number of clusters
  vec vn = VN(Kmax + 1, N, GAMMA, K_prior);  // Precompute VN terms
  arma_rng::set_seed(seed);
  
  // Store the parameters
  vec K_iter = zeros<vec>(max_iters);
  mat group_iter = zeros<mat>(N, max_iters);
  mat theta_iter = zeros<mat>(N, max_iters);
  mat mu_iter = zeros<mat>(Kmax, max_iters);
  mat Sigma_iter = zeros<mat>(Kmax, max_iters);
  
  // Initialize the parameters
  vec clusters = unique(group_t);
  int K = clusters.n_elem;
  K_iter(0) = K;
  group_iter.col(0) = group_t;
  theta_iter.col(0) = theta_t;
  mu_iter.col(0).head(K) = mu_t;
  Sigma_iter.col(0).head(K) = Sigma_t;
  
  for (int t = 1; t < max_iters; t++) {
    double temperature = 3 * pow(0.99, t);  // Annealing temperature
    
    // Update theta for each data point
    for (int i = 0; i < N; i++) {
      theta_t(i) = update_theta(theta_t(i), X[i], group_t(i), mu_t, Sigma_t);
    }
    
    mu_t = update_mu(mu_t, group_t, theta_t, Sigma_t, mu0, tau);
    
    Sigma_t = update_sigma(Sigma_t, group_t, theta_t, mu_t, tau, alpha, beta);
    
    update_group(group_t, mu_t, Sigma_t, theta_t, vn, G, f, GAMMA, mu0, 
                 tau, alpha, beta, temperature);
    // Store current iteration values
    clusters = unique(group_t);
    K = clusters.n_elem;
    K_iter(t) = K;
    group_iter.col(t) = group_t;
    theta_iter.col(t) = theta_t;
    mu_iter.col(t).head(K) = mu_t;
    Sigma_iter.col(t).head(K) = Sigma_t;
  }
  
  return List::create(Named("K_iter") = K_iter,
                      Named("group_iter") = group_iter,
                      Named("theta_iter") = theta_iter,
                      Named("mu_iter") = mu_iter,
                      Named("Sigma_iter") = Sigma_iter);
}







