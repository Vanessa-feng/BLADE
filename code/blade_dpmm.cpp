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
double dpois(int x, double lambda, bool log_value = true){
  if(lambda <= 0){
    throw std::logic_error("lambda should be positive");
  }
  double r;
  if(x == 0){
    r = - lambda;
  }else{
    vec s1 = arma::linspace(1, x, x);
    r = -sum(log(s1)) + x * log(lambda) - lambda;
  }
  if(log_value){return r;}
  else{return exp(r);}
}


// four parameter density function
vec dGBeta(vec x, double theta, double lambda,
           double alpha, double beta, bool log_value = true){
  if(lambda <= 0){
    throw std::logic_error("lambda should be positive");
  }
  int n = x.n_elem;
  vec output =vec(n, fill::zeros);
  double r = alpha / (alpha + beta);
  double term1 = lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta) - (alpha + beta - 1)*log(lambda);
  for(int i=0; i < n; i++ ){
    
    double term2 = (alpha - 1)*log(x(i) - theta + r * lambda);
    double term3 = (beta - 1)* log(-x(i) + theta + (1-r)*lambda);
    
    if(log_value){output(i) = term1 + term2 + term3;}
    else{output(i) =  exp(term1 + term2 + term3);}
  }
  return output;
}


// MH update theta
double update_theta(const double theta_old, const vec d, const double lambda, 
                    const double alpha, const double beta, const int z, 
                    const vec mu, const vec Sigma){
  double r = alpha / (alpha + beta);
  vec low = {min(d), max(d)-(1-r)*lambda};
  vec high = {min(d) + r*lambda, max(d)};
  double theta_prop = arma::randu(distr_param(low.max(), high.min()));
  double log_ratio = sum(dGBeta(d, theta_prop,lambda, alpha, beta) - dGBeta(d, theta_old,lambda, alpha, beta));
  log_ratio += log_normpdf(theta_prop, mu(z), sqrt(Sigma(z))) - log_normpdf(theta_old, mu(z), sqrt(Sigma(z)));
  //log_ratio += - 0.5 *(theta_prop - theta_old) * (theta_prop + theta_old - 2*mu(z)) / Sigma(z);
  double theta_new;  
  if(log_ratio > log(randu())){theta_new = theta_prop;}
  else{theta_new = theta_old;}
  return theta_new;
}


// update alpha and beta
vec update_r_total(double r_old, double total_old, vec d, double lambda, double theta){
  double r_min = (theta - min(d))/lambda;
  double r_max = 1 - (max(d) - theta) / lambda;
  
  double r_prop = randu(distr_param(r_min, r_max));
  double total_prop = randg(distr_param(1.0, total_old));
  double log_ratio = sum(dGBeta(d, theta,lambda, total_prop*r_prop, total_prop * (1-r_prop)) 
                           - dGBeta(d, theta,lambda, total_old*r_old, total_old * (1-r_old)));
  vec output = vec(2, fill::zeros);                          
  if(log_ratio > log(randu())){
    output(0) = r_prop;
    output(1) = total_prop;
  }
  else{
    output(0) = r_old;
    output(1) = total_old;
  }
  return output;
}


// update mu
vec update_mu(vec &mu_t, const vec group_t, const vec theta, const vec Sigma, const double mu0, const double tau){
  int K = mu_t.n_elem;
  for(int k=0; k < K; k++){
    uvec index = find(group_t == k);
    int n_k = index.n_elem;
    double mu_star = (sum(theta(index)) + mu0)/ (n_k + tau);
    double Sigma_star = Sigma(k) / (n_k + tau);
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


vec crp_weights(const vec &group_t, const double gamma){
  int N = group_t.n_elem; 
  arma::vec unique_clusters = unique(group_t);
  int K = unique_clusters.n_elem;
  
  // Compute histogram: counts of each cluster
  arma::Col<unsigned int> cluster_counts = hist(group_t, unique_clusters);
  arma::vec counts = arma::conv_to<arma::vec>::from(cluster_counts);
  
  // calculate probabilities
  arma::vec probs(K + 1);
  
  // Probability for existing clusters
  for (int k = 0; k < K; k++) {
    probs(k) = counts(k) / (gamma + N - 1);
  }
  
  // Probability for a new cluster
  probs(K) = gamma / (gamma + N - 1);
  return log(probs);
}

// integration
double integral(double x, double mu0, double tau, double alpha, double beta){
  double s = 1 + 1.0 / tau; 
  double term1 = -0.5 * (log(2*3.14159) +  log(s)) + alpha * log(beta) - 
    (alpha + 0.5)*log(beta + (x -mu0)*(x - mu0)/(2*s) );
  double term2 = lgamma(alpha + 0.5) - lgamma(alpha);
  return term1 + term2;
}

void update_group(vec &group_t, vec &mu_t, vec &Sigma_t, const vec &Y_t, 
                  const double &gamma, const double mu0, const double &tau, 
                  const double alpha, const double beta, const double &temperature){
  
  int N = Y_t.n_elem;
  vec clusters_name = unique(group_t);
  int K = clusters_name.n_elem;
  
  
  for(int i = 0; i < N; i++){
    int label_old = group_t(i);
    uvec group_i = find(group_t == label_old);
    int c_size = group_i.n_elem;
    
    // CPR prababilities
    vec probs = crp_weights(group_t, gamma);
    
    // prob of existing clusters
    for(int k = 0; k < K; k++){
      probs(k) += log_normpdf(Y_t(i), mu_t(k), sqrt(Sigma_t(k)));
    }
    
    // prob of new cluster
    probs(K) += integral(Y_t(i), mu0, tau, alpha, beta);
    
    // normalized
    probs = (probs - max(probs))/temperature;
    probs = exp(probs) / sum(exp(probs));
    
    // sample a new label
    int label_new = as_scalar(Rcpp::RcppArmadillo::sample(regspace(0, K), 1, false, probs));
  
    if(c_size > 1){
      if(label_new >= K){
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
      }
      group_t(i) = label_new;
      clusters_name = unique(group_t);
      K = clusters_name.n_elem;
    }else{
      //The cluster is a singleton, only K choices
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


// [[Rcpp::export]]
Rcpp::List runMCMC_dpmm(vec &group_t, vec &alpha_t, vec &beta_t, vec &theta_t, 
                   vec &mu_t, vec &Sigma_t, const vec &lambda, const Rcpp::List &X, 
                   const double &tau, const double &mu0, const double &alpha, const double &beta, 
                   const double gamma=1, const int max_iters = 100, const int seed = 12569){
  int N = group_t.n_elem;
  arma_rng::set_seed(seed);
  
  //store the parameters
  vec K_iter = zeros<vec>(max_iters);
  mat group_iter = zeros<mat>(N, max_iters);
  mat theta_iter = zeros<mat>(N, max_iters);
  mat alpha_iter = zeros<mat>(N, max_iters);
  mat beta_iter = zeros<mat>(N, max_iters);
  std::vector<vec> mu_history;
  std::vector<vec> Sigma_history;

  
  //store the initialized pars
  vec clusters = unique(group_t);
  int K = clusters.n_elem;
  K_iter(0) = K;
  group_iter.col(0) = group_t;
  theta_iter.col(0) = theta_t;
  alpha_iter.col(0) = alpha_t;
  beta_iter.col(0) = beta_t;
  mu_history.push_back(mu_t);
  Sigma_history.push_back(Sigma_t);

  for(int t = 1; t < max_iters; t++){
    double temperature = 3*pow(0.99, t);
    
    // update parameters for each cell
    for(int i=0; i < N; i++){
      theta_t(i) = update_theta(theta_t(i), X[i], lambda(i), alpha_t(i), beta_t(i), group_t(i), mu_t, Sigma_t);
      double total_old = alpha_t(i) + beta_t(i);
      double r_old = alpha_t(i) / total_old; 
      vec TEMP = update_r_total(r_old, total_old, X[i], lambda(i), theta_t(i));
      alpha_t(i) = TEMP(0) * TEMP(1);
      beta_t(i) = TEMP(1) - alpha_t(i);
    }
    mu_t = update_mu(mu_t, group_t, theta_t, Sigma_t, mu0, tau);
    Sigma_t = update_sigma(Sigma_t, group_t, theta_t, mu_t, tau, alpha, beta);
   
    update_group(group_t, mu_t, Sigma_t, theta_t, gamma, mu0, tau, alpha, beta, temperature);
    
    clusters = unique(group_t);
    K = clusters.n_elem;
    K_iter(t) = K;
    group_iter.col(t) = group_t;
    theta_iter.col(t) = theta_t;
    alpha_iter.col(t) = alpha_t;
    beta_iter.col(t) = beta_t;
    
    mu_history.push_back(mu_t);
    Sigma_history.push_back(Sigma_t);
  }
  
  return List::create(Named("K_iter") = K_iter,
                      Named("group_iter") = group_iter, 
                      Named("theta_iter") = theta_iter, 
                      Named("alpha_iter") = alpha_iter,
                      Named("beta_iter") = beta_iter,
                      Named("mu_iter") = mu_history, 
                      Named("Sigma_iter") = Sigma_history);
}







