# include <RcppParallel.h>
# include <RcppArmadillo.h>
# include <sys/time.h>
# include <vector>
# include <algorithm>
# include <cstdlib>
# include <iostream>
#include <RcppArmadilloExtensions/sample.h>
#include <iomanip>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppParallel, RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]

struct boosting_linear_model{
  arma::mat X_all;
  arma::mat X;
  arma::vec y;
  arma::vec beta;
  int num_parameters;
  int total_steps;
  double step_size;
  arma::vec resid;
  arma::vec xj_resid;
  arma::vec xk_resid;
  int j_star;
  int p;
  int n;
  double AIC;
  double BIC;
  double ProfLik;
  double sigma2;
  double gamma;
  arma::cube hat_mat_list;
  arma::mat I_B_mat;
  double df;
  arma::mat I_n;
  int selected_idx;
  arma::vec x_size;
  arma::uvec adj_var;
  arma::uvec X_loc;
  int p_adj;
  arma::mat X_adj;
  arma::vec index_select;
  arma::mat P_A;
};

// [[Rcpp::export]]
arma::cube compute_hat_mat_list(arma::mat& X){
  arma::cube hat_mat_list = zeros<arma::cube>(X.n_rows,X.n_rows,X.n_cols);
  for(int j=0; j<X.n_cols;j++){
    hat_mat_list.slice(j) = X.col(j)*X.col(j).t()/norm(X.col(j));
  }
  return hat_mat_list;
}

void update_num_parameters(boosting_linear_model& model){
  arma::uvec nonzero_beta = find(model.beta != 0);
  model.num_parameters = nonzero_beta.size();
}

void udpate_I_B_mat(boosting_linear_model& model){
  model.I_B_mat *= model.I_n - model.step_size*model.hat_mat_list.slice(model.selected_idx);
}

void update_df(boosting_linear_model& model){
  model.df = trace(model.I_n - model.I_B_mat);
}

void update_sigma2(boosting_linear_model& model){
    model.sigma2 = arma::norm(model.resid);
}

void update_AIC(boosting_linear_model& model){
  model.AIC = log(model.sigma2) + (1.0 + model.df/model.n)/(1.0-(model.df+2.0)/model.n);
}

void update_BIC(boosting_linear_model& model){
  model.BIC = -(-2*log(model.sigma2) + model.num_parameters*log(model.n) - 2*model.gamma*lgamma(model.num_parameters + 1) + 2*model.gamma*lgamma(model.X.n_rows - model.num_parameters + 1));
}

void update_ProfLik(boosting_linear_model& model){
  arma::mat temp = model.resid.t()*model.P_A*model.resid;
  model.ProfLik = -(model.n*log(model.sigma2) + (1/model.sigma2)*temp(0,0));
}

void update_resid(boosting_linear_model& model){
  model.resid = model.y - model.X*model.beta;
}

void initialize_model(boosting_linear_model& model, arma::mat& X_all, arma::vec& y, arma::vec& x_size,
                      int& total_steps, double& step_size, arma::uvec& adj_var, double& stop_tol, double& gamma, Rcpp::String stop_method){
  model.n = X_all.n_rows;
  model.X_all = X_all;
  model.y = y;
  model.x_size = x_size;
  model.num_parameters = 0;
  model.total_steps = total_steps;
  model.step_size = step_size;
  model.I_n = eye<mat>(model.n,model.n);
  model.p_adj = adj_var.n_elem;
  if(adj_var(0) != 999){
    model.adj_var =  adj_var - 1;
    model.X_adj = X_all.cols(model.adj_var);
    model.P_A = model.I_n - model.X_adj*inv(model.X_adj.t()*model.X_adj)*model.X_adj.t() ; // 1 - hat matrix
  }else{
    model.adj_var = adj_var;
    model.p_adj = 0;
    model.P_A = model.I_n ;
  }
  model.p = X_all.n_cols - model.p_adj;
  model.beta = zeros<vec>(model.p);
  model.X = model.X_all.cols(model.p_adj, model.X_all.n_cols - 1);
  model.hat_mat_list = compute_hat_mat_list(model.X);
  model.I_B_mat = eye<mat>(model.n, model.n);
  model.xj_resid = zeros<vec>(model.p);
  model.selected_idx = -1;
  model.gamma = gamma;
  update_resid(model);
  update_sigma2(model);
  update_num_parameters(model);
  update_AIC(model);
  update_BIC(model);
  update_ProfLik(model);
}

void update_xj_resid(boosting_linear_model& model){
  if(model.p_adj != 0){
    model.xj_resid = model.X.t()*model.P_A*(model.y - model.X*model.beta);
  }else{
    for(int j=0;j<model.p;j++)
      model.xj_resid(j) = sum(model.X.col(j) % model.resid);
  }
}

void update_beta(boosting_linear_model& model){
  model.selected_idx = index_max(model.xj_resid);
  if(model.p_adj != 0){
    mat L2 = 2*model.X.t()*model.P_A*model.X;
    double L2_inv = 1/L2(model.selected_idx, model.selected_idx);
    model.beta(model.selected_idx) += model.step_size*L2_inv*(1/model.n)*model.xj_resid(model.selected_idx)*model.x_size(model.selected_idx);
  }else{
    model.beta(model.selected_idx) += model.step_size*model.xj_resid(model.selected_idx)*model.x_size(model.selected_idx);
  }
}

// [[Rcpp::export]]
List model_fitting(arma::mat& X,
                   arma::vec& y,
                   arma::vec& x_size,
                   int& total_steps,
                   double& step_size, arma::uvec& adj_var, double stop_tol = 1e-5, double gamma = 0, Rcpp::String stop_method = "AIC"){

  boosting_linear_model model;
  initialize_model(model, X, y, x_size, total_steps, step_size, adj_var, stop_tol, gamma, stop_method);
  vec AIC_list = zeros<vec>(total_steps);
  vec BIC_list = zeros<vec>(total_steps);
  vec ProfLik_list = zeros<vec>(total_steps);
  vec beta = model.beta;
  int step = 0;
  for(step=0;step<total_steps;step++){
    update_xj_resid(model);
    update_beta(model);
    update_resid(model);
    update_sigma2(model);
    if(stop_method == "BIC"){
      update_BIC(model);
      BIC_list(step) = model.BIC;
      if(step>100){
        if(BIC_list(step)-BIC_list(step-10)>stop_tol){
          break;
        } else{
          beta = model.beta;
          update_ProfLik(model);
          ProfLik_list(step) = model.ProfLik;
        }
      }else{
        beta = model.beta;
        update_ProfLik(model);
        ProfLik_list(step) = model.ProfLik;
      }
    }else{
      udpate_I_B_mat(model);
      update_df(model);
      update_AIC(model);
      AIC_list(step) = model.AIC;
      if(step>100){
        if(AIC_list(step)-AIC_list(step-10)>stop_tol){
          break;
        } else{
          beta = model.beta;
          update_ProfLik(model);
          ProfLik_list(step) = model.ProfLik;
        }
      }else{
        beta = model.beta;
        update_ProfLik(model);
        ProfLik_list(step) = model.ProfLik;
      }
    }
  }
  vec ProfLik_out = ProfLik_list.elem(linspace<uvec>(0,step-1,step));
  if(stop_method == "BIC"){
    vec BIC_out = BIC_list.elem(linspace<uvec>(0,step-1,step));
    return List::create(Named("beta") = beta,
                        Named("BIC") = BIC_out,
                        Named("steps") = step,
                        Named("ProfLik") = ProfLik_out);
  }else{
    vec AIC_out = AIC_list.elem(linspace<uvec>(0,step-1,step));
    return List::create(Named("beta") = beta,
                        Named("AIC") = AIC_out,
                        Named("steps") = step,
                        Named("ProfLik") = ProfLik_out);
  }
}

