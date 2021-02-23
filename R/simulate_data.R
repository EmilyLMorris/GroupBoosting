omega.pd.func2 <- function(x, nrow){ # input x is the current omega for one subject (one row)
  x.mat <- matrix(NA, nrow = nrow, ncol = nrow, byrow = T)
  x.mat[lower.tri(x.mat)] <- x
  x.mat[upper.tri(x.mat)] = t(x.mat)[upper.tri(x.mat)]
  diag(x.mat) <- 1
  lambda_R <- min(eigen(x.mat)$values)
  epsilon <- stats::runif(1, 0, 0.001)
  x.pd <- x.mat + (abs(lambda_R) + epsilon)*diag(nrow)
  return(x.pd)
  # return(x.mat)
}

#' Simulate grouped data
#'
#' @param nodes integer value of the number of nodes in the network
#' @param n integer, number of subjects
#' @param num.groups integer, number of groups that the nodes in the network belong to
#' @param q.groups integer, number of group with non-zero signal
#' @param sparse_g ineger, number of groups with sparse signal
#' @param dense_g integer, number of groups with dense signal among the edges
#' @param effect_size double value indicating the effect size for the non-zero signals
#' @param SNR double value for the desired signal to noise ratio
#' @param structure default is "block" which defines the grouping among edges in the network, "nodes" is also available if the grouping should be by nodes instead
#' @param adj_var integer for the number of adjustment variables
#' @param rho_adj_var double value of the average correlation between connectivity of the edges and the adjustment variables
#' @param seed integer value to set.seed
#' @return Y contains the simulated response vector, X is a matrix of predictors including any adjustment variables, and group is a vector of the group membership of each edg
#' @export
simul_dat_group <- function(nodes = 25, n = 100, num.groups = 5 ,
                            q.groups = 3, sparse_g = 0.1,  dense_g = 0.9, # density of signal among
                            effect_size = 1, SNR = 0.9, structure = "block", adj_var = 999, rho_adj_var, seed = 2020){
  set.seed(seed)
  num.edges <- choose(nodes,2)
  components <- sample(1:2, prob=c(0.7,0.3), size=num.edges*n, replace=TRUE)
  mus <- c(1, -1)
  sds <- c(.1, .1)
  samples <- stats::rnorm(n=num.edges*n,mean=mus[components],sd=sds[components]) # use indicator to ensure large magnitude for omega nonzero values
  omega.sim <- matrix(samples, nrow = n, ncol = num.edges)

  if(structure == "nodes"){
    p = nodes^2
    omega.pd <- t(apply(omega.sim, 1, omega.pd.func2, nodes)) # ensure positive definite
    group.mat <- matrix(rep(c(1:nodes), each = nodes), nrow = nodes, ncol = nodes, byrow = T)
    group.vec <- as.factor(as.vector(group.mat))
    group.sel <- c(1:ceiling(q.groups/2))
    group.sel.sparse <- c((ceiling(q.groups/2)+1):q.groups)

  }else{ # block group structure
    p = num.edges
    mat.ind <- matrix(NA, ncol = nodes, nrow = nodes)
    omega.pd <- omega.sim

    group.mat <- matrix(NA, nrow = nodes, ncol = nodes)
    row.ind <- rep(1:num.groups, each = nodes/num.groups)
    col.ind <- rep(1:num.groups, each = nodes/num.groups)
    for(ind in 1:length(row.ind)){
      for(ind2 in 1:length(col.ind)){
        group.mat[ind, ind2] <- paste(row.ind[ind], col.ind[ind2], sep =":")
      }
    }
    group.mat[lower.tri(group.mat)]  <- t(group.mat)[lower.tri(group.mat)]
    group.vec <- as.factor(as.vector(group.mat[lower.tri(group.mat)]))
    group.sel <- paste0(c(1:q.groups),":", c(1:q.groups)) # all on diagonal groups
    group.sel.sparse <- as.list(outer(c(1:q.groups), c(1:q.groups), paste, sep = ":")) # off diagonal groups
    group.sel.sparse <- unlist(group.sel.sparse[!group.sel.sparse %in% group.sel])
  }
  # simulate true beta
  true_beta <- rep(0, p)
  for(b in 1:length(group.sel)){ # assign magnitude to dense groups
    signal <- sample(which(group.vec == group.sel[b]), size = floor(dense_g*length(which(group.vec == group.sel[b]))))
    true_beta[signal] <- effect_size
  }
  for(b2 in 1:length(group.sel.sparse)){ # assign magnitude to sparse groups
    signal2 <- sample(which(group.vec == group.sel.sparse[b2]), size = ceiling(sparse_g*length(which(group.vec == group.sel.sparse[b2]))))
    true_beta[signal2] <- effect_size
  }
  sigma_eps = sqrt(stats::var(omega.pd%*%true_beta)*(1-SNR)/SNR )
  if(adj_var != 999){
    beta_adj_var <- stats::rnorm(adj_var, effect_size, effect_size/2)
    X.adj_var <- matrix(stats::rnorm(adj_var*n, 0, 1), nrow = n, ncol = adj_var)

    X = cbind(X.adj_var, omega.pd)
    beta_final = c(beta_adj_var, true_beta)
    sigma_eps = sqrt(stats::var(X%*%beta_final)*(1-SNR)/SNR )
    Y = X.adj_var%*%beta_adj_var + omega.pd%*%true_beta + stats::rnorm(n, mean = 0, sd = sigma_eps)
    SNR = stats::var(cbind(X.adj_var,omega.pd)%*%c(beta_adj_var, true_beta))/stats::var(Y)
  }else{
    Y = omega.pd%*%true_beta + stats::rnorm(n, mean = 0, sd = sigma_eps)
    SNR = stats::var(omega.pd%*%true_beta)/stats::var(Y)
    X = omega.pd
    beta_final = true_beta
  }
  return(list(Y=Y, X=X, beta=beta_final, group = group.vec))
}
