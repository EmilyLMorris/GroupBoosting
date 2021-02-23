#' Group boosting fitting
#'
#' @param X matrix of predictors including any adjustment variables, which may not be grouped, of dimension n $/times$ p.
#' @param Y vector of the response variable, should be length n.
#' @param group vector containing an indicator of group membership of the p predictors in \code{X}
#' @param total_steps integer value of the maximum number of iterations to run the boosting algorithm
#' @param step_size a double value used to control how much to increase the estimate by at each step
#' @param adj_var integer number for the number of adjustment variables, 999 used to denote no adjustment variables
#' @param stop_tol double value used to determine what change in AIC/BIC/EBIC is small enough to stop iterating
#' @param gamma double value used to determine the penalty when using EBIC as stopping criteria
#' @param lasso_lambda double value of lambda used for lasso if regressing out the effect of adjustment variables
#' @param weighted default value is FALSE, set to TRUE if you want to use a weighted version to summarize within groups for the first stage
#' @return final_beta contains a vector with the final coefficient estimates, update1 is a vector that contains the coefficients from the first stage estimation, AIC1 contains the AIC path of the iterations from the first stage update, and BIC2 contains the BIC path from the second stage iterations
#' @export
group_boosting_fit <- function(X, Y, group, total_steps=50000, step_size=1e-6, adj_var = 999, stop_tol=-1e-7, gamma = 0, lasso_lambda = 0.0314, weighted = FALSE){
  if(!weighted){
    output = boosting_fitting_regress(X, Y, group, total_steps, step_size, adj_var, stop_tol, gamma, lasso_lambda)
    return(list(beta = output$final_beta, group_beta = output$update1, AIC = output$AIC1, BIC = output$BIC2, Z = output$Z))
  }else{
    output = boosting_fitting_weighted(X, Y, group, total_steps, step_size, adj_var, stop_tol, gamma, lasso_lambda)
    return(list(beta = output$final_beta, group_beta = output$update1, AIC = output$AIC1, BIC = output$BIC2, groups_weighted = output$groups.weighted))
  }
}

# weighted group boosting function
boosting_fitting_weighted <- function(X, Y, group, total_steps=50000, step_size=1e-6, adj_var = 999, stop_tol=-1e-7, gamma = 0, lasso_lambda = 0.0314){
  num.groups <- length(unique(group))
  groups.unique <- unique(group)
  group.size <- NULL
  sd.group <- NULL
  if(adj_var[1] != 999){
    fit.lasso1 = glmnet::glmnet(X, Y, penalty.factor=c(rep(0,length(adj_var)), rep(1, ncol(X) - length(adj_var))), lambda = lasso_lambda)
    fitted_values <- X[,adj_var]%*%fit.lasso1$beta[adj_var]
    residuals1 <- Y - fitted_values
    X.removed_adj = X[,-adj_var]
    index.select <- which(!(c(1:ncol(X)) %in% adj_var))
  }else{
    X.removed_adj = X
    index.select <- c(1:ncol(X))
    residuals1 = Y
  }
  X.group <- matrix(NA, nrow = nrow(X.removed_adj), ncol = num.groups*2)
  group.vec.new <- rep(NA, nrow(X.removed_adj))
  group.unique.new <- NULL
  for(g in 1:num.groups){
    group.temp <- groups.unique[g]
    group.idx <- which(group == group.temp)
    x.temp <- X.removed_adj[,which(group == group.temp)] # combine (sum or average) over groups
    if(is.vector(x.temp)){
      X.group[,(2*g-1)] = x.temp
      group.size = c(group.size, 1)
      group.vec.new[which(group == group.temp)] <- groups.unique[g]
      group.unique.new <- c(group.unique.new, groups.unique[g])
    }else{
      glmnet.temp <- glmnet::cv.glmnet(x.temp, Y, alpha = 1) # fit lasso to the edges in the group
      glmnet.temp.beta <- as.vector(stats::coef(glmnet.temp, s=glmnet.temp$lambda.min))[-1]
      group.neg <- glmnet.temp.beta[which(glmnet.temp.beta < 0)]
      group.pos <- glmnet.temp.beta[which(glmnet.temp.beta >= 0)]
      if(length(group.pos) > 1){
        strong.pos.signal <- max(group.pos)
        strong.signals.pos <- which.max(group.pos)
        weak.signals.pos <- which(group.pos < max(group.pos))
        if(length(weak.signals.pos) == 1){
          weighted.pos <- rowMeans(cbind(x.temp[,weak.signals.pos], x.temp[,strong.signals.pos]))
        }else{
          weighted.pos <- apply(x.temp[,weak.signals.pos], 2, function(x){rowMeans(cbind(x, x.temp[,strong.signals.pos]))})
        }
        x.temp.pos <- cbind(weighted.pos, x.temp[,strong.signals.pos])
        X.group[,(2*g-1)] = rowMeans(x.temp.pos)
        group.size = c(group.size, ncol(x.temp.pos))
        group.vec.new[group.idx[which(glmnet.temp.beta >= 0)]] <- paste(groups.unique[g], ":positive", sep = "")
        group.unique.new <- c(group.unique.new, paste(groups.unique[g], ":positive", sep = ""))
      }else if(length(group.pos) == 1){
        strong.pos.signal <- max(group.pos)
        strong.signals.pos <- which.max(group.pos)
        X.group[,(2*g-1)] <- x.temp[,strong.signals.pos]
        group.size = c(group.size, 1)
        group.vec.new[group.idx[which(glmnet.temp.beta >= 0)]] <- paste(groups.unique[g], ":positive", sep = "")
        group.unique.new <- c(group.unique.new, paste(groups.unique[g], ":positive", sep = ""))
      }else{
        X.group[,(2*g-1)] <- NA
      }
      if(length(group.neg) > 1){
        strong.neg.signal <- min(group.neg) # quantile(group.neg, 0.1)
        strong.signals.neg <- which.min(group.neg) # index of strong negative signal
        weak.signals.neg <- which(group.neg > min(group.neg))
        if(length(weak.signals.neg) == 1){
          weighted.neg <- rowMeans(cbind(x.temp[,weak.signals.neg], x.temp[,strong.signals.neg]))
        }else{
          weighted.neg <- apply(x.temp[,weak.signals.neg], 2, function(x){rowMeans(cbind(x, x.temp[,strong.signals.neg]))})
        }
        x.temp.neg <- cbind(weighted.neg, x.temp[,strong.signals.neg])
        X.group[,2*g] = rowMeans(x.temp.neg)
        group.size = c(group.size, ncol(x.temp.neg))
        group.vec.new[group.idx[which(glmnet.temp.beta < 0)]] <- paste(groups.unique[g], ":negative", sep = "")
        group.unique.new <- c(group.unique.new, paste(groups.unique[g], ":negative", sep = ""))
      }else if(length(group.neg) == 1){
        strong.neg.signal <- min(group.neg) # quantile(group.neg, 0.1)
        strong.signals.neg <- which.min(group.neg) # index of strong negative signal
        X.group[,2*g] <- x.temp[,strong.signals.neg]
        group.size = c(group.size, 1)
        group.vec.new[group.idx[which(glmnet.temp.beta < 0)]] <- paste(groups.unique[g], ":negative", sep = "")
        group.unique.new <- c(group.unique.new, paste(groups.unique[g], ":negative", sep = ""))
      }else{
        X.group[,2*g] <- NA
      }
    }
  }
  X.group <- X.group[,colSums(is.na(X.group)) < nrow(X.group)]
  X.std <- apply(X.group, 2, scale, center = TRUE, scale = TRUE)

  time.start <- proc.time()
  boosting.group.fit <- model_fitting(X.std, Y, x_size = group.size, total_steps, step_size, adj_var, stop_tol, gamma, stop_method = "AIC")
  time.groupsel = (proc.time() - time.start)[3]
  beta.update1 <- rep(NA, length = ncol(X))
  for(g in 1:num.groups){
    beta.update1[index.select[which(group == groups.unique[g])]] <- boosting.group.fit$beta[g]
  }
  if(adj_var[1] != 999){
    beta.update1[adj_var] <- 0
  }
  if(all(boosting.group.fit$beta == 0)){
    print("No groups selected.")
    X.red = X
  }else{
    groups.selected <- groups.unique[which(boosting.group.fit$beta != 0)]
    idx.keep <- which(group %in% groups.selected)
    if(adj_var[1] != 999){
      X.red <- as.matrix(X[,idx.keep])
    }else{
      X.red <- X[,idx.keep]
    }
  }
  step_factor = mean(group.size)
  boosting.ind.fit <- model_fitting(X.red, Y, x_size = rep(1, ncol(X.red)), total_steps*step_factor, step_size, adj_var, stop_tol, gamma, stop_method = "BIC")
  if(adj_var[1] != 999){
    beta.update2 <- rep(0, length = (ncol(X) - length(adj_var)))
    beta.update2[idx.keep] <- boosting.ind.fit$beta
  }else{
    beta.update2 <- rep(0, length = ncol(X))
    beta.update2[idx.keep] <- boosting.ind.fit$beta
  }
  return(list(final_beta = beta.update2, update1 = beta.update1, AIC1 = boosting.group.fit$AIC, BIC2 = boosting.ind.fit$BIC, groups.weighted = group.vec.new))
}

# function to fit group boosting algorithm and first regressing out the affect of adjustment variables with lasso
boosting_fitting_regress <- function(X, Y, group, total_steps=50000, step_size=1e-6, adj_var = 999, stop_tol=-1e-7, gamma = 0, lasso_lambda = 0.0314){
  num.groups <- length(unique(group))
  groups.unique <- unique(group)
  group.size <- NULL
  sd.group <- NULL
  if(adj_var[1] != 999){
    fit.lasso1 = glmnet::glmnet(X, Y, penalty.factor=c(rep(0,length(adj_var)), rep(1, ncol(X) - length(adj_var))), lambda = lasso_lambda)
    fitted_values <- X[,adj_var]%*%fit.lasso1$beta[adj_var]
    residuals1 <- Y - fitted_values
    X.removed_adj = X[,-adj_var]
    index.select <- which(!(c(1:ncol(X)) %in% adj_var))
  }else{
    X.removed_adj = X
    index.select <- c(1:ncol(X))
    residuals1 = Y
  }
  X.group <- matrix(NA, nrow = nrow(X.removed_adj), ncol = num.groups)
  for(g in 1:num.groups){
    group.temp <- groups.unique[g]
    x.temp <- X.removed_adj[,which(group == group.temp)]
    if(is.vector(x.temp)){
      X.group[,g] = x.temp
    }else{
      X.group[,g] = rowMeans(x.temp)
    }
    group.size = c(group.size, sum(group == groups.unique[g]))
  }
  adj_var1 = 999
  time.start <- proc.time()
  boosting.group.fit <- model_fitting(X.group, residuals1, x_size = group.size, total_steps, step_size, adj_var=adj_var1, stop_tol, gamma, stop_method = "AIC")
  time.groupsel = (proc.time() - time.start)[3] # takes about 20 minutes for real data with 1000 steps
  if(adj_var[1] != 999){
    beta.update1 <- rep(0, length = (ncol(X) - length(adj_var)))
    for(g in 1:num.groups){
      beta.update1[which(group == groups.unique[g])] <- boosting.group.fit$beta[g]
    }
  }else{
    beta.update1 <- rep(0, length = ncol(X))
    for(g in 1:num.groups){
      beta.update1[which(group == groups.unique[g])] <- boosting.group.fit$beta[g]
    }
  }
  if(all(boosting.group.fit$beta == 0)){
    print("No groups selected.")
    X.red = X
  }else{
    groups.selected <- groups.unique[which(boosting.group.fit$beta != 0)]
    idx.keep <- which(group %in% groups.selected)
    if(adj_var[1] != 999){
      X.red <- as.matrix(X.removed_adj[,idx.keep])
    }else{
      X.red <- X[,idx.keep]
    }
  }
  step_factor = mean(group.size)
  if(is.vector(X.red)){
    print("Only one variable selected")
    boosting.ind.fit <- model_fitting(X = X.red, y = residuals1, x_size = 1, total_steps = total_steps*step_factor, step_size = step_size, adj_var = adj_var1, stop_tol, gamma, stop_method = "BIC") # don't use residuals - use Y directly
  }else{
    boosting.ind.fit <- model_fitting(X = X.red, y = residuals1, x_size = rep(1, ncol(X.red)), total_steps = total_steps*step_factor, step_size = step_size, adj_var = adj_var1, stop_tol, gamma, stop_method = "BIC") # don't use residuals - use Y directly
  }
  if(adj_var[1] != 999){ # needed if adjustment variables
    beta.update2 <- rep(0, length = (ncol(X) - length(adj_var)))
    beta.update2[idx.keep] <- boosting.ind.fit$beta
  }else{
    beta.update2 <- rep(0, length = ncol(X))
    beta.update2[idx.keep] <- boosting.ind.fit$beta
  }
  return(list(final_beta = beta.update2, update1 = beta.update1, AIC1 = boosting.group.fit$AIC, BIC2 = boosting.ind.fit$BIC, Z = X.group))
}



# testing cpp function with example #### (put this in a different file)
# set.seed(2020)
# dat <- simul_dat_linear_model(n = 100, p = 10, rho_X = 0.5, effect_size = 0.5)
# step_size <- 0.01
# total_steps <- 5000
# stop_tol <- -1e-7
# beta_start <- rep(0, ncol(dat$X))
# adj_var = c(1,2,3)
# boosting_res = model_fitting(dat$X, dat$Y, rep(1,ncol(dat$X)), total_steps, step_size= 0.001, adj_var = 999, stop_tol, stop_method = "AIC")
# lasso_res = lasso_fitting(dat$X,dat$Y)
# enet_res = enet_fitting(dat$X,dat$Y,alpha=0.5)
# boosting_res$beta
# acctab = rbind(
#   compute_acc(boosting_res$beta,dat$beta),
#   compute_acc(lasso_res$beta,dat$beta),
#   compute_acc(enet_res$beta,dat$beta))
# rownames(acctab) = c("boosting no adj","lasso","enet")
# print(acctab,digit=2)
#
#
# compute_acc <- function(beta_est,true_beta){
#   true_idx = which(true_beta!=0)
#   select = factor(ifelse(beta_est!=0,1,0),levels=c(1,0))
#   true_select = factor(ifelse(true_beta!=0,1,0),levels=c(1,0))
#   tab = table(select,true_select)
#   TP = tab[1,1]
#   FP = tab[1,2]
#   FN = tab[2,1]
#   TN = tab[2,2]
#   Sensitivity  = TP/(TP+FN)
#   Specificity = TN/(TN+FP)
#   FDR = FP/(TP+FP)
#
#   return(c(null_mse = mean((beta_est[-true_idx]-true_beta[-true_idx])^2),
#            nonnull_mse = mean((beta_est[true_idx]-true_beta[true_idx])^2),
#            mse = mean((beta_est-true_beta)^2),
#            Sensitivity = Sensitivity,
#            Specificity = Specificity,
#            FDR = FDR))
# }

