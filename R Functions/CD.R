CD = function(K, train, index = "CD_mean(v2)",
              var_G = 10, # sigma_G^2
              var_ax1 = c( rep(10, length(train)-1), 20 ), # sigma_(ax1)^2
              var_E = NULL, # sigma_E^2
              trial = 0){ # calculate the CD value in which trials
  
  env = length(train) # number of trials
  
  # variance components of MGE model
  Omega_G = matrix(var_G, env, env) # Omega_G
  diag(Omega_G) = diag(Omega_G) + var_ax1 
  
  if( is.null(var_E) ){var_E = diag(Omega_G)} # diagonal elements of Omega_E

  nc = nrow(K) # the number of the whole candidate population  
  n = c() # the number of training set varieties
  for(e in seq(env)){n[e] = length(train[[e]])}

  # the matrices Gt, Gct and Mt
  Gt = Mt = matrix(0, nrow = sum(n), ncol = sum(n))
  Gct = matrix(0, nrow = nc*env, ncol = sum(n))
  
  # fill in the elements of the matrices Gt, Gct and Mt
  for(i in seq(env)){
    idx_i = sum(n[seq(i)]) - n[i] + seq(n[i])
    
    for(j in seq(env)){
      idx_j = sum(n[seq(j)]) - n[j] + seq(n[j])
      idx_all = nc*(i-1) + seq(nc)
      
      Gt[idx_i, idx_j] = Omega_G[i,j] * K[train[[i]], train[[j]]] # Gt
      Gct[idx_all, idx_j] = Omega_G[i,j] * K[, train[[j]]] # Gct
    }
    
    ni = n[i]
    Mt[idx_i, idx_i] = 1/var_E[i] * (diag(1, ni, ni) - matrix(1/ni, ni, ni)) # Mt
  }
  
  # matrices A and B
  A = Gct %*% solve(Mt %*% Gt + diag(1, sum(n), sum(n))) %*% Mt %*% t(Gct) # A: variance matrix of g_hat
  B = Omega_G %x% K # B: variance matrix of g
  
  if (trial == 0) { # all of the trials
    
    if(index == "CD_mean(v2)"){ # CD criterion is CD_mean(v2)
      
      CD_value = mean(diag(A) / diag(B)) # diagonal elements of A and B
      
    } else if (index == "CD_mean.MET"){ # CD criterion is CD_mean.MET
      
      A_all = B_all = 0
      
      # sum the corresponding element of the same variety in every sub-matrix
      for (i in seq(env)) {
        i_all = nc*(i-1) + seq(nc)
        for (j in seq(env)) {
          j_all = nc*(j-1) + seq(nc)
          
          A_sub = diag(A[i_all, j_all])
          B_sub = diag(B[i_all, j_all])
          
          A_all = A_all + A_sub
          B_all = B_all + B_sub
        }
      }
      
      CD_value = mean(A_all / B_all)
    }
    
  } else if (trial %in% seq(env)) { # only calculate the CD value of a single trial
    
    idx_trial = seq(nc) + nc*(trial-1)
    CD_value = mean( c(diag(A) / diag(B))[idx_trial] )
    
  }
  
  return(CD_value)
}
