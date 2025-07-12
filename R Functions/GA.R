GA = function(K, # kinship matrix
              DatasetName, # which dataset
              n = c(50, 100), # training set size
              option = "Diff", # choose different training sets among trials
              n_iter = 12000, # number of minimum iterations
              index, # CD criterion
              new = T, # whether there's an ongoing GA process
              file_path = paste0("data/", DatasetName, "/", index, " ", "Number = ", paste(n, collapse = " & ") ) # the path for saving each iteration of GA
              ){  # n : Number of training set
  
  cat("\n'", "Dataset :", DatasetName, ";",
      "Training Set Size =", paste(n, collapse = " & "), ";", 
      "Index =", index, "'\n")
  
  ## 1. Initialization
  
  N = nrow(K) # number of all varieties
  cand = seq(N) # candidates
  env = length(n) # number of trials
  nt = mean(n)
  
  # number of training sets
  if (option == "Same" & length(unique(n)) == 1) { 
    list_len = 1 # same optimal training set across all trials
  } else if (option == "Diff") {
    list_len = seq(env) # different optimal training sets across all trials
  } else {
    stop("Invalid option. Please choose either 'Same' or 'Diff'.")
  }
  
  # Set Solution Size
  if (nt < 50){ ss = 50 } else if (nt >= 50 && nt <= 200){ ss = nt } else (ss = 200)
  
  
  if(new == T){ # start GA from the beginning
    
    # randomly sample optimal training set as a data frame
    sol = lapply(seq(ss), function(i) {
      sapply(list_len, function(j) {
        sample(c(rep(1, n[j]), rep(0, N - n[j]))) # 1 represents the varieties selected
      }) %>% as.data.frame()
    })
    
    max.score = data.frame(matrix(NA, 1, env + 1)) # the max CD value of each iteration
    colnames(max.score) = c("Overall", paste0("Trial", seq_len(env)))
    
  } else if (new == F){      # if GA was interrupted, start from where it was paused.
    load( paste0(file_path, ".RData") )
  }
  
  # the solutions and scores from the last iteration
  prev_sol = lapply(seq(ss), function(i) {
    sapply(list_len, function(j) {rep(0, N)}) %>% as.data.frame()
  })
  prev_score = matrix(0, nrow = ss, ncol = env+1)

  
  
  ## Progress Bar
  suppressPackageStartupMessages({
    library(httr)
    library(progress)
  })
  
  
  pb <- progress_bar$new(format = paste0(" [:bar] :current/", n_iter, " :percent ; Time: :elapsedfull ; ","Rate: :tick_rate iter/sec", "     :spin"), 
                         clear = FALSE, width = 100, total = (n_iter+100))
  pb$tick(iter)
  
  ## Parallel Computing
  library(parallel)
  
  myCoreNums = detectCores() # core numbers of the computer
  cl <- makeCluster(myCoreNums*0.6) # only use 60% of the cores
  invisible(clusterEvalQ(cl, suppressPackageStartupMessages( # don't show warning messages
    c(library(parallel), library(tidyverse), library(magrittr), 
      library("glmnet"), library("compiler", include.only = "cmpfun")))))
  clusterExport(cl, c("CD", "ss", "env", "cand", "K", "N", "n", "list_len", "index"), envir = environment()) # export items to every core
  
  
  ### Start GA iteration
  
  stop = F # criteria for GA to stop
  while(stop == F){
    iter = iter + 1 #
    
    # list out the solutions that are changed last iteration
    changed_sol = c()
    for(i in seq(ss)){
      if(!(all(sol[[i]] == prev_sol[[i]]))){
        changed_sol = c(changed_sol, i)
      }
    }
    
    clusterExport(cl, c("sol", "ss", "changed_sol", "prev_score", "max.score"), envir = environment()) # export the items that are changed every iteration
    
    ## 2. Evaluation: calculate CD value
    score = parLapply(cl, seq(ss), function(i) {
      
      # only calculate the training sets that are different from last iteration to reduce computing time
      if (i %in% changed_sol) { 
        train = lapply(seq(env), function(e) {
          cand[which(sol[[i]][, (e - 1) %% ncol(sol[[i]]) + 1] == 1)]
        })
        res_all = sapply(0:env, function(trial_id) { # calculate the overall CD value and that from each trial
          CD(K = K, train = train, index = index, trial = trial_id)
        })
        res_all
      } else {
        prev_score[i, ] # don't calculate the unchanged training sets 
      }
    }) %>% 
      do.call(rbind, .) %>% 
      as.data.frame()
    
    colnames(score) = c("Overall", paste0("Trial", 1:env) )
    
    # check if I didn't calculate the changed solutions
    check = c()
    for(i in seq(ss)){ if( all (prev_score[i,]==score[i,])  ){ check = c(check, i) } }
    if( all(setdiff(seq(ss),changed_sol) != check) & length(check)!=0 ){ stop("Didn't calculate the changed solutions.") } # pause GA

    
    prev_sol = sol ; prev_score = score # the training set and score from the previous iteration is replaced by that of this iteration
    
    max.score[iter, ] = score[which.max(score$Overall), ] # the max CD value of this iteration
    row.names(max.score)[iter] = iter
    
    elite = which(rank(-score$Overall, ties.method = "max") <= ceiling(ss * 0.1)) # the top 10% training sets as "elite"
    del = which(rank(-score$Overall, ties.method = "min") >= ceiling(ss * 0.6)) # the bottom 60% training sets as "delete"
    
    # print CD value
    opt.sol = which(score$Overall == max(score$Overall)) # the training set with the highest overall CD value
    opt.sol = sol[[ opt.sol[1] ]]
    opt = lapply(seq(env), function(e) {
      cand[which(opt.sol[, (e - 1) %% ncol(sol[[i]]) + 1] == 1)]
    })
    
    clusterExport(cl, c("opt"), envir = environment())
    
    ### Stop criteria
    
    # when the number of iterations exceed the number of minimum iterations
    th = 0.0001 # threshold
    if(iter > n_iter){ 
      if (all(sapply(env+1, function(e) { 
        # converge: stop GA when the difference between this iteration and the last 500 iteration didn't exceed the threshold
        score[max.row, e] - max.score[iter - 500, e] < th}
      ) )) { stop = T } 
    }
    
    if(iter > 2){
      if (all(sapply(1, function(e) { 
        score[max.row, e] < max.score[iter-1, e]} # stop GA when the max CD value decreases
      ) ))  {stop("Max score drops!")} 
    }
    
    
    
    ## 3. Crossover
    for(i in del){ # delete training sets
      repeat {
        chr = sample(seq(ss)[elite], 2) # sample two elite training sets
        pos = sample(seq(N - 1), 1) # sample the location where crossover occurs
        for(e in seq(list_len)){ # replace the delete training sets with combinations of two elite training sets
          sol[[i]][,e] = c(sol[[ chr[1] ]][1:pos, e], sol[[ chr[2] ]][ (pos+1):N, e])
          }
        
        
        # in case the outcome training set size is wrong during crossover
        n.sol = sapply(seq(env), function(e){ length(which(sol[[i]][, (e - 1) %% ncol(sol[[i]]) + 1] == 1)) })
        
        if (all(abs(n.sol-n) == 0)) { break } # don't stop when the training set size is wrong 
      }
    }
    
    ## 4. Mutation    
    sol.new = sol 
    r = 0.04*n # the number of varieties to be mutated
    
    for(i in seq(ss)){
      for(e in seq(list_len)){
        # sample two varieties which aren't selected (0) and two that are selected (1)
        pos = c( sample(which(sol[[i]][,e] == 0), ceiling(r[e])), sample(which(sol[[i]][,e] == 1), ceiling(r[e])) ) 
        
        # replace 0 into 1 and 1 into 0
        sol.new[[i]][pos, e] = abs(sol.new[[i]][pos, e] - 1) 
      }
    }
    
    clusterExport(cl, c("ss", "sol.new", "sol", "elite", "del", "score", "n"), envir = environment())   
    sol = parLapply(cl, seq(ss), function(i){
      
      if (!(i %in% elite)) {
        sol.new[[i]]  # directly output the mutated solution for non-elite training sets
        
      } else { # output the mutated training set for elite training sets only when the mutated training set has higher CD value
        old = score[i, ]
        train.new = lapply(seq(env), function(e) {
          cand[which(sol.new[[i]][, (e - 1) %% ncol(sol[[i]]) + 1] == 1)]
        })
        
        new = sapply(0:env, function(e) {
          CD(K = K, train = train.new, index = index, trial = e)
        })
        
        if ( all(new > old) ) {sol.new[[i]]} else {sol[[i]]} # output the original training set if new CD value isn't higher
        
      }
    })
    
    pb$tick()  # update progress bar
    
    # save the iteration number, the training set, and the max CD score of each iteration "in case GA is interrupted"
    save(iter, sol, max.score, file = paste0(file_path, ".RData")) 
    
    if(stop == T) {cat("\nGenetic Algorithm has ended!\n\n")} # stop GA when the stop criteria is met
  } ### End GA
  
  stopCluster(cl) # stop parallel computing
  if (!pb$finished) { pb$terminate() } # terminate progress bar
  
  ## Output optimized training set
  
  # the training set with the highest CD value in the last iteration is selected to be the optimized one
  opt.sol = which(score$Overall == max(score$Overall)) 
  opt.sol = sol[[ opt.sol[1] ]]
  opt.set = sapply(seq(env), function(e) {
    c(cand[which(opt.sol[, (e - 1) %% ncol(sol[[i]]) + 1] == 1)], rep(NA, max(n) - n[e]))
  }) %>% as.data.frame() #save the optimized training set into a data frame
  
  # retuen optimized training set as the first element, and the highest CD values of each iteration as the second element 
  return(list(Opt.set = opt.set, CD = max.score)) 
  
}
