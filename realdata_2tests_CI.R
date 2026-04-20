#2 tests, 1 population, conditional independence assumption holds
set.seed(0720)



rindependent_2_1 <- function(api,bpi,as1,bs1,ac1,bc1,as2,bs2,ac2,bc2,result,n,burn){
  #api,bpi are the hyperparameters in the prior beta distribution for prevalence
  #as1,bs1 are the hyperparameters in the prior beta distribution for sensitivity of test1
  #ac1,bc1 are the hyperparameters in the prior beta distribution for specificity of test1
  #as2,bs2 are the hyperparameters in the prior beta distribution for sensitivity of test2
  #ac2,bc2 are the hyperparameters in the prior beta distribution for specificity of test2
  #result is the observed data, which is a 2*2 matrix, the first row is n_11,n_10, the second row is n_01,n_00
  #n_10=N(T_1=1,T_2=0), n_01=N(T_1=0,T_2=1)
  #n is the number of cycles
  #burn is the number of cycles which are used to achieve convergence
  X <- matrix(0, n, 5)  
  pi0 <- rbeta(1,api,bpi)
  s10 <- rbeta(1,as1,bs1)
  c10 <- rbeta(1,ac1,bc1)
  s20 <- rbeta(1,as2,bs2)
  c20 <- rbeta(1,ac2,bc2)
  X[1, ] <- c(pi0, s10, c10, s20, c20)  
  n_11 <- result[1,1]  
  n_10 <- result[1,2]  
  n_01 <- result[2,1] 
  n_00 <- result[2,2]  
  N <- n_11+n_10+n_01+n_00  
  for (i in 2:n){
    pi <- X[i-1, 1]  
    S1 <- X[i-1, 2]  
    C1 <- X[i-1, 3]  
    S2 <- X[i-1, 4] 
    C2 <- X[i-1, 5]  
    y_11 <- rbinom(1, n_11, (pi*S1*S2)/((pi*S1*S2)+((1-pi)*(1-C1)*(1-C2))))  
    y_10 <- rbinom(1, n_10, (pi*S1*(1-S2))/((pi*S1*(1-S2))+((1-pi)*(1-C1)*C2))) 
    y_01 <- rbinom(1, n_01, (pi*(1-S1)*S2)/((pi*(1-S1)*S2)+((1-pi)*C1*(1-C2)))) 
    y_00 <- rbinom(1, n_00, (pi*(1-S1)*(1-S2))/((pi*(1-S1)*(1-S2))+((1-pi)*C1*C2))) 
    X[i, 1] <- rbeta(1, (y_11+y_10+y_01+y_00+api), (N-(y_11+y_10+y_01+y_00)+bpi)) 
    X[i, 2] <- rbeta(1, (y_11+y_10+as1), (y_01+y_00+bs1))  
    X[i, 3] <- rbeta(1, (n_01+n_00-y_01-y_00+ac1), (n_11+n_10-y_11-y_10+bc1))  
    X[i, 4] <- rbeta(1, (y_11+y_01+as2), (y_10+y_00+bs2)) 
    X[i, 5] <- rbeta(1, (n_10+n_00-y_10-y_00+ac2), (n_11+n_01-y_11-y_01+bc2)) 
    
  }
  b <- burn + 1
  x <- X[b:n, ]  
  return(x)
}


observed <- matrix(c(38,2,87,35),nrow=2,byrow=TRUE)
result <- rindependent_2_1(1,1,
                           4.44,13.31,
                           71.25,3.75,
                           21.96,5.49,
                           4.1,1.76,
                           observed,12000,2000)


mean(result[,1])
mean(result[,2])
mean(result[,3])
mean(result[,4])
mean(result[,5])






