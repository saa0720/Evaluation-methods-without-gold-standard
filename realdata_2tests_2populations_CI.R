#2 tests and 2 populations, conditional independence assumption and equal accuracy assumption holds
library(ggplot2)
set.seed(0720)

rindependent_2_2 <- function(api1,bpi1,api2,bpi2,as1,bs1,ac1,bc1,as2,bs2,ac2,bc2,result1,result2,n,burn){
  #api1,bpi1 are the hyperparameters in the prior beta distribution for prevalence
  #api1,bpi1 are the hyperparameters in the prior beta distribution for prevalence
  #as1,bs1 are the hyperparameters in the prior beta distribution for sensitivity of test1
  #ac1,bc1 are the hyperparameters in the prior beta distribution for specificity of test1
  #as2,bs2 are the hyperparameters in the prior beta distribution for sensitivity of test2
  #ac2,bc2 are the hyperparameters in the prior beta distribution for specificity of test2
  #result1 is the observed data in the first population, which is a 2*2 matrix
  #result2 is the observed data in the second population, which is a 2*2 matrix
  #n is the number of cycles
  #burn is the number of cycles which are used to achieve convergence
  X <- matrix(0, n, 6) 
  pi10 <- rbeta(1,api1,bpi1)
  pi20 <- rbeta(1,api2,bpi2)
  s10 <- rbeta(1,as1,bs1)
  c10 <- rbeta(1,ac1,bc1)
  s20 <- rbeta(1,as2,bs2)
  c20 <- rbeta(1,ac2,bc2)
  X[1, ] <- c(pi10, pi20, s10, c10, s20, c20)  
  n_111 <- result1[1,1]  
  n_101 <- result1[1,2]  
  n_011 <- result1[2,1] 
  n_001 <- result1[2,2] 
  N1 <- n_111+n_101+n_011+n_001  
  n_112 <- result2[1,1]  
  n_102 <- result2[1,2]  
  n_012 <- result2[2,1]  
  n_002 <- result2[2,2]  
  N2 <- n_112+n_102+n_012+n_002  
  for (i in 2:n){
    pi1 <- X[i-1, 1]  
    pi2 <- X[i-1, 2] 
    S1 <- X[i-1, 3]  
    C1 <- X[i-1, 4]  
    S2 <- X[i-1, 5]  
    C2 <- X[i-1, 6]  
    y_111 <- rbinom(1, n_111, (pi1*S1*S2)/((pi1*S1*S2)+((1-pi1)*(1-C1)*(1-C2))))  
    y_101 <- rbinom(1, n_101, (pi1*S1*(1-S2))/((pi1*S1*(1-S2))+((1-pi1)*(1-C1)*C2)))  
    y_011 <- rbinom(1, n_011, (pi1*(1-S1)*S2)/((pi1*(1-S1)*S2)+((1-pi1)*C1*(1-C2))))  
    y_001 <- rbinom(1, n_001, (pi1*(1-S1)*(1-S2))/((pi1*(1-S1)*(1-S2))+((1-pi1)*C1*C2)))  
    y_112 <- rbinom(1, n_112, (pi2*S1*S2)/((pi2*S1*S2)+((1-pi2)*(1-C1)*(1-C2))))  
    y_102 <- rbinom(1, n_102, (pi2*S1*(1-S2))/((pi2*S1*(1-S2))+((1-pi2)*(1-C1)*C2)))  
    y_012 <- rbinom(1, n_012, (pi2*(1-S1)*S2)/((pi2*(1-S1)*S2)+((1-pi2)*C1*(1-C2))))  
    y_002 <- rbinom(1, n_002, (pi2*(1-S1)*(1-S2))/((pi2*(1-S1)*(1-S2))+((1-pi2)*C1*C2)))  
    X[i, 1] <- rbeta(1, (y_111+y_101+y_011+y_001+api1), (N1-(y_111+y_101+y_011+y_001)+bpi1))  
    X[i, 2] <- rbeta(1, (y_112+y_102+y_012+y_002+api2), (N2-(y_112+y_102+y_012+y_002)+bpi2)) 
    X[i, 3] <- rbeta(1, (y_111+y_101+y_112+y_102+as1), (y_011+y_001+y_012+y_002+bs1)) 
    X[i, 4] <- rbeta(1, (n_011+n_001-y_011-y_001+n_012+n_002-y_012-y_002+ac1), (n_111+n_101-y_111-y_101+n_112+n_102-y_112-y_102+bc1))  
    X[i, 5] <- rbeta(1, (y_111+y_011+y_112+y_012+as2), (y_101+y_001+y_102+y_002+bs2))  
    X[i, 6] <- rbeta(1, (n_101+n_001-y_101-y_001+n_102+n_002-y_102-y_002+ac2), (n_111+n_011-y_111-y_011+n_112+n_012-y_112-y_012+bc2)) 
  }
  b <- burn + 1
  x <- X[b:n, ]  
  return(x)
}


observed1 <- matrix(c(14,4,9,528),nrow=2,byrow=TRUE)
observed2 <- matrix(c(887,31,37,367),nrow=2,byrow=TRUE)


result <- rindependent_2_2(1,1,
                           1,1,
                           1,1,
                           1,1,
                           1,1,
                           1,1,
                           observed1,observed2,
                           12000,2000)


mean(result[,1])
mean(result[,2])
mean(result[,3])
mean(result[,4])
mean(result[,5])
mean(result[,6])






