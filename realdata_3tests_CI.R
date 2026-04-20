#3 tests and conditional independence assumption holds
library(ggplot2)
set.seed(0720)



rindependent_3_1 <- function(api,bpi,as1,bs1,ac1,bc1,as2,bs2,ac2,bc2,as3,bs3,ac3,bc3,result,n,burn){
  #api,bpi are the hyperparameters in the prior beta distribution for prevalence
  #as1,bs1 are the hyperparameters in the prior beta distribution for sensitivity of test1
  #ac1,bc1 are the hyperparameters in the prior beta distribution for specificity of test1
  #as2,bs2 are the hyperparameters in the prior beta distribution for sensitivity of test2
  #ac2,bc2 are the hyperparameters in the prior beta distribution for specificity of test2
  #as3,bs3 are the hyperparameters in the prior beta distribution for sensitivity of test3
  #ac3,bc3 are the hyperparameters in the prior beta distribution for specificity of test3
  #result is the observed data, which is a N*3 matrix
  #n is the number of cycles
  #burn is the number of cycles which are used to achieve convergence
  X <- matrix(0, n, 7)  
  pi0 <- rbeta(1,api,bpi)
  s10 <- rbeta(1,as1,bs1)
  c10 <- rbeta(1,ac1,bc1)
  s20 <- rbeta(1,as2,bs2)
  c20 <- rbeta(1,ac2,bc2)
  s30 <- rbeta(1,as3,bs3)
  c30 <- rbeta(1,ac3,bc3)
  X[1, ] <- c(pi0, s10, c10, s20, c20, s30, c30)  
  N <- dim(result)[1]
  D <- rep(0, N)
  for (i in 2:n){
    pi <- X[i-1, 1]  
    S1 <- X[i-1, 2]  
    C1 <- X[i-1, 3]  
    S2 <- X[i-1, 4]  
    C2 <- X[i-1, 5]  
    S3 <- X[i-1, 6]  
    C3 <- X[i-1, 7]  
    for (k in 1:N){
      p1 <- pi * (S1^(result[k,1])) * ((1-S1)^(1-result[k,1])) * 
        (S2^(result[k,2])) * ((1-S2)^(1-result[k,2])) * 
        (S3^(result[k,3])) * ((1-S3)^(1-result[k,3]))
      
      p2 <- (1-pi) * (C1^(1-result[k,1])) * ((1-C1)^(result[k,1])) * 
        (C2^(1-result[k,2])) * ((1-C2)^(result[k,2])) * 
        (C3^(1-result[k,3])) * ((1-C3)^(result[k,3]))
      
      pk <- p1/(p1+p2)
      D[k] <- rbinom(1,1,pk)
    }
    X[i, 1] <- rbeta(1, (sum(D)+api), (N-sum(D)+bpi))  
    X[i, 2] <- rbeta(1, (sum(D*result[,1])+as1), (sum(D*(1-result[,1]))+bs1))  
    X[i, 3] <- rbeta(1, (sum((1-D)*(1-result[,1]))+ac1), (sum((1-D)*result[,1])+bc1))  
    X[i, 4] <- rbeta(1, (sum(D*result[,2])+as2), (sum(D*(1-result[,2]))+bs2))  
    X[i, 5] <- rbeta(1, (sum((1-D)*(1-result[,2]))+ac2), (sum((1-D)*result[,2])+bc2))  
    X[i, 6] <- rbeta(1, (sum(D*result[,3])+as3), (sum(D*(1-result[,3]))+bs3))  
    X[i, 7] <- rbeta(1, (sum((1-D)*(1-result[,3]))+ac3), (sum((1-D)*result[,3])+bc3))  
  }
  b <- burn + 1
  x <- X[b:n, ]  
  return(x)
}



observed_1 <- matrix(rep(c(0, 0, 0), 1513), nrow = 1513, byrow = TRUE)
observed_2 <- matrix(rep(c(0, 0, 1), 21), nrow = 21, byrow = TRUE)
observed_3 <- matrix(rep(c(0, 1, 0), 59), nrow = 59, byrow = TRUE)
observed_4 <- matrix(rep(c(0, 1, 1), 11), nrow = 11, byrow = TRUE)
observed_5 <- matrix(rep(c(1, 0, 0), 23), nrow = 23, byrow = TRUE)
observed_6 <- matrix(rep(c(1, 0, 1), 19), nrow = 19, byrow = TRUE)
observed_7 <- matrix(rep(c(1, 1, 0), 12), nrow = 12, byrow = TRUE)
observed_8 <- matrix(rep(c(1, 1, 1), 34), nrow = 34, byrow = TRUE)
observed <- rbind(observed_1,observed_2,observed_3,observed_4,observed_5,observed_6,observed_7,observed_8)


result <- rindependent_3_1(1,1,
                           1,1,
                           1,1,
                           1,1,
                           1,1,
                           1,1,
                           1,1,
                           observed,12000,2000)

mean(result[,1])
mean(result[,2])
mean(result[,3])
mean(result[,4])
mean(result[,5])
mean(result[,6])
mean(result[,7])

