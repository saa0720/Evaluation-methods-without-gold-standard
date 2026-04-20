#random effect model for two tests with equal b assumption
#To utilize the prior given in the original text, we assume here thatP(T_{jk}=0|D_k=0,I_k=i_k)=\Phi(a_{j0}+b_0i_k)

set.seed(0720)


#Calculate the hyperparameter in Bernoulli distribution when generate true disease status
truedisease <- function(para, observed, measure){
  #para=(pi,a10,a20,b0,a11,a21,b1)
  #observed is the observed value for the ith individual, a 2-dim vector
  #measure is the random effect
  pi <- para[1]
  a10 <- para[2]
  a20 <- para[3]
  b0 <- para[4]
  a11 <- para[5]
  a21 <- para[6]
  b1 <- para[7]
  t1 <- observed[1]  
  t2 <- observed[2]
  index1 <- pi * ((pnorm(a11+b1*measure))^t1) * ((1-(pnorm(a11+b1*measure)))^(1-t1)) * 
    ((pnorm(a21+b1*measure))^t2) * ((1-(pnorm(a21+b1*measure)))^(1-t2))
  index0 <- (1-pi) * ((pnorm(a10+b0*measure))^(1-t1)) * ((1-(pnorm(a10+b0*measure)))^(t1)) *
    ((pnorm(a20+b0*measure))^(1-t2)) * ((1-(pnorm(a20+b0*measure)))^(t2))
  p <- index1/(index1+index0)  
  return(p)
}




#Calculate the joint posterior distribution for (a11,a21,b1)
Lik_a1b1 <- function(a11, a21, b1, result, truestatus, rf, mua1, sigmaa1, mua2, sigmaa2, mub1, sigmab1){
  N <- dim(result)[1]
  
  li <- rep(0, N)
  for (i in 1:N){
    li[i] <- ((pnorm(a11+b1*rf[i]))^(truestatus[i]*result[i,1])) * 
      ((1-(pnorm(a11+b1*rf[i])))^(truestatus[i]*(1-result[i,1]))) * 
      ((pnorm(a21+b1*rf[i]))^(truestatus[i]*result[i,2])) * 
      ((1-(pnorm(a21+b1*rf[i])))^(truestatus[i]*(1-result[i,2])))
  }
  lik <- prod(li)*dnorm(a11, mua1, sigmaa1)*dnorm(a21, mua2, sigmaa2)*dnorm(b1, mub1, sigmab1)
  return(lik)
}



#Generate samples of (a11,a21,b1) using sampling-importance resampling
SIR_a1b1 <- function(result, truestatus, rf, mua1, sigmaa1, mua2, sigmaa2, mub1, sigmab1, m){
  samplea11 <- rnorm(m, -0.5, 1)  
  samplea21 <- rnorm(m, 1.5, 1)
  sampleb1 <- rnorm(m, 1.5, 1)
  li <- rep(0, m)
  g <- rep(0, m)
  weight <- rep(0,m)
  for (i in 1:m){
    li[i] <- Lik_a1b1(samplea11[i], samplea21[i], sampleb1[i], result, truestatus, rf, mua1, sigmaa1, mua2, sigmaa2, mub1, sigmab1)
    g[i] <- dnorm(samplea11[i], -0.5, 1)*dnorm(samplea21[i], 1.5, 1)*dnorm(sampleb1[i], 1.5, 1)
    weight[i] <- li[i]/g[i]  
  }
  
  weight <- weight/sum(weight)
  position <- sample(1:m, 1, prob=weight)
  a11 <- samplea11[position]
  a21 <- samplea21[position]
  b1 <- sampleb1[position]
  return(c(a11,a21,b1))
}









#Calculate the joint posterior distribution for (a10,a20,b0)
Lik_a0b0 <- function(a10, a20, b0, result, truestatus, rf, mua1, sigmaa1, mua2, sigmaa2, mub0, sigmab0){
  N <- dim(result)[1]
  li <- rep(0, N)
  for (i in 1:N){
    li[i] <- ((pnorm(a10+b0*rf[i]))^(truestatus[i]*(1-result[i,1]))) * 
      ((1-(pnorm(a10+b0*rf[i])))^(truestatus[i]*(result[i,1]))) * 
      ((pnorm(a20+b0*rf[i]))^(truestatus[i]*(1-result[i,2]))) * 
      ((1-(pnorm(a20+b0*rf[i])))^(truestatus[i]*(result[i,2])))
  }
  lik <- prod(li)*dnorm(a10, mua1, sigmaa1)*dnorm(a20, mua2, sigmaa2)*dnorm(b0, mub0, sigmab0)
  return(lik)
}



#Generate samples of (a10,a20,b0) using sampling-importance resampling
SIR_a0b0 <- function(result, truestatus, rf, mua1, sigmaa1, mua2, sigmaa2, mub0, sigmab0, m){
  samplea10 <- rnorm(m, 3, 1)  
  samplea20 <- rnorm(m, 0, 1)
  sampleb0 <- rnorm(m, 1, 1)
  li <- rep(0, m)
  g <- rep(0, m)
  weight <- rep(0,m)
  for (i in 1:m){
    li[i] <- Lik_a0b0(samplea10[i], samplea20[i], sampleb0[i], result, truestatus, rf, mua1, sigmaa1, mua2, sigmaa2, mub0, sigmab0)
    g[i] <- dnorm(samplea10[i], 3, 1)*dnorm(samplea20[i], 0, 1)*dnorm(sampleb0[i], 1, 1)
    weight[i] <- li[i]/g[i]  
  }
  
  weight <- weight/sum(weight) 
  position <- sample(1:m, 1, prob=weight)
  a10 <- samplea10[position]
  a20 <- samplea20[position]
  b0 <- sampleb0[position]
  return(c(a10,a20,b0))
}





#Calculate the posterior distribution for random effect
Lik_r <- function(r, t1, t2, a10, a20, b0, a11, a21, b1, d){
  li <- ((pnorm(a11+b1*r))^(d*t1)) * 
    ((1-(pnorm(a11+b1*r)))^(d*(1-t1))) *
    ((pnorm(a21+b1*r))^(d*t2)) * 
    ((1-(pnorm(a21+b1*r)))^(d*(1-t2))) *
    ((pnorm(a10+b0*r))^((1-d)*(1-t1))) * 
    ((1-(pnorm(a10+b0*r)))^((1-d)*(t1))) *
    ((pnorm(a20+b0*r))^((1-d)*(1-t2))) * 
    ((1-(pnorm(a20+b0*r)))^((1-d)*(t2))) * dnorm(r)
  return(li)
}


#Generate samples of random effect using sampling-importance resampling
SIR_r <- function(t1, t2, a10, a20, b0, a11, a21, b1, d, m){
  sample1 <- rnorm(m, 0, 1)  
  li <- rep(0, m)
  g <- rep(0, m)
  weight <- rep(0,m)
  
  for (i in 1:m){
    li[i] <- Lik_r(sample1[i], t1, t2, a10, a20, b0, a11, a21, b1, d)
    g[i] <- dnorm(sample1[i], 0, 1)
    weight[i] <- li[i]/g[i]   
  }
  
  
  weight <- weight/sum(weight) 
  r <- sample(sample1, 1, prob=weight)
  return(r)
}








randomeffect_2 <- function(prior,observe,n,burn){
  #prior is a 7*2 matrix
  #the first row is the hyperparameter in beta prior distribution for prevalence
  #the second row is the hyperparameter in normal prior distribution for a10
  #the second row is the hyperparameter in normal prior distribution for a20
  #the second row is the hyperparameter in normal prior distribution for b0
  #the second row is the hyperparameter in normal prior distribution for a11
  #the second row is the hyperparameter in normal prior distribution for a21
  #the second row is the hyperparameter in normal prior distribution for b1
  #observe is a 2*2 matrix
  #n is the number of cycle
  #burn is the number of cycle for convergence
  
  X <- matrix(0, n, 7)  
  n_11 <- observe[1,1]  
  n_10 <- observe[1,2]  
  n_01 <- observe[2,1] 
  n_00 <- observe[2,2]  
  N <- n_11+n_10+n_01+n_00
  
  
  result1 <- matrix(0, nrow = N, ncol = 2)
  result1[1:(n_11+n_10), 1] <- 1 
  result1[1:n_11, 2] <- 1
  result1[(n_11+n_10+1):(n_11+n_10+n_01), 2] <- 1
  
 
  X[1, 1] <- rbeta(1,prior[1,1],prior[1,2])
  for (p in 2:7) {
    X[1, p] <- rnorm(1, prior[p,1], prior[p,2])
  }
  

  rf <- rnorm(N)
  
  
  for (i in 2:n){
    pi <- X[i-1, 1]  
    a10 <- X[i-1, 2]
    a20 <- X[i-1, 3]
    b0 <- X[i-1, 4]
    a11 <- X[i-1, 5]
    a21 <- X[i-1, 6]
    b1 <- X[i-1, 7]
  
    truestatus <- rep(0, N)
    for (j in 1:N){
      truestatus[j] <- rbinom(1, 1, truedisease(X[i-1,], result1[j,], rf[j]))
    }
    
    X[i, 1] <- rbeta(1, prior[1,1]+sum(truestatus), prior[1,2]+N-sum(truestatus))
    
    
    
    para0 <- SIR_a0b0(result1, truestatus, rf, prior[2,1], prior[2,2], prior[3,1], prior[3,2],
                      prior[4,1], prior[4,2], 100)
    X[i, 2:4] <- para0
    para1 <- SIR_a1b1(result1, truestatus, rf, prior[5,1], prior[5,2], prior[6,1], prior[6,2],
                      prior[7,1], prior[7,2], 100)
    X[i, 5:7] <- para1
    
    
    for (j in 1:N){
      rf[j] <- SIR_r(result1[j,1], result1[j,2], X[i,2], X[i,3], X[i,4], X[i,5], X[i,6], X[i,7], truestatus[j], 100)
    }
    cat(i,"completed\n")
  }
  
  b <- burn + 1
  x <- X[b:n, ] 

  return(x)
}

prior <- matrix(c(1,1,
                  2.171, 0.261,
                  0.692, 0.560,
                  0.861, 0.5,
                  -0.811,0.380,
                  1.012, 0.268,
                  0.668, 0.5),ncol=2,byrow = TRUE)
observe <- matrix(c(38,2,87,35),nrow=2,byrow = TRUE)                 
n <- 12000
burn <- 2000

result1 <- randomeffect_2(prior,observe,n,burn) 




mean(result1[,1])
mean(result1[,2])
mean(result1[,3])
mean(result1[,4])
mean(result1[,5])
mean(result1[,6])
mean(result1[,7])



se_1 <- pnorm(result1[,5]/(sqrt(1+result1[,7]^2)))
mean(se_1)
se_2 <- pnorm(result1[,6]/(sqrt(1+result1[,7]^2)))
mean(se_2)
sp_1 <- pnorm(result1[,2]/(sqrt(1+result1[,4]^2)))
mean(sp_1)
sp_2 <- pnorm(result1[,3]/(sqrt(1+result1[,4]^2)))
mean(sp_2)

