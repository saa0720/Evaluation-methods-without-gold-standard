#fixed effects model
set.seed(0720)

#Lik_Se1 is the likelihood function of Se1
Lik_Se1 <- function(S1, S2,covs,y_11,y_10,y_01,y_00,as1,bs1,us,bcovs){
  li1 <- (S1*S2+covs)^y_11
  li2 <- (S1*(1-S2)-covs)^y_10*10^200
  li3 <- ((1-S1)*S2-covs)^y_01
  li4 <- ((1-S1)*(1-S2)+covs)^y_00
  li5 <- (S1^(as1-1))*((1-S1)^(bs1-1))*((us-covs)^(bcovs-1))
  li <- li1*li2*li3*li4*li5
  return(li)
}



#SIR_Se1 is the function to generate samples of Se1 using sampling-importance resampling
SIR_Se1 <- function(S2,covs,y_11,y_10,y_01,y_00,as1,bs1,us,bcovs,m){
  #m is the number of samples
  sample1 <- runif(m, covs/(1-S2), 1-covs/S2)
  li <- rep(0, m)
  g <- rep(0, m)
  weight <- rep(0,m)
  for (j in 1:m){
    li[j] <- Lik_Se1(sample1[j],S2,covs,y_11,y_10,y_01,y_00,as1,bs1,us,bcovs)
    g[j] <- dbeta(sample1[j], 1, 1)
    if ((li[j]<=0)||(is.na(li[j]))){
      weight[j] <- 0
    }
    else {
      weight[j] <- li[j]/g[j]
    }
  }
  
  weight <- weight/sum(weight) 
  S1 <- sample(sample1, 1, prob=weight)
  return(S1)
}





#Lik_Se2 is the likelihood function of Se2
Lik_Se2 <- function(S1, S2,covs,y_11,y_10,y_01,y_00,as2,bs2,us,bcovs){
  li1 <- (S1*S2+covs)^y_11
  li2 <- (S1*(1-S2)-covs)^y_10
  li3 <- ((1-S1)*S2-covs)^y_01
  li4 <- ((1-S1)*(1-S2)+covs)^y_00*10^200
  li5 <- (S2^(as2-1))*((1-S2)^(bs2-1))*((us-covs)^(bcovs-1))
  li <- li1*li2*li3*li4*li5
  return(li)
}



#SIR_Se2 is the function to generate samples of Se2 using sampling-importance resampling
SIR_Se2 <- function(S1,covs,y_11,y_10,y_01,y_00,as2,bs2,us,bcovs,m){
  sample1 <- runif(m, covs/(1-S1), 1-covs/S1)
  li <- rep(0, m)
  g <- rep(0, m)
  weight <- rep(0,m)
  for (i in 1:m){
    li[i] <- Lik_Se2(S1,sample1[i],covs,y_11,y_10,y_01,y_00,as2,bs2,us,bcovs)
    g[i] <- dbeta(sample1[i], 1, 1)
    if (li[i]<=0){
      weight[i] <- 0
    }
    else {
      weight[i] <- li[i]/g[i]
    }
  }
  
  weight <- weight/sum(weight) 
  S2 <- sample(sample1, 1, prob=weight)
  return(S2)
}





#Lik_Sp1 is the likelihood function of Sp1
Lik_Sp1 <- function(C1, C2, covc,n_11,n_10,n_01,n_00,y_11,y_10,y_01,y_00,ac1,bc1,uc,bcovc){
  li1 <- ((1-C1)*(1-C2)+covc)^(n_11-y_11)
  li2 <- ((1-C1)*C2-covc)^(n_10-y_10)
  li3 <- (C1*(1-C2)-covc)^(n_01-y_01)
  li4 <- (C1*C2+covc)^(n_00-y_00)
  li5 <- (C1^(ac1-1))*((1-C1)^(bc1-1))*((uc-covc)^(bcovc-1))
  li <- li1*li2*li3*li4*li5
  return(li)
}



#SIR_Sp1 is the function to generate samples of Sp1 using sampling-importance resampling
SIR_Sp1 <- function(C2,covc,n_11,n_10,n_01,n_00,y_11,y_10,y_01,y_00,ac1,bc1,uc,bcovc,m){
  #m代表生成的sample个数
  #假设sampling importance resampling的g函数为beta(ac1,bc1)，保证了定义域相同
  #从g函数中sample m个样本
  sample1 <- rbeta(m, 1, 1)
  #计算每个sample的似然函数
  li <- rep(0, m)
  #计算每个sample在生成函数g中的似然
  g <- rep(0, m)
  #计算每个sample对应的weight
  weight <- rep(0,m)
  for (i in 1:m){
    li[i] <- Lik_Sp1(sample1[i],C2,covc,n_11,n_10,n_01,n_00,y_11,y_10,y_01,y_00,ac1,bc1,uc,bcovc)
    g[i] <- dbeta(sample1[i], 1, 1)
    if (li[i]<=0){
      weight[i] <- 0
    }
    else {
      weight[i] <- li[i]/g[i]
    }
  }
  
  weight <- weight/sum(weight) #将weight归一化，作为sample的prob
  #按照weight从原来的y中sample出一个样本，作为Se1的值
  C1 <- sample(sample1, 1, prob=weight)
  return(C1)
  
}





#Lik_Sp2 is the likelihood function of Sp2
Lik_Sp2 <- function(C1, C2,covc,n_11,n_10,n_01,n_00,y_11,y_10,y_01,y_00,ac2,bc2,uc,bcovc){
  li1 <- ((1-C1)*(1-C2)+covc)^(n_11-y_11)
  li2 <- ((1-C1)*C2-covc)^(n_10-y_10)
  li3 <- (C1*(1-C2)-covc)^(n_01-y_01)
  li4 <- (C1*C2+covc)^(n_00-y_00)
  li5 <- (C2^(ac2-1))*((1-C2)^(bc2-1))*((uc-covc)^(bcovc-1))
  li <- li1*li2*li3*li4*li5
  return(li)
}



#SIR_Sp2 is the function to generate samples of Sp2 using sampling-importance resampling
SIR_Sp2 <- function(C1,covc,n_11,n_10,n_01,n_00,y_11,y_10,y_01,y_00,ac2,bc2,uc,bcovc,m){
  sample1 <- rbeta(m, 1, 1)
  li <- rep(0, m)
  g <- rep(0, m)
  weight <- rep(0,m)
  for (i in 1:m){
    li[i] <- Lik_Sp2(C1,sample1[i],covc,n_11,n_10,n_01,n_00,y_11,y_10,y_01,y_00,ac2,bc2,uc,bcovc)
    g[i] <- dbeta(sample1[i], 1, 1)
    if (li[i]<=0){
      weight[i] <- 0
    }
    else {
      weight[i] <- li[i]/g[i]
    }
  }
  
  weight <- weight/sum(weight) 
  C2 <- sample(sample1, 1, prob=weight)
  return(C2)
}





#Lik_covs is the likelihood function of covs
Lik_covs <- function(S1, S2, covs, y_11,y_10,y_01,y_00,us,acovs,bcovs){
  li1 <- (S1*S2+covs)^(y_11)
  li2 <- (S1*(1-S2)-covs)^(y_10)
  li3 <- ((1-S1)*S2-covs)^(y_01)
  li4 <- ((1-S1)*(1-S2)+covs)^(y_00)
  li5 <- (covs^(acovs-1))*((us-covs)^(bcovs-1))
  li <- li1*li2*li3*li4*li5
  return(li)
}



#SIR_covs is the function to generate samples of covs using sampling-importance resampling
SIR_covs <- function(S1, S2, y_11,y_10,y_01,y_00,us,acovs,bcovs,m){
  sample1 <- runif(m, 0, us)
  li <- rep(0, m)
  g <- rep(0, m)
  weight <- rep(0,m)
  for (i in 1:m){
    li[i] <- Lik_covs(S1,S2,sample1[i],y_11,y_10,y_01,y_00,us,acovs,bcovs)
    g[i] <- dunif(sample1[i], 0, us)
    if (li[i]<=0){
      weight[i] <- 0
    }
    else {
      weight[i] <- li[i]/g[i]
    }
  }
  
  weight <- weight/sum(weight) 
  covs <- sample(sample1, 1, prob=weight)
  return(covs)
}





#Lik_covc is the likelihood function of covc
Lik_covc <- function(C1, C2, covc, n_11,n_10,n_01,n_00,y_11,y_10,y_01,y_00,uc,acovc,bcovc){
  li1 <- ((1-C1)*(1-C2)+covc)^(n_11-y_11)*10^200
  li2 <- ((1-C1)*C2-covc)^(n_10-y_10)
  li3 <- (C1*(1-C2)-covc)^(n_01-y_01)
  li4 <- (C1*C2+covc)^(n_00-y_00)
  li5 <- (covc^(acovc-1))*((uc-covc)^(bcovc-1))
  li <- li1*li2*li3*li4*li5
  return(li)
}



#SIR_covc is the function to generate samples of covc using sampling-importance resampling
SIR_covc <- function(C1, C2, n_11,n_10,n_01,n_00,y_11,y_10,y_01,y_00,uc,acovc,bcovc,m){
  sample1 <- runif(m, 0, uc)
  li <- rep(0, m)
  g <- rep(0, m)
  weight <- rep(0,m)
  for (i in 1:m){
    li[i] <- Lik_covc(C1,C2,sample1[i],n_11,n_10,n_01,n_00,y_11,y_10,y_01,y_00,uc,acovc,bcovc)
    g[i] <- dunif(sample1[i], 0, uc)
    if (li[i]<=0){
      weight[i] <- 0
    }
    else {
      weight[i] <- li[i]/g[i]
    }
  }
  
  weight <- weight/sum(weight) 
  covc <- sample(sample1, 1, prob=weight)
  return(covc)
}





fixedeffect_2_1 <- function(api,bpi,as1,bs1,as2,bs2,ac1,bc1,ac2,bc2,acovs,bcovs,acovc,bcovc,result,n,burn){
  #api,bpi are the hyperparameters in the prior beta distribution for prevalence
  #as1,bs1 are the hyperparameters in the prior beta distribution for sensitivity of test1
  #ac1,bc1 are the hyperparameters in the prior beta distribution for specificity of test1
  #as2,bs2 are the hyperparameters in the prior beta distribution for sensitivity of test2
  #ac2,bc2 are the hyperparameters in the prior beta distribution for specificity of test2
  #acovs,bcovs are the hyperparameters in the prior beta distribution for covariance in the diseased population
  #acovc,bcovc are the hyperparameters in the prior beta distribution for covariance in the nondiseased population
  #result is the observed data, which is a 2*2 matrix, the first row is n_11,n_10, the second row is n_01,n_00
  #n_10=N(T_1=1,T_2=0), n_01=N(T_1=0,T_2=1)
  #n is the number of cycles
  #burn is the number of cycles which are used to achieve convergence
  X <- matrix(0, n, 7)  
  X[1, 1] <- rbeta(1, api, bpi)
  X[1, 2] <- rbeta(1, as1, bs2)
  X[1, 3] <- rbeta(1, ac1, bc1)
  X[1, 4] <- rbeta(1, as2, bs2)
  X[1, 5] <- rbeta(1, ac2, bc2)
  us <- min(X[1,2], X[1,4])-X[1,2]*X[1,4]
  X[1, 6] <- runif(1, 0, us)
  uc <- min(X[1,3], X[1,5])-X[1,3]*X[1,5]
  X[1, 7] <- runif(1, 0, uc)
  n_11 <- result[1,1]  
  n_10 <- result[1,2]  
  n_01 <- result[2,1]  
  n_00 <- result[2,2]  
  N <- n_11+n_10+n_01+n_00  
  for (i in 2:n){
    pi1 <- X[i-1, 1]  
    S1 <- X[i-1, 2]  
    C1 <- X[i-1, 3]  
    S2 <- X[i-1, 4]  
    C2 <- X[i-1, 5]  
    covs <- X[i-1, 6]  
    covc <- X[i-1, 7]  
    
    
    p_11_num <- pi1*(S1*S2+covs)
    p_11_den2 <- (1-pi1)*((1-C1)*(1-C2)+covc)
    p_11 <- p_11_num/(p_11_num+p_11_den2)
    y_11 <- rbinom(1, n_11, p_11)  
    
    p_10_num <- pi1*(S1*(1-S2)-covs)
    p_10_den2 <- (1-pi1)*((1-C1)*C2-covc)
    p_10 <- p_10_num/(p_10_num+p_10_den2)
    y_10 <- rbinom(1, n_10, p_10) 
    
    p_01_num <- pi1*((1-S1)*S2-covs)
    p_01_den2 <- (1-pi1)*(C1*(1-C2)-covc)
    p_01 <- p_01_num/(p_01_num+p_01_den2)
    y_01 <- rbinom(1, n_01, p_01) 
    
    p_00_num <- pi1*((1-S1)*(1-S2)+covs)
    p_00_den2 <- (1-pi1)*(C1*C2+covc)
    p_00 <- p_00_num/(p_00_num+p_00_den2) 
    y_00 <- rbinom(1, n_00, p_00)  
    
    
    X[i, 1] <- rbeta(1, (api+y_11+y_10+y_01+y_00), (bpi+N-(y_11+y_10+y_01+y_00))) 
    
   
    X[i, 2] <- SIR_Se1(X[(i-1),4],X[(i-1),6],y_11,y_10,y_01,y_00,as1,bs1,us,bcovs,500) 
    
    
  
    X[i, 4] <- SIR_Se2(X[i-1,2], X[i-1,6], y_11, y_10, y_01, y_00, as2, bs2, us, bcovs, 500)
    
   
    X[i, 3] <- SIR_Sp1(X[i-1,5], X[i-1,7], n_11,n_10,n_01,n_00,y_11,y_10,y_01,y_00,ac1,bc1,uc,bcovc,500)
    
   
    X[i, 5] <- SIR_Sp2(X[i-1,3], X[i-1,7], n_11,n_10,n_01,n_00,y_11,y_10,y_01,y_00,ac2,bc2,uc,bcovc,500)
    
   
    us <- (min(X[i,2], X[i,4]))-(X[i,2]*X[i,4])
    X[i, 6] <- SIR_covs(X[i,2], X[i,4], y_11,y_10,y_01,y_00,us,acovs,bcovs,500)
    
    uc <- (min(X[i,3],X[i,5]))-(X[i,3]*X[i,5])
    X[i, 7] <- SIR_covc(X[i,3], X[i,5], n_11,n_10,n_01,n_00,y_11,y_10,y_01,y_00,uc,acovc,bcovc,500)
   
  }
  b <- burn + 1
  x <- X[b:n, ]  
  return(x)
}


observed <- matrix(c(38,2,87,35), nrow=2,byrow = TRUE)


result <- fixedeffect_2_1(1,1,
                          4.44,13.31,
                          21.96,5.49,
                          71.25,3.75,
                          4.1,1.76,
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



