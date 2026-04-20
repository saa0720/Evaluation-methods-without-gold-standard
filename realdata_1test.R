#This code is for one test scenario without covariates.
library(ggplot2)
set.seed(0720)




r_1_1 <- function(api,bpi,as,bs,ac,bc,result,n,burn){
  #api,bpi are the hyperparameters in the prior beta distribution for prevalence
  #as,bs are the hyperparameters in the prior beta distribution for sensitivity
  #ac,bc are the hyperparameters in the prior beta distribution for specificity
  #result is the observed data, which is a two-dimensional vector, the first component is the number of individuals tested positive
  #n is the number of cycles
  #burn is the number of cycles which are used to achieve convergence
  X <- matrix(0, n, 3)  
  #generate initial values
  pi0 <- rbeta(1,api,bpi)
  se0 <- rbeta(1,as,bs)
  sp0 <- rbeta(1,ac,bc)
  X[1, ] <- c(pi0, se0, sp0)  
  n_1 <- result[1]  
  n_0 <- result[2]  
  N <- n_1+n_0  
  for (i in 2:n){
    pi <- X[i-1, 1]  
    se <- X[i-1, 2]  
    sp <- X[i-1, 3]  
    s_1 <- rbinom(1, n_1, (pi*se)/(pi*se+(1-pi)*(1-sp)))  #s_1=N(T=1,D=1), latent data
    s_0 <- rbinom(1, n_0, (pi*(1-se))/(pi*(1-se)+(1-pi)*sp))  #s_0=N(T=0,D=1), latent data
    X[i, 1] <- rbeta(1, (s_1+s_0+api), (N-(s_1+s_0)+bpi))  
    X[i, 2] <- rbeta(1, (s_1+as), (s_0+bs))  
    X[i, 3] <- rbeta(1, (n_0-s_0+ac), (n_1-s_1+bc))  
  }
  b <- burn + 1
  x <- X[b:n, ]  
  return(x)
}




result <- r_1_1(1,1,
                4.44,13.31,
                71.25,3.75,
                c(40,122),
                12000,2000)
mean(result[,1])
mean(result[,2])
mean(result[,3])














#Plotting the prior and posterior distributions
set.seed(0720)
pi_prior <- rbeta(10000,1,1)
pi_plot <- cbind(pi_prior,result[,1])
colnames(pi_plot) <- c("prior","posterior")
pi_plot <- as.data.frame(pi_plot)
p1 <- ggplot() +
  geom_density(data = pi_plot, aes(x = prior, y = ..density..), linewidth = 0.8, alpha = 0.15, fill = "#e18d82", color = "#e18d82", show.legend = FALSE) +
  geom_density(data = pi_plot, aes(x = posterior, y = ..density..), linewidth = 0.8, fill = "#4b72db", alpha = 0.15, color = "#4b72db", show.legend = FALSE) +
  geom_vline(aes(xintercept = mean(pi_plot[,2]), color = "posterior"), linetype = "dashed", linewidth = 0.8) + 
  geom_vline(aes(xintercept = mean(pi_plot[,1]), color = "prior"), linetype = "dashed", linewidth = 0.8) +
  theme_bw() +
  labs(x = "prevalence", y = "density") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +  
  scale_color_manual(name = "distribution",
                     values = c('prior' = '#e18d82', 'posterior' = '#4b72db')) + 
  guides(color = guide_legend(override.aes = list(linetype = "solid", size = 1.5)),
         linetype = guide_legend(override.aes = list(size = 1.5))) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.14, 0.90),
    legend.key.size = unit(0.5, "lines"),  
    legend.text = element_text(size = 4)  
  )






set.seed(0720)
se_prior <- rbeta(10000,4.44,13.31)
se_plot <- cbind(se_prior,result[,2])
colnames(se_plot) <- c("prior","posterior")
se_plot <- as.data.frame(se_plot)
p2 <- ggplot() +
  geom_density(data = se_plot, aes(x = prior, y = ..density.., color = "#e18d82"), linewidth = 0.8, alpha = 0.15, fill = "#e18d82") +
  geom_density(data = se_plot, aes(x = posterior, y = ..density.., color = "#4b72db"), linewidth = 0.8, fill = "#4b72db", alpha = 0.15) +
  geom_vline(xintercept = mean(se_plot[,2]), colour = "#4b72db", linetype = "dashed", linewidth = 0.8) + 
  geom_vline(xintercept = mean(se_plot[,1]), colour = "#e18d82", linetype = "dashed", linewidth = 0.8) +
  theme_bw() +
  labs(x = "sensitivity", y = "density") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +  
  scale_color_manual(name = "distribution",
                     values = c('#e18d82' = '#e18d82', "#4b72db" = '#4b72db'), 
                     aesthetics = c("colour"),
                     labels = c('prior', 'posterior')) + 
  theme(legend.position = "none")








set.seed(0720)
sp_prior <- rbeta(10000,71.25,3.75)
sp_plot <- cbind(sp_prior,result[,3])
colnames(sp_plot) <- c("prior","posterior")
sp_plot <- as.data.frame(sp_plot)
p3 <- ggplot() +
  geom_density(data = sp_plot, aes(x = prior, y = ..density.., color = "#e18d82"), linewidth = 0.8, alpha = 0.15, fill = "#e18d82") +
  geom_density(data = sp_plot, aes(x = posterior, y = ..density.., color = "#4b72db"), linewidth = 0.8, fill = "#4b72db", alpha = 0.15) +
  geom_vline(xintercept = mean(sp_plot[,2]), colour = "#4b72db", linetype = "dashed", linewidth = 0.8) + 
  geom_vline(xintercept = mean(sp_plot[,1]), colour = "#e18d82", linetype = "dashed", linewidth = 0.8) +
  theme_bw() +
  labs(x = "specificity", y = "density") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) + 
  scale_color_manual(name = "distribution",
                     values = c('#e18d82' = '#e18d82', "#4b72db" = '#4b72db'), 
                     aesthetics = c("colour"),
                     labels = c('prior', 'posterior')) + 
  theme(legend.position = "none")





p_1 <- ggpubr::ggarrange(p1,p2,p3,nrow=1,ncol=3)
p_1
