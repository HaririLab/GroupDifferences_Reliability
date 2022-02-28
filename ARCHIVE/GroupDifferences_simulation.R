# library(shiny)
# library(cowplot)
library(ggplot2)
library(ggpubr)
library("psych") # inc cohen.d
library("MASS")
library(pwr)
library(lme4)

N <- 50
N_iter <- 500
input_gp1_mean <- 0
input_gp2_mean <- 1
sds <- seq(.25,5,.25)
iccs <- seq(0,1,.2)

#### FIRST PASS ####

simulations <- data.frame(ICC=double(), cor1=double(), cor2=double(), sd=double(), iter=integer(), p_t1=double(), p_t2=double())

for (sd in sds) {
  for (icc in iccs) {
    for (i in 1:N_iter) {
    
      input_gp1_sd <- sd
      input_gp2_sd <- sd
      input_ICC <- icc
      
      # run calculations
      mu_gp1 <- rep(input_gp1_mean, 2)
      var_gp1 <- input_gp1_sd*input_gp1_sd # use the same sd/variance for t1 and t2
      covar_gp1 <- var_gp1*input_ICC # use formula r(x1,x2) = cov(x1,x2)/sd(x1)*sd(x2) to calculate what covariance we need to achieve desired correlation / ICC
      sigma_gp1 <- matrix(c(var_gp1, covar_gp1,
                            covar_gp1, var_gp1), ncol=2, byrow=TRUE)
      sim_gp1 <- mvrnorm(n=N, mu=mu_gp1, Sigma=sigma_gp1)
      
      mu_gp2 <- rep(input_gp2_mean, 2)
      var_gp2 <- input_gp2_sd*input_gp2_sd # use the same sd/variance for t1 and t2
      covar_gp2 <- var_gp2*input_ICC # use formula r(x1,x2) = cov(x1,x2)/sd(x1)*sd(x2) to calculate what covariance we need to achieve desired correlation / ICC
      sigma_gp2 <- matrix(c(var_gp2, covar_gp2,
                            covar_gp2, var_gp2), ncol=2, byrow=TRUE)
      sim_gp2 <- mvrnorm(n=N, mu=mu_gp2, Sigma=sigma_gp2)
      
      df1 <- data.frame(value_t1=sim_gp1[,1], value_t2=sim_gp1[,2], group="Group 1")
      df2 <- data.frame(value_t1=sim_gp2[,1], value_t2=sim_gp2[,2], group="Group 2")
      df <- rbind(df1,df2)
      
      # ggplot(df, aes(x=group, y=value_t1)) +
      #   geom_violin()
      # 
      # ggplot(df, aes(x=group, y=value_t2)) +
      #   geom_violin()
      
      c1 <- cor(df[which(df$group=="Group 1"),"value_t1"], df[which(df$group=="Group 1"),"value_t2"])
      c2 <- cor(df[which(df$group=="Group 2"),"value_t1"], df[which(df$group=="Group 2"),"value_t2"])
      
      test_t1 <- t.test(df[which(df$group=="Group 1"),"value_t1"], df[which(df$group=="Group 2"),"value_t1"])
      test_t2 <- t.test(df[which(df$group=="Group 1"),"value_t2"], df[which(df$group=="Group 2"),"value_t2"])
      
      simulations <- rbind( simulations, data.frame(ICC=icc, cor1=c1, cor2=c2, sd=sd, iter=i, p_t1=test_t1$p.value, p_t2=test_t2$p.value) )
      
    }
  }
  print(paste("Finished sd",sd))
}

simulations$sd <- as.factor(simulations$sd)
ggplot(simulations, aes(x=sd, y=p_t1)) + facet_grid(ICC ~ .) + geom_violin()

simulations_mean_t1ps <- array(NA, c(length(iccs),length(sds)) )
for (i in 1:length(iccs)) {
  for (j in 1:length(sds)) {
    simulations_mean_t1ps[i,j] <- mean(simulations[which(simulations$ICC==iccs[i] & simulations$sd==sds[j]), "p_t1"])
  }
}
rownames(simulations_mean_t1ps) <- iccs
colnames(simulations_mean_t1ps) <- sds


#### SECOND PASS ####

simulations2 <- data.frame(cor1=double(), cor2=double(), true_sd=double(), noise_sd=double(), iter=integer(), p_t1=double(), p_t2=double(), d_true=double(), d_t1=double(), d_t2=double())
sds_noise <- seq(.5,4,.5)

for (true_sd in sds) {
  for(noise_sd in sds_noise) {
    for (i in 1:N_iter) {
      
      input_gp1_sd <- true_sd
      input_gp2_sd <- true_sd

      true_gp1 <- rnorm(n=N, mean=input_gp1_mean, sd=input_gp1_sd)
      true_gp2 <- rnorm(n=N, mean=input_gp2_mean, sd=input_gp2_sd)
      
      obs_gp1_t1 <- true_gp1 + rnorm(n=N, mean=0, sd=noise_sd)
      obs_gp2_t1 <- true_gp2 + rnorm(n=N, mean=0, sd=noise_sd)
      obs_gp1_t2 <- true_gp1 + rnorm(n=N, mean=0, sd=noise_sd)
      obs_gp2_t2 <- true_gp2 + rnorm(n=N, mean=0, sd=noise_sd)      
      
      c1 <- cor(obs_gp1_t1, obs_gp1_t2)
      c2 <- cor(obs_gp2_t1, obs_gp2_t2)
      # ICC(cbind(obs_gp2_t1, obs_gp2_t2)) # or could use this
      
      d_true <- cohen.d(c(true_gp1,true_gp2), c(rep('gp1',N), rep('gp2',N)))$estimate
      d_t1 <- cohen.d(c(obs_gp1_t1,obs_gp2_t1), c(rep('gp1',N), rep('gp2',N)))$estimate
      d_t2 <- cohen.d(c(obs_gp1_t2,obs_gp2_t2), c(rep('gp1',N), rep('gp2',N)))$estimate
      
      test_t1 <- t.test(obs_gp1_t1, obs_gp2_t1)
      test_t2 <- t.test(obs_gp1_t2, obs_gp2_t2)
      
      simulations2 <- rbind( simulations2, data.frame(cor1=c1, cor2=c2, true_sd=true_sd, noise_sd=noise_sd, iter=i, p_t1=test_t1$p.value, p_t2=test_t2$p.value, d_true=d_true, d_t1=d_t1, d_t2=d_t2) )
      
    }
  }
  print(paste("Finished sd",true_sd))
}

# simulations$sd <- as.factor(simulations$sd)
# ggplot(simulations, aes(x=sd, y=p_t1)) + facet_grid(ICC ~ .) + geom_violin()
# 
simulations2_mean_t1ps <- array(NA, c(length(sds),length(sds_noise)) )
simulations2_mean_c1s <- array(NA, c(length(sds),length(sds_noise)) )
simulations2_mean_d.true <- array(NA, c(length(sds),length(sds_noise)) )
simulations2_mean_d.t1 <- array(NA, c(length(sds),length(sds_noise)) )
for (i in 1:length(sds)) {
  for (j in 1:length(sds_noise)) {
    simulations2_mean_t1ps[i,j] <- mean(simulations2[which(simulations2$true_sd==sds[i] & simulations2$noise_sd==sds_noise[j]), "p_t1"])
    simulations2_mean_c1s[i,j] <- mean(simulations2[which(simulations2$true_sd==sds[i] & simulations2$noise_sd==sds_noise[j]), "cor1"])
    simulations2_mean_d.true[i,j] <- mean(simulations2[which(simulations2$true_sd==sds[i] & simulations2$noise_sd==sds_noise[j]), "d_true"])
    simulations2_mean_d.t1[i,j] <- mean(simulations2[which(simulations2$true_sd==sds[i] & simulations2$noise_sd==sds_noise[j]), "d_t1"])
  }
}
rownames(simulations2_mean_t1ps) <- sds;    colnames(simulations2_mean_t1ps) <- sds_noise
rownames(simulations2_mean_c1s) <- sds;     colnames(simulations2_mean_c1s) <- sds_noise
rownames(simulations2_mean_d.true) <- sds;  colnames(simulations2_mean_d.true) <- sds_noise
rownames(simulations2_mean_d.t1) <- sds;    colnames(simulations2_mean_d.t1) <- sds_noise

