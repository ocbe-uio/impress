design1 <- matrix(c(0,0,0,
                    1,0,0,
                    1,1,0,
                    1,1,1,
                    0,0,0),nrow = 5, ncol = 3, byrow=TRUE)

design2 <- matrix(c(0,0,0,
                    1,1,1,
                    1,1,0,
                    1,0,0,
                    0,0,0),nrow = 5, ncol = 3, byrow=TRUE)

ngrp <- 3 # number of treatment groups
mu_trt <- 20 #treatment effect
mu_j <- c(0,0,0,0,0) #time effect
t <- length(mu_j) #number of timepoints
sd_j <- 10.6 #SD of intercept
sd_res <- 10.6 #SD of residual
n.iter <- 500 #number of simulated trials
ssize <- c(2,3,4) # Sample sizes per intervention group to form the results
N <- 50

library(tidyverse)
simpop <- tibble(pid = rep(1:(ngrp*N),t)) %>%
  arrange(pid) %>%
  mutate(group = rep(1:ngrp,each=N*t)) %>%
  group_by(pid) %>%
  mutate(time = 1:t-1) %>%
  mutate(trt = design1[,max(group)]) %>%
  mutate(mu_j = mu_j) %>%
  modify_at(c("pid", "time", "group", "trt"),as.factor) %>%
  mutate(a_i=rnorm(1, sd=sd_j)) %>%
  mutate(y=a_i + mu_j + as.numeric(trt)*mu_trt + rnorm(t, sd=sd_res)) %>%
  ungroup 
  
res1 <- nlme::lme(y ~ time + trt, data=simpop, random = ~ 1|pid)
res1

