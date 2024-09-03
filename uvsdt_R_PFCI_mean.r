library(tidyverse)

dat <- read.table("subdata.csv", header=T,sep=",") 

result = dat %>%             
  filter(TorF1 == 0 | TorF1 == 1) %>% # dual-task trials
  select(subnum, exposure, maskwidth, colorcond, TorF1, colorres, TorF2, gray_gray, gray_mix, gray_color, mix_gray, mix_mix, mix_color, color_gray, color_mix, color_color) %>%    
  group_by(subnum, colorcond, exposure, maskwidth) %>%                
  summarize(gray_gray= sum(gray_gray), gray_mix= sum(gray_mix), gray_color= sum(gray_color), mix_gray= sum(mix_gray), mix_mix= sum(mix_mix), mix_color= sum(mix_color), color_gray= sum(color_gray), color_mix= sum(color_mix), color_color= sum(color_color)) #%>%       

#基準となるパラメータ
sigma83ms <- 1
mu83ms_gray <- 0
sigma117ms <- 1
sigma150ms <- 1

### Function for model fitting
fit_uvsdt_mle <- function(data, add_constant = TRUE) {
  
  # add_constant = TRUE adds a small value to the response frequency vectors.
  # correction against extreme estimates 
  if (add_constant) {
    data <- data+0.5
  }
  
  # initial guess for parameter values d-primeベースで出す
  mu83ms_9deg  <- qnorm(data[4,2]/(data[4,1]+data[4,2]))-qnorm(data[1,2]/(data[1,1]+data[1,2]))
  mu83ms_13deg <- qnorm(data[5,2]/(data[5,1]+data[5,2]))-qnorm(data[1,2]/(data[1,1]+data[1,2]))
  mu83ms_17deg <- qnorm(data[6,2]/(data[6,1]+data[6,2]))-qnorm(data[1,2]/(data[1,1]+data[1,2]))
  mu83ms_21deg <- qnorm(data[7,2]/(data[7,1]+data[7,2]))-qnorm(data[1,2]/(data[1,1]+data[1,2]))
  mu83ms_25deg <- qnorm(data[8,2]/(data[8,1]+data[8,2]))-qnorm(data[1,2]/(data[1,1]+data[1,2]))
  mu83ms_color <- qnorm(data[19,2]/(data[19,1]+data[19,2]))-qnorm(data[1,2]/(data[1,1]+data[1,2]))
  mu117ms <- 1
  mu150ms <- 1

  cri <- 2
  
  guess <- c(mu83ms_9deg, mu83ms_13deg, mu83ms_17deg, mu83ms_21deg, mu83ms_25deg, mu83ms_color, mu117ms, mu150ms, cri)
  
  # model fit
  fit <- suppressWarnings(optim(uvsdt_logL, 
                                par = guess, 
                                lower =      c(0.1, 0.1, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5,  0.5),
                                upper =      c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 2.0, 2.0,  5.0), 
                                gr = NULL, method = "BFGS", 
                                control = list("maxit" = 10000, 
                                               "parscale" = c(1,   1,    1,   1,    1,   1,  0.001,  0.001,   1))))
  
  # outputs
  mu83ms_9deg  <- fit$par[1]
  mu83ms_13deg <- fit$par[2] 
  mu83ms_17deg <- fit$par[3] 
  mu83ms_21deg <- fit$par[4] 
  mu83ms_25deg <- fit$par[5] 
  mu83ms_color <- fit$par[6] 
  mu117ms<-  fit$par[7] 
  mu150ms <- fit$par[8] 
  cri <- fit$par[9] 
  
  logL <- -fit$value
  
  AIC <- -2 * logL + 2 * 9
  
  predicted_data <- matrix(c(rep(1,21)), nrow=21, ncol=2)
  ##############################################################################
  # model predictions
  predicted_data <- matrix(NA, nrow=21, ncol=2)
  #gray
  mean_one_mat <- c(mu83ms_gray, mu83ms_gray*mu117ms, mu83ms_gray*mu150ms)
  sigmamat <- c(sigma83ms,sigma117ms,sigma150ms)
  for (cond in 1:3) {
    pred_gray_others_rate <- pnorm(cri, mean= mean_one_mat[cond], sd=sigmamat[cond])
    pred_gray_color_rate <- 1-pred_gray_others_rate
    
    #反応数
    pred_nr_gray_color <- sum(data[cond,1]+data[cond,2]) * pred_gray_color_rate
    pred_nr_gray_others <- sum(data[cond,1]+data[cond,2]) * pred_gray_others_rate
    predicted_data[cond,1] <- pred_gray_others_rate
    predicted_data[cond,2] <- pred_gray_color_rate
  }
  
  #chimera
  mean_one_mat <- c(mu83ms_9deg,mu83ms_13deg,mu83ms_17deg,mu83ms_21deg,mu83ms_25deg, mu83ms_9deg*mu117ms,mu83ms_13deg*mu117ms,mu83ms_17deg*mu117ms,mu83ms_21deg*mu117ms,mu83ms_25deg*mu117ms, mu83ms_9deg*mu150ms,mu83ms_13deg*mu150ms,mu83ms_17deg*mu150ms,mu83ms_21deg*mu150ms,mu83ms_25deg*mu150ms)
  sigmamat <- c(rep(sigma83ms,5),rep(sigma117ms,5),rep(sigma150ms,5))
  for (cond in 4:18) {
    pred_chimera_others_rate <- pnorm(cri, mean=mean_one_mat[cond-3], sd=sigmamat[cond-3])
    pred_chimera_color_rate <-1-pred_chimera_others_rate
    #反応数
    pred_nr_chimera_color <- sum(data[cond,1]+data[cond,2]) *pred_chimera_color_rate
    pred_nr_chimera_others <- sum(data[cond,1]+data[cond,2]) * pred_chimera_others_rate
    predicted_data[cond,1] <- pred_chimera_others_rate
    predicted_data[cond,2] <- pred_chimera_color_rate
  }
  
  #color
  mean_two_mat <- c(mu83ms_color,mu83ms_color*mu117ms,mu83ms_color*mu150ms)
  sigmamat <- c(sigma83ms,sigma117ms,sigma150ms)
  for (cond in 19:21) {
    pred_color_others_rate <- pnorm(cri, mean=mean_two_mat[cond-18], sd=sigmamat[cond-18])
    pred_color_color_rate <- 1-pred_color_others_rate
    #反応数
    pred_nr_color_color <- sum(data[cond,1]+data[cond,2]) * pred_color_color_rate
    pred_nr_color_others <- sum(data[cond,1]+data[cond,2]) * pred_color_others_rate
    predicted_data[cond,1] <- pred_color_others_rate
    predicted_data[cond,2] <- pred_color_color_rate
  }
  ##############################################################################
  
  
  est <- data.frame(mu83ms_9deg   = mu83ms_9deg, 
                    mu83ms_13deg  = mu83ms_13deg, 
                    mu83ms_17deg  = mu83ms_17deg, 
                    mu83ms_21deg  = mu83ms_21deg, 
                    mu83ms_25deg  = mu83ms_25deg, 
                    mu83ms_color = mu83ms_color,
                    mu117ms = mu117ms,
                    mu150ms = mu150ms,
                    sigma117ms = sigma117ms,
                    sigma150ms = sigma150ms,
                    cri = cri, 
                    logL = logL,
                    AIC = AIC)
  return(list(est,predicted_data))
}


### Likelihood function
uvsdt_logL <- function(x, inputs) {
  
  # target parameters
  mu83ms_9deg  <- x[1]
  mu83ms_13deg <- x[2] 
  mu83ms_17deg <- x[3] 
  mu83ms_21deg <- x[4] 
  mu83ms_25deg <- x[5] 
  mu83ms_color <- x[6] 
  mu117ms<-  x[7] 
  mu150ms <- x[8] 
  cri <- x[9] 
  
  ##############################################################################
  # model predictions
  predicted_data <- matrix(NA, nrow=21, ncol=2)
  #gray
  mean_one_mat <- c(mu83ms_gray, mu83ms_gray*mu117ms, mu83ms_gray*mu150ms)
  sigmamat <- c(sigma83ms,sigma117ms,sigma150ms)
  for (cond in 1:3) {
    pred_gray_others_rate <- pnorm(cri, mean= mean_one_mat[cond], sd=sigmamat[cond])
    pred_gray_color_rate <- 1-pred_gray_others_rate
    
    #反応数
    pred_nr_gray_color <- sum(data[cond,1]+data[cond,2]) * pred_gray_color_rate
    pred_nr_gray_others <- sum(data[cond,1]+data[cond,2]) * pred_gray_others_rate
    predicted_data[cond,1] <- pred_gray_others_rate
    predicted_data[cond,2] <- pred_gray_color_rate
  }
  
  #chimera
  mean_one_mat <- c(mu83ms_9deg,mu83ms_13deg,mu83ms_17deg,mu83ms_21deg,mu83ms_25deg, mu83ms_9deg*mu117ms,mu83ms_13deg*mu117ms,mu83ms_17deg*mu117ms,mu83ms_21deg*mu117ms,mu83ms_25deg*mu117ms, mu83ms_9deg*mu150ms,mu83ms_13deg*mu150ms,mu83ms_17deg*mu150ms,mu83ms_21deg*mu150ms,mu83ms_25deg*mu150ms)
  sigmamat <- c(rep(sigma83ms,5),rep(sigma117ms,5),rep(sigma150ms,5))
  for (cond in 4:18) {
    pred_chimera_others_rate <- pnorm(cri, mean=mean_one_mat[cond-3], sd=sigmamat[cond-3])
    pred_chimera_color_rate <-1-pred_chimera_others_rate
    #反応数
    pred_nr_chimera_color <- sum(data[cond,1]+data[cond,2]) *pred_chimera_color_rate
    pred_nr_chimera_others <- sum(data[cond,1]+data[cond,2]) * pred_chimera_others_rate
    predicted_data[cond,1] <- pred_chimera_others_rate
    predicted_data[cond,2] <- pred_chimera_color_rate
  }
  
  #color
  mean_two_mat <- c(mu83ms_color,mu83ms_color*mu117ms,mu83ms_color*mu150ms)
  sigmamat <- c(sigma83ms,sigma117ms,sigma150ms)
  for (cond in 19:21) {
    pred_color_others_rate <- pnorm(cri, mean=mean_two_mat[cond-18], sd=sigmamat[cond-18])
    pred_color_color_rate <- 1-pred_color_others_rate
    #反応数
    pred_nr_color_color <- sum(data[cond,1]+data[cond,2]) * pred_color_color_rate
    pred_nr_color_others <- sum(data[cond,1]+data[cond,2]) * pred_color_others_rate
    predicted_data[cond,1] <- pred_color_others_rate
    predicted_data[cond,2] <- pred_color_color_rate
  }
  ##############################################################################
  
  # log likelihood
  logL <- sum(data * log(predicted_data))
  if (is.nan(logL)) {
    logL <- -Inf
  }
  
  logL <- -logL
  return(logL)
  #if (1 < par[3] & par[3] < par[4] & par[4] < par[5] & par[5] < par[6] & par[7] < par[8] & par[8] < par[9] & par[9]< par[10] & par[10]< par[11] & par[12] < par[13] & par[13] < par[14] & par[14] < par[15] & par[15] < par[16]) {
  #logL <- -logL　#対数尤度を対数損失にする
  #return(logL)
  #} 
  #else {
  # return(999999)
  #}
  
}


estimates <- c()
predicted_array <- array(NA, dim = c(21, 2, 44))
data_rate_array <- array(NA, dim = c(21, 2, 44))

for (i in 4:47) {
  
  sub <- i
  
  #データ成型
  subdata <- result[((sub-4)*21+1):((sub-4)*21+21), ]
  data <- matrix(NA, nrow=21, ncol=2)
  #gray
  data[1:3,1]<-unlist(subdata[1:3,5]+subdata[1:3,6])#_others
  data[1:3,2]<-unlist(subdata[1:3,7])#_color
  #chimera
  data[4:18,1]<-unlist(subdata[4:18,8]+subdata[4:18,9])#_others
  data[4:18,2]<-unlist(subdata[4:18,10])#_color
  #color
  data[19:21,1]<-unlist(subdata[19:21,11]+subdata[19:21,12])#_others
  data[19:21,2]<-unlist(subdata[19:21,13])#_color
  
  data_rate_array[,1,(sub-3)]<-data[,1]/(data[,1]+data[,2])
  data_rate_array[,2,(sub-3)]<-data[,2]/(data[,1]+data[,2])

  
  ### Fitting

  fit <- fit_uvsdt_mle(data, add_constant = TRUE)
  df <- fit[[1]]
  df$sub <- i
  estimates <- rbind(estimates, df)
  a <- fit[[2]]
  predicted_array[,,(sub-3)] <- a
}

###################################
#mean+variance,mean,variance
###################################
predicted_array[,,1]
data_rate_array[,,1]


