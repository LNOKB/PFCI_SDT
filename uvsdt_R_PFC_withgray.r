##color, chimera, grayscale #途中
library(tidyverse)

dat <- read.table("subdata.csv", header=T,sep=",") 

result = dat %>%             
  filter(TorF1 == 0 | TorF1 == 1) %>% # dual-task trials
  select(subnum, exposure, maskwidth, colorcond, TorF1, colorres, TorF2, gray_gray, gray_mix, gray_color, mix_gray, mix_mix, mix_color, color_gray, color_mix, color_color) %>%    
  group_by(subnum, colorcond, exposure, maskwidth) %>%                
  summarize(gray_gray= sum(gray_gray), gray_mix= sum(gray_mix), gray_color= sum(gray_color), mix_gray= sum(mix_gray), mix_mix= sum(mix_mix), mix_color= sum(mix_color), color_gray= sum(color_gray), color_mix= sum(color_mix), color_color= sum(color_color)) #%>%       
#print()  

#for (sub in 4:47) {

sub <- 1

#データ成型
subdata <- result[((sub-1)*18+1):((sub-1)*18+18), ]
data <- matrix(NA, nrow=21, ncol=2)
#gray
data[1:3,1]<-unlist(result[1:3,5]+result[1:3,6])#_others
data[1:3,2]<-unlist(result[1:3,7])#_color
#chimera
data[4:18,1]<-unlist(result[4:18,8]+result[4:18,9])#_others
data[4:18,2]<-unlist(result[4:18,10])#_color
#color
data[19:21,1]<-unlist(result[19:21,11]+result[19:21,12])#_others
data[19:21,2]<-unlist(result[19:21,13])#_color

#基準となるパラメータ
sigma83ms <- 1
mu83ms_gray <- 0



### Function for model fitting
fit_uvsdt_mle <- function(data, add_constant = TRUE) {
  
  # add_constant = TRUE adds a small value to the response frequency vectors.
  # correction against extreme estimates 
  if (add_constant) {
    data <- data+0.5
  }
  
  # initial guess for parameter values d-primeベースで出す
  #mu83ms_9deg <- 1.1
  #mu83ms_13deg <- 2
  #mu83ms_17deg <- 3 
  #mu83ms_21deg <- 4
  #mu83ms_25deg <- 5
  #mu83ms_color <- 6
  #mu117ms<-  2
  #mu150ms <- 3
  #sigma117ms <- 1.5
  #sigma150ms <- 2
  #cri <- 1
  
  mu83ms_9deg <-qnorm(data[4,2]/(data[4,1]+data[4,2]))-qnorm(data[1,2]/(data[1,1]+data[1,2]))
  mu83ms_13deg <-qnorm(data[5,2]/(data[5,1]+data[5,2]))-qnorm(data[1,2]/(data[1,1]+data[1,2]))
  mu83ms_17deg <- qnorm(data[6,2]/(data[6,1]+data[6,2]))-qnorm(data[1,2]/(data[1,1]+data[1,2]))
  mu83ms_21deg <- qnorm(data[7,2]/(data[7,1]+data[7,2]))-qnorm(data[1,2]/(data[1,1]+data[1,2]))
  mu83ms_25deg <- qnorm(data[8,2]/(data[8,1]+data[8,2]))-qnorm(data[1,2]/(data[1,1]+data[1,2]))
  mu83ms_color <- qnorm(data[19,2]/(data[19,1]+data[19,2]))-qnorm(data[1,2]/(data[1,1]+data[1,2]))
  mu117ms<-  mu83ms_9deg*2
  mu150ms <- mu83ms_9deg*3
  sigma117ms <- sigma83ms*1.5
  sigma150ms <- sigma83ms*2
  cri <- 1
  
  guess <- c(mu83ms_9deg, mu83ms_13deg, mu83ms_17deg, mu83ms_21deg, mu83ms_25deg, mu83ms_color, mu117ms, mu150ms, sigma117ms, sigma150ms, cri)
  
  # model fit
  fit <- suppressWarnings(optim(uvsdt_logL, 
                                par = guess, 
                                inputs = list("gray_others" = data[1:3,1],  "gray_color" = data[1:3,2], "chimera_others" = data[4:18,1],  "chimera_color" = data[4:18,2], "color_others" = data[19:21,1], "color_color" = data[19:21,2]), 
                                gr = NULL, method = "BFGS", control = list("maxit" = 10000)))
  
  # outputs
  mu83ms_9deg <- fit$par[1]
  mu83ms_13deg <- fit$par[2] 
  mu83ms_17deg <- fit$par[3] 
  mu83ms_21deg <- fit$par[4] 
  mu83ms_25deg <- fit$par[5] 
  mu83ms_color <- fit$par[6] 
  mu117ms<-  fit$par[7] 
  mu150ms <- fit$par[8] 
  sigma117ms <- fit$par[9] 
  sigma150ms <- fit$par[10]
  cri <- fit$par[11] 
  
  logL <- -fit$value
  
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
                    logL = logL)
  return(est)
}


### Likelihood function
uvsdt_logL <- function(x, inputs) {
  
  # target parameters
  mu83ms_9deg <- x[1]
  mu83ms_13deg <- x[2] 
  mu83ms_17deg <- x[3] 
  mu83ms_21deg <- x[4] 
  mu83ms_25deg <- x[5] 
  mu83ms_color <- x[6] 
  mu117ms<-  x[7] 
  mu150ms <- x[8] 
  sigma117ms <- x[9] 
  sigma150ms <- x[10]
  cri <- x[11] 
  
  # model predictions
  predicted_data <- matrix(NA, nrow=21, ncol=2)
  #gray
  mean_one_mat <- c(mu83ms_gray, mu83ms_gray*mu117ms, mu83ms_gray*mu150ms)
  sigmamat <- c(sigma83ms,sigma117ms,sigma150ms)
  for (cond in 1:3) {
    mean_inverted <- -mean_one_mat[cond]+2*cri
    pred_gray_color_rate <- pnorm(cri, mean=mean_inverted, sd=sigmamat[cond])
    pred_gray_others_rate <- 1-pred_gray_color_rate
    #反応数
    pred_nr_gray_color <- sum(data[cond,1]+data[cond,2]) * pred_gray_color_rate
    pred_nr_gray_others <- sum(data[cond,1]+data[cond,2]) * pred_gray_others_rate
    predicted_data[cond,2] <-pred_nr_gray_color
    predicted_data[cond,1] <- pred_nr_gray_others
  }
  
  #chimera
  mean_one_mat <- c(mu83ms_9deg,mu83ms_13deg,mu83ms_17deg,mu83ms_21deg,mu83ms_25deg, mu83ms_9deg*mu117ms,mu83ms_13deg*mu117ms,mu83ms_17deg*mu117ms,mu83ms_21deg*mu117ms,mu83ms_25deg*mu117ms, mu83ms_9deg*mu150ms,mu83ms_13deg*mu150ms,mu83ms_17deg*mu150ms,mu83ms_21deg*mu150ms,mu83ms_25deg*mu150ms)
  sigmamat <- c(rep(sigma83ms,5),rep(sigma117ms,5),rep(sigma150ms,5))
  for (cond in 4:18) {
    mean_inverted <- -mean_one_mat[cond-3]+2*cri
    pred_chimera_color_rate <- pnorm(cri, mean=mean_inverted, sd=sigmamat[cond-3])
    pred_chimera_others_rate <- 1-pred_chimera_color_rate
    #反応数
    pred_nr_chimera_color <- sum(data[cond,1]+data[cond,2]) *pred_chimera_color_rate
    pred_nr_chimera_others <- sum(data[cond,1]+data[cond,2]) * pred_chimera_others_rate
    predicted_data[cond,2] <- pred_nr_chimera_color
    predicted_data[cond,1] <- pred_nr_chimera_others
  }
  
  #color
  mean_two_mat <- c(mu83ms_color,mu83ms_color*mu117ms,mu83ms_color*mu150ms)
  sigmamat <- c(sigma83ms,sigma117ms,sigma150ms)
  for (cond in 19:21) {
    mean_inverted <- -mean_one_mat[cond-18]+2*cri
    pred_color_color_rate <- pnorm(cri, mean=mean_inverted, sd=sigmamat[cond-18])
    pred_color_others_rate <- 1-pred_color_color_rate
    #反応数
    pred_nr_color_color <- sum(data[cond,1]+data[cond,2]) * pred_color_color_rate
    pred_nr_color_others <- sum(data[cond,1]+data[cond,2]) * pred_color_others_rate
    predicted_data[cond,2] <-pred_nr_color_color
    predicted_data[cond,1] <- pred_nr_color_others
  }
  
  
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

### Fitting
fit <- fit_uvsdt_mle(data, add_constant = TRUE)
fit

###################################
#mean+variance,mean,variance
###################################

