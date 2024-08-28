##color, chimera, grayscale #途中
library(tidyverse)

dat <- read.table("subdata.csv", header=T,sep=",") 

result = dat %>%             
  filter(TorF1 == 0 | TorF1 == 1) %>% # dual-task trials
  select(subnum, exposure, maskwidth, colorcond, TorF1, colorres, TorF2, gray_gray, gray_color, mix_mix, mix_color, color_gray, color_mix, color_color) %>%    
  group_by(subnum, colorcond, exposure, maskwidth) %>%                
  summarize(gray_gray= sum(gray_gray), gray_color= sum(gray_color), mix_mix= sum(mix_mix), mix_color= sum(mix_color), color_gray= sum(color_gray), color_mix= sum(color_mix), color_color= sum(color_color)) #%>%       
  #print()  

#for (sub in 4:47) {

sub <- 1

sigma83ms <- 1
mu83ms_gray <- 0

### Function for model fitting
fit_uvsdt_mle <- function(data, add_constant = TRUE) {
  
  # add_constant = TRUE adds a small value to the response frequency vectors.
  # correction against extreme estimates 
  constant <- 0
  if (add_constant) {
    constant <- 0.5
  }
  data<-cbind(result[,1:4],result[,5:10]+constant)
  
  # initial guess for parameter values d-primeベースで出す
  a <- 1.1
  b <- 1.2
  c <- 1
  d <- 2
  e <- 3
  f <- 4
  g <- 1
  h <- 2
  i <- 3
  j <- 4
  k <- 5
  l <- 1
  m <- 2
  n <- 3
  o <- 4
  p <- 5
  q <- 6
  r <- 6
  s <- 6
  cri <- 1 
  
  guess <- c(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, cri)
  
  # model fit
  fit <- suppressWarnings(optim(uvsdt_logL, 
                                par = guess, 
                                inputs = list("gray_gray" = data$gray_gray, "gray_color" = data$gray_color, "mix_mix" = data$mix_mix, "mix_color" = data$mix_color, "color_mix" = data$color_mix, "color_color" = data$color_color), 
                                gr = NULL, method = "BFGS", control = list("maxit" = 10000)))
  
  ########################################################
  #a <- 1.1
  #b <- 1.2
  #c <- 1
  #d <- 2
  #e <- 3
  #f <- 4
  #g <- 1
  #h <- 2
  #i <- 3
  #j <- 4
  #k <- 5
  #l <- 1
  #m <- 2
  #n <- 3
  #o <- 4
  #p <- 5
  #q <- 6
  #r <- 6
  #s <- 6
  #fit$par <- c(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, cri)
  ########################################################
  
  
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
  
  #反応率　criterionより上側の面積を出す
  predicted_data <- matrix(NA, nrow=20, ncol=2)
  
  #gray&color    gray_gray, gray_color color_gray, , 
  mean_one_mat <- c(mu83ms_gray, mu83ms_gray*mu117ms, mu83ms_gray*mu117ms)
  sigmamat <- c(sigma83ms,sigma117ms,sigma150ms)
  for (cond in 1:3) {
    pred_gray_gray_rate <- pnorm(cri, mean=mean_one_mat[cond], sd=sigmamat[cond])
    pred_gray_color_rate <- 1-pred_chimera_chimera_rate
    #反応数
    pred_nr_gray_gray <- sum(data$gray_gray[cond]+data$gray_color[cond]) * pred_gray_gray_rate
    pred_nr_gray_color <- sum(data$gray_gray[cond]+data$gray_color[cond]) * pred_gray_color_rate
    predicted_data[cond,1] <-pred_nr_gray_gray
    predicted_data[cond,2] <- pred_nr_gray_color
  }
  
  
  #chimera&color
  mean_one_mat <- c(mu83ms_9deg,mu83ms_13deg,mu83ms_17deg,mu83ms_21deg,mu83ms_25deg, mu83ms_9deg*mu117ms,mu83ms_13deg*mu117ms,mu83ms_17deg*mu117ms,mu83ms_21deg*mu117ms,mu83ms_25deg*mu117ms, mu83ms_9deg*mu150ms,mu83ms_13deg*mu150ms,mu83ms_17deg*mu150ms,mu83ms_21deg*mu150ms,mu83ms_25deg*mu150ms)
  sigmamat <- c(rep(sigma83ms,5),rep(sigma117ms,5),rep(sigma150ms,5))
  for (cond in 4:18) {
    pred_chimera_chimera_rate <- pnorm(cri, mean=mean_one_mat[cond-3], sd=sigmamat[cond-3])
    pred_chimera_color_rate <- 1-pred_chimera_chimera_rate
    #反応数
    pred_nr_chimera_chimera <- sum(data$chimera_chimera[cond]+data$chimera_color[cond]) *pred_chimera_chimera_rate
    pred_nr_chimera_color <- sum(data$chimera_chimera[cond]+data$chimera_color[cond]) * pred_chimera_color_rate
    predicted_data[cond,1] <- pred_nr_chimera_chimera
    predicted_data[cond,2] <- pred_nr_chimera_color
  }
  
  #color&chimera
  mean_two_mat <- c(mu83ms_color,mu83ms_color*mu117ms,mu83ms_color*mu150ms)
  sigmamat <- c(sigma83ms,sigma117ms,sigma150ms)
  for (cond in 19:21) {
    pred_color_chimera_rate <- pnorm(cri, mean=mean_two_mat[cond-18], sd=sigmamat[cond-18])
    pred_color_color_rate <- 1-pred_color_chimera_rate
    #反応数
    pred_nr_color_chimera <- sum(data$color_chimera[cond]+data$color_color[cond]) * pred_color_chimera_rate
    pred_nr_color_color <- sum(data$color_chimera[cond]+data$color_color[cond]) * pred_color_color_rate
    predicted_data[cond,1] <-pred_nr_color_chimera
    predicted_data[cond,2] <- pred_nr_color_color
  }
     

  
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
  return(list(est,predicted_data))
#(est)
  
}

#}

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
  
  # empirical data
  hit <- inputs$hit
  miss <- inputs$miss 
  cr <- inputs$cr
  fa <- inputs$fa
  
  # model predictions
  #反応率　criterionより上側の面積を出す
  sigmamat <- c(rep(sigma83ms,5),rep(sigma117ms,5),rep(sigma150ms,5))
  mean_one_mat <- c(mu83ms_9deg,mu83ms_13deg,mu83ms_17deg,mu83ms_21deg,mu83ms_25deg,mu117ms_9deg,mu117ms_13deg,mu117ms_17deg,mu117ms_21deg,mu117ms_25deg,mu150ms_9deg,mu150ms_13deg,mu150ms_17deg,mu150ms_21deg,mu150ms_25deg)
  mean_two_mat <- c(rep(mu83ms_color,5),rep(mu117ms_color,5),rep(mu150ms_color,5))
  
  predicted_data <- matrix(NA, nrow=15, ncol=4)
  
  for (cond in 1:15) {
    
    pred_missr <- pnorm(cri, mean=mean_two_mat[cond], sd=sigmamat[cond])
    pred_hitr <- 1-pred_missr
    pred_crr <- pnorm(cri, mean=mean_one_mat[cond], sd=sigmamat[cond])
    pred_far <- 1-pred_crr
    
    #反応数
    pred_nr_miss <- sum(hit+miss) * pred_missr
    pred_nr_hit <- sum(hit+miss) * pred_hitr
    pred_nr_cr <- sum(cr+fa) * pred_crr
    pred_nr_fa <- sum(cr+fa) * pred_far
    
    predicted_data[cond,1] <- pred_hitr
    predicted_data[cond,2] <- pred_missr
    predicted_data[cond,3] <- pred_crr
    predicted_data[cond,4] <- pred_far
  }
  
  # log likelihood ここがまわってない
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
fit <- fit_uvsdt_mle(data, add_constant = FALSE)
fit

###################################
#mean+variance,mean,variance
###################################

