##color, chimera
library(tidyverse)

dat <- read.table("subdata.csv", header=T,sep=",") 

result = dat %>%             
  filter(TorF1 == 0 | TorF1 == 1) %>% # dual-task trials
  filter(colorcond == 1 | colorcond == 2) %>%  # chimera, full-color trials
  select(subnum, exposure, maskwidth, colorcond, TorF1, colorres, TorF2, mix_mix, mix_color, color_mix, color_color) %>%    
  group_by(subnum, exposure, colorcond, maskwidth) %>%                
  summarize(hit= sum(color_color), miss= sum(color_mix), cr= sum(mix_mix), fa= sum(mix_color)) %>%       
  print()  

#for (sub in 4:47) {

sub <- 1

#データ成型
subresult <- result[((sub-1)*18+1):((sub-1)*18+18), ]
subresult[1:5,5:6]　<- subresult[6,5:6] 
subresult[7:11,5:6]　<- subresult[12,5:6] 
subresult[13:17,5:6]　<- subresult[18,5:6] 
data <- subresult[c(-6,-12,-18),]

sigma83ms <- 1
mu83ms_9deg <- 0

### Function for model fitting
fit_uvsdt_mle <- function(data, add_constant = TRUE) {
  
  # add_constant = TRUE adds a small value to the response frequency vectors.
  # correction against extreme estimates 
  if (add_constant) {
    data$hit <- data$hit + 0.5
    data$miss <- data$miss + 0.5
    data$cr <- data$cr + 0.5
    data$fa <- data$fa + 0.5
  }
  
  # initial guess for parameter values
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
                                inputs = list("hit" = data$hit, "miss" = data$miss, "cr" = data$cr, "fa" = data$fa), 
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
  sigma117ms <- fit$par[1] * sigma83ms
  sigma150ms <- fit$par[2] * sigma83ms
  
  mu83ms_13deg <- fit$par[3] #* mu83ms_9deg
  mu83ms_17deg <- fit$par[4] #* mu83ms_9deg
  mu83ms_21deg <- fit$par[5] #* mu83ms_9deg
  mu83ms_25deg <- fit$par[6] #* mu83ms_9deg
  
  mu117ms_9deg <-  fit$par[7] #* mu83ms_9deg
  mu117ms_13deg <- fit$par[8] #* mu83ms_9deg
  mu117ms_17deg <- fit$par[9] #* mu83ms_9deg
  mu117ms_21deg <- fit$par[10] #* mu83ms_9deg
  mu117ms_25deg <- fit$par[11] #* mu83ms_9deg
  
  mu150ms_9deg <-  fit$par[12] #* mu83ms_9deg
  mu150ms_13deg <- fit$par[13] #* mu83ms_9deg
  mu150ms_17deg <- fit$par[14] #* mu83ms_9deg
  mu150ms_21deg <- fit$par[15] #* mu83ms_9deg
  mu150ms_25deg <- fit$par[16] #* mu83ms_9deg
  
  mu83ms_color <- fit$par[17] #* mu83ms_9deg
  mu117ms_color <- fit$par[18] #* mu83ms_9deg
  mu150ms_color <- fit$par[19] #* mu83ms_9deg
  
  cri <- fit$par[20] 
  
  logL <- -fit$value
  
  #反応率　criterionより上側の面積を出す
  sigmamat <- c(rep(sigma83ms,5),rep(sigma117ms,5),rep(sigma150ms,5))
  mean_one_mat <- c(mu83ms_9deg,mu83ms_13deg,mu83ms_17deg,mu83ms_21deg,mu83ms_25deg,mu117ms_9deg,mu117ms_13deg,mu117ms_17deg,mu117ms_21deg,mu117ms_25deg,mu150ms_9deg,mu150ms_13deg,mu150ms_17deg,mu150ms_21deg,mu150ms_25deg)
  mean_two_mat <- c(rep(mu83ms_color,5),rep(mu117ms_color,5),rep(mu150ms_color,5))
  
  predicted_data2 <- matrix(NA, nrow=15, ncol=4)
  for (cond in 1:15) {
    
    pred_missr <- pnorm(cri, mean=mean_two_mat[cond], sd=sigmamat[cond])
    pred_hitr <- 1-pred_missr
    pred_crr <- pnorm(cri, mean=mean_one_mat[cond], sd=sigmamat[cond])
    pred_far <- 1-pred_crr
    
    #反応数
    pred_nr_miss <- sum(data$hit[cond]+data$miss[cond]) * pred_missr
    pred_nr_hit <- sum(data$hit[cond]+data$miss[cond]) * pred_hitr
    pred_nr_cr <- sum(data$cr[cond]+data$fa[cond]) * pred_crr
    pred_nr_fa <- sum(data$cr[cond]+data$fa[cond]) * pred_far
    
    predicted_data2[cond,1] <-pred_nr_hit
    predicted_data2[cond,2] <- pred_nr_miss
    predicted_data2[cond,3] <- pred_nr_cr
    predicted_data2[cond,4] <- pred_nr_fa
  }
  
  
  
  est <- data.frame(sigma117ms = sigma117ms,
                    sigma150ms = sigma150ms,
                    mu83ms_9deg   = mu83ms_9deg, 
                    mu83ms_13deg  = mu83ms_13deg, 
                    mu83ms_17deg  = mu83ms_17deg, 
                    mu83ms_21deg  = mu83ms_21deg, 
                    mu83ms_25deg  = mu83ms_25deg, 
                    mu117ms_9deg  = mu117ms_9deg, 
                    mu117ms_13deg = mu117ms_13deg, 
                    mu117ms_17deg = mu117ms_17deg, 
                    mu117ms_21deg = mu117ms_21deg, 
                    mu117ms_25deg = mu117ms_25deg, 
                    mu150ms_9deg  = mu150ms_9deg, 
                    mu150ms_13deg = mu150ms_13deg, 
                    mu150ms_17deg = mu150ms_17deg,
                    mu150ms_21deg = mu150ms_21deg, 
                    mu150ms_25deg = mu150ms_25deg,   
                    mu83ms_color = mu83ms_color,
                    mu117ms_color = mu117ms_color,
                    mu150ms_color = mu150ms_color,
                    cri = cri, 
                    logL = logL)
  return(predicted_data2)
#(est)
  
}

#}

### Likelihood function
uvsdt_logL <- function(x, inputs) {
  
  # target parameters
  
  sigma117ms <- sigma83ms * x[1]
  sigma150ms <- sigma83ms * x[2]
  
  mu83ms_13deg <-  x[3]
  mu83ms_17deg <-  x[4]
  mu83ms_21deg <-  x[5]
  mu83ms_25deg <-  x[6]
  
  mu117ms_9deg <-  x[7] 
  mu117ms_13deg <-  x[8]
  mu117ms_17deg <-  x[9]
  mu117ms_21deg <-  x[10] 
  mu117ms_25deg <-  x[11]
  
  mu150ms_9deg <-  x[12]
  mu150ms_13deg <-  x[13]
  mu150ms_17deg <-  x[14]
  mu150ms_21deg <-  x[15]
  mu150ms_25deg <-  x[16]
  
  mu83ms_color <-  x[17]
  mu117ms_color <-  x[18]
  mu150ms_color <-  x[19]
  
  cri <- x[20] 
  
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

