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
  sigma117ms <- 1
  sigma150ms <- 1
  cri <- 2
  
  guess <- c(mu83ms_9deg, mu83ms_13deg, mu83ms_17deg, mu83ms_21deg, mu83ms_25deg, mu83ms_color, mu117ms, mu150ms, sigma117ms, sigma150ms, cri)
  
  # model fit
  fit <- suppressWarnings(optim(uvsdt_logL, 
                                par = guess, 
                                lower =      c(0.1, 0.1, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5),
                                upper =      c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 2.0, 2.0, 2.0, 2.0, 5.0), 
                                gr = NULL, method = "BFGS", 
                                control = list("maxit" = 10000, 
                                               "parscale" = c(1,   1,    1,   1,    1,   1,  0.001,  0.001,  0.001,  0.001,  1))))
  
  # outputs
  mu83ms_9deg  <- fit$par[1]
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
  
  AIC <- -2 * logL + 2 * 11
  
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
  sigma117ms <- x[9] 
  sigma150ms <- x[10]
  cri <- x[11] 
  
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

data_rate_add <- array(NA, dim = c(21, 2, 44))
dprime_mat <- matrix(c(rep(NA,18*44)), nrow=18, ncol=44)
C_mat <- matrix(c(rep(NA,18*44)), nrow=18, ncol=44)

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
  
  data_rate_add[,1,(sub-3)]<-(data[,1]+0.5)/((data[,1]+0.5)+(data[,2]+0.5))
  data_rate_add[,2,(sub-3)]<-(data[,2]+0.5)/((data[,1]+0.5)+(data[,2]+0.5))
  
  #########カラーをS2分布としたd'やCの算出[グレーをS1分布とすべき？]############################################################
  
  dprime_gr_83 <-  qnorm(data_rate_add[19,2,(sub-3)])-qnorm(data_rate_add[1,2,(sub-3)])
  dprime_gr_117 <- qnorm(data_rate_add[20,2,(sub-3)])-qnorm(data_rate_add[2,2,(sub-3)])
  dprime_gr_150 <- qnorm(data_rate_add[21,2,(sub-3)])-qnorm(data_rate_add[3,2,(sub-3)])
  
  dprime_9d_83<-  qnorm(data_rate_add[19,2,(sub-3)])-qnorm(data_rate_add[4,2,(sub-3)])
  dprime_13d_83<- qnorm(data_rate_add[19,2,(sub-3)])-qnorm(data_rate_add[5,2,(sub-3)])
  dprime_17d_83<- qnorm(data_rate_add[19,2,(sub-3)])-qnorm(data_rate_add[6,2,(sub-3)])
  dprime_21d_83<- qnorm(data_rate_add[19,2,(sub-3)])-qnorm(data_rate_add[7,2,(sub-3)])
  dprime_25d_83<- qnorm(data_rate_add[19,2,(sub-3)])-qnorm(data_rate_add[8,2,(sub-3)])
  
  dprime_9d_117<-  qnorm(data_rate_add[20,2,(sub-3)])-qnorm(data_rate_add[9,2,(sub-3)])
  dprime_13d_117<- qnorm(data_rate_add[20,2,(sub-3)])-qnorm(data_rate_add[10,2,(sub-3)])
  dprime_17d_117<- qnorm(data_rate_add[20,2,(sub-3)])-qnorm(data_rate_add[11,2,(sub-3)])
  dprime_21d_117<- qnorm(data_rate_add[20,2,(sub-3)])-qnorm(data_rate_add[12,2,(sub-3)])
  dprime_25d_117<- qnorm(data_rate_add[20,2,(sub-3)])-qnorm(data_rate_add[13,2,(sub-3)])
  
  dprime_9d_150<-  qnorm(data_rate_add[21,2,(sub-3)])-qnorm(data_rate_add[14,2,(sub-3)])
  dprime_13d_150<- qnorm(data_rate_add[21,2,(sub-3)])-qnorm(data_rate_add[15,2,(sub-3)])
  dprime_17d_150<- qnorm(data_rate_add[21,2,(sub-3)])-qnorm(data_rate_add[16,2,(sub-3)])
  dprime_21d_150<- qnorm(data_rate_add[21,2,(sub-3)])-qnorm(data_rate_add[17,2,(sub-3)])
  dprime_25d_150<- qnorm(data_rate_add[21,2,(sub-3)])-qnorm(data_rate_add[18,2,(sub-3)])
  
  C_gr_83 <- -0.5 * (qnorm(data_rate_add[1,2,(sub-3)])+qnorm(data_rate_add[19,2,(sub-3)]))
  C_gr_117 <- -0.5 * (qnorm(data_rate_add[2,2,(sub-3)])+qnorm(data_rate_add[20,2,(sub-3)]))
  C_gr_150 <- -0.5 * (qnorm(data_rate_add[3,2,(sub-3)])+qnorm(data_rate_add[21,2,(sub-3)]))
  
  C_9d_83<- -0.5 * (qnorm(data_rate_add[4,2,(sub-3)])+qnorm(data_rate_add[19,2,(sub-3)]))
  C_13d_83<- -0.5 * (qnorm(data_rate_add[5,2,(sub-3)])+qnorm(data_rate_add[19,2,(sub-3)]))
  C_17d_83<- -0.5 * (qnorm(data_rate_add[6,2,(sub-3)])+qnorm(data_rate_add[19,2,(sub-3)]))
  C_21d_83<- -0.5 * (qnorm(data_rate_add[7,2,(sub-3)])+qnorm(data_rate_add[19,2,(sub-3)]))
  C_25d_83<- -0.5 * (qnorm(data_rate_add[8,2,(sub-3)])+qnorm(data_rate_add[19,2,(sub-3)]))
  
  C_9d_117<- -0.5 * (qnorm(data_rate_add[9,2,(sub-3)])+qnorm(data_rate_add[20,2,(sub-3)]))
  C_13d_117<- -0.5 * (qnorm(data_rate_add[10,2,(sub-3)])+qnorm(data_rate_add[20,2,(sub-3)]))
  C_17d_117<- -0.5 * (qnorm(data_rate_add[11,2,(sub-3)])+qnorm(data_rate_add[20,2,(sub-3)]))
  C_21d_117<- -0.5 * (qnorm(data_rate_add[12,2,(sub-3)])+qnorm(data_rate_add[20,2,(sub-3)]))
  C_25d_117<- -0.5 * (qnorm(data_rate_add[13,2,(sub-3)])+qnorm(data_rate_add[20,2,(sub-3)]))
  
  C_9d_150<- -0.5 * (qnorm(data_rate_add[14,2,(sub-3)])+qnorm(data_rate_add[21,2,(sub-3)]))
  C_13d_150<- -0.5 * (qnorm(data_rate_add[15,2,(sub-3)])+qnorm(data_rate_add[21,2,(sub-3)]))
  C_17d_150<- -0.5 * (qnorm(data_rate_add[16,2,(sub-3)])+qnorm(data_rate_add[21,2,(sub-3)]))
  C_21d_150<- -0.5 * (qnorm(data_rate_add[17,2,(sub-3)])+qnorm(data_rate_add[21,2,(sub-3)]))
  C_25d_150<- -0.5 * (qnorm(data_rate_add[18,2,(sub-3)])+qnorm(data_rate_add[21,2,(sub-3)]))
  
  dprime <- c(dprime_gr_83,dprime_gr_117,dprime_gr_150,dprime_9d_83,dprime_13d_83,dprime_17d_83,dprime_21d_83,dprime_25d_83,dprime_9d_117,dprime_13d_117,dprime_17d_117,dprime_21d_117,dprime_25d_117,dprime_9d_150,dprime_13d_150,dprime_17d_150,dprime_21d_150,dprime_25d_150)
  C <- c(C_gr_83,C_gr_117,C_gr_150,C_9d_83,C_13d_83,C_17d_83,C_21d_83,C_25d_83,C_9d_117,C_13d_117,C_17d_117,C_21d_117,C_25d_117,C_9d_150,C_13d_150,C_17d_150,C_21d_150,C_25d_150)
  
  dprime_mat[,(sub-3)] <- dprime
  C_mat[,(sub-3)] <- C
  ################################################################################
  

  
  ### Fitting

  fit <- fit_uvsdt_mle(data, add_constant = TRUE)
  df <- fit[[1]]
  df$sub <- i
  estimates <- rbind(estimates, df)
  a <- fit[[2]]
  predicted_array[,,(sub-3)] <- a
}

#predicted_array[,,1]
#data_rate_array[,,1]

mean_predicted <- apply(predicted_array, c(1, 2), mean)*100
mean_predicted2 <- c(mean_predicted[1,2],mean_predicted[4:8,2],mean_predicted[19,2],mean_predicted[2,2],mean_predicted[9:13,2],mean_predicted[20,2],mean_predicted[3,2],mean_predicted[14:18,2],mean_predicted[21,2])
mean_data <-  apply(data_rate_array, c(1, 2), mean)*100
mean_data2 <- rbind(mean_data[1,],mean_data[4:8,],mean_data[19,],mean_data[2,],mean_data[9:13,],mean_data[20,],mean_data[3,],mean_data[14:18,],mean_data[21,])

################################################################################mean+variance,mean,variance,nullの検定
for (i in 1:4){
  now <- log(estimates[,i+6])
  #now <- estimates[,i+6]
  t_test_below1 <- t.test(now, mu = 0, alternative = "less")
  t_test_above1 <- t.test(now, mu = 0, alternative = "greater")
  print(t_test_below1)
  print(t_test_above1)
}

library(ggplot2)
data_parameter_plot <- data.frame(
  Parameters = rep(c("λ117ms", "λ150ms", "σ117ms", "σ150ms"), each = 44),
  Value = c(estimates[,7], estimates[,8], estimates[,9], estimates[,10])
)
data_parameter_plot$Value <- log(data_parameter_plot$Value)

# グラフの作成
ggplot(data_parameter_plot, aes(x = Parameters, y = Value)) +
  geom_violin(fill = "skyblue", color = "black") +  # ヴァイオリンプロットの描画
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # y = 1 の位置に横線を追加
  labs(x = "Parameters",
       y = "Value") +
  coord_cartesian(ylim = c(-1, 1)) +  # y軸の範囲を0から1.2に)設定
  stat_summary(fun = mean, geom = "point", 
               shape =16, size = 2, color = "black")+
  theme_classic() +  # テーマの設定
  theme(
    plot.title = element_text(size = 20 * 2),    # タイトルのサイズを5倍に
    axis.title.x = element_text(size = 14 * 2),  # x軸ラベルのサイズを5倍に
    axis.title.y = element_text(size = 14 * 2),  # y軸ラベルのサイズを5倍に
    axis.text.x = element_text(size = 11 * 2),   # x軸目盛りのサイズを5倍に
    axis.text.y = element_text(size = 11 * 2)    # y軸目盛りのサイズを5倍に
  )



################################################################################ 参加者ごとのSDTの図示
pickup_sub <- 1

plot_sdt_distributions <- function(means, sds, attention_levels, image_types, colors) {
  data_SDT_plot <- data.frame()
  
  for (i in 1:length(attention_levels)) {
    for (j in 1:length(image_types)) {
      x <- seq(means[i,j] - 3 * sds[i,j], means[i,j] + 3 * sds[i,j], length.out = 100)
      y <- dnorm(x, mean = means[i,j], sd = sds[i,j])

      
      # データフレームに追加
      data_SDT_plot <- rbind(data_SDT_plot, data.frame(
        x = x,
        y = y,
        Attention = attention_levels[i],
        ImageType = image_types[j],
        color = colors[j]
      ))
    }
  }
  
  p <- ggplot(data_SDT_plot, aes(x = x, y = y, color = ImageType)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = colors) +
    labs(x = "Signal strength",
         y = "Probability Density") +
    scale_y_continuous(breaks=seq(0,1,length=5),limits=c(0,1))+
    geom_vline(xintercept = estimates[pickup_sub,11], linetype = "dashed", color = "black") + 
    facet_wrap(~ Attention, nrow = 3, scales = "free_y") +  # 各注意条件で分布を重ねる
    theme_minimal() +
    theme(legend.position = "top")
  
  print(p)
}

# 平均値と標準偏差
means_83ms <- c(mu83ms_gray, estimates[pickup_sub,1], estimates[pickup_sub,2], estimates[pickup_sub,3], estimates[pickup_sub,4], estimates[pickup_sub,5], estimates[pickup_sub,6])  
means_117ms <- means_83ms * estimates[pickup_sub,7]
means_150ms <- means_83ms * estimates[pickup_sub,8]
means <- rbind(means_83ms,means_117ms,means_150ms)
sd_83ms <- c(rep(sigma83ms, 7))  
sd_117ms <- c(rep(estimates[pickup_sub,9], 7))  
sd_150ms <- c(rep(estimates[pickup_sub,10], 7))  
sds <- rbind(sd_83ms,sd_117ms,sd_150ms)

attention_levels <- c("0_83ms", "1_117ms", "2_150ms")
image_types <- c("0_gray", "1_9deg", "2_13deg", "3_17deg", "4_21deg", "5_25deg", "6_color")
colors <- c("grey", "mistyrose", "pink", "salmon", "orangered2","red", "brown")

# 関数を呼び出し
plot_sdt_distributions(means, sds, attention_levels, image_types, colors)



################################################################################ 参加者平均SDTの図示

plot_sdt_distributions <- function(means, sds, attention_levels, image_types, colors) {
  data_SDT_plot <- data.frame()
  
  for (i in 1:length(attention_levels)) {
    for (j in 1:length(image_types)) {
      x <- seq(means[i,j] - 3 * sds[i,j], means[i,j] + 3 * sds[i,j], length.out = 100)
      y <- dnorm(x, mean = means[i,j], sd = sds[i,j]) #最大

      
      # データフレームに追加
      data_SDT_plot <- rbind(data_SDT_plot, data.frame(
        x = x,
        y = y,
        Attention = attention_levels[i],
        ImageType = image_types[j],
        color = colors[j]
      ))
    }
  }
  
  p <- ggplot(data_SDT_plot, aes(x = x, y = y, color = ImageType)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = colors) +
    labs(x = "Signal strength",
         y = "Probability Density") +
    scale_y_continuous(breaks=seq(0,1,length=5),limits=c(0,1))+
    geom_vline(xintercept = estimates[pickup_sub,11], linetype = "dashed", color = "black") + 
    facet_wrap(~ Attention, nrow = 3, scales = "free_y") +  # 各注意条件で分布を重ねる
    theme_minimal() +
    theme(legend.position = "top")
  
  print(p)
}

# 平均値と標準偏差
 means_83ms <- c(mu83ms_gray, mean(estimates[,1]), mean(estimates[,2]), mean(estimates[,3]), mean(estimates[,4]), mean(estimates[,5]),mean(estimates[,6]))  
 means_117ms <- means_83ms * mean(estimates[,7])
 means_150ms <- means_83ms * mean(estimates[,8])
 means <- rbind(means_83ms,means_117ms,means_150ms)
 sd_83ms <- c(rep(sigma83ms, 7))  
 sd_117ms <- c(rep(mean(estimates[,9]), 7))  
 sd_150ms <- c(rep(mean(estimates[,10]), 7))  
 sds <- rbind(sd_83ms,sd_117ms,sd_150ms)
 
attention_levels <- c("0_83ms", "1_117ms", "2_150ms")
image_types <- c("0_gray", "1_9deg", "2_13deg", "3_17deg", "4_21deg", "5_25deg", "6_color")
colors <- c("grey", "mistyrose", "pink", "salmon", "orangered2","red", "brown")

# 関数を呼び出し
plot_sdt_distributions(means, sds, attention_levels, image_types, colors)



#########積み上げグラフでの比較######################################
# のデータを作成
data_stackedbar <- data.frame(
  Degree = factor(rep(c("gray","9deg", "13deg", "17deg", "21deg", "25deg", "color"), each = 2), 
                  levels = c("gray", "9deg", "13deg", "17deg", "21deg", "25deg", "color")),
  Type = rep(c("0_others","1_full-color"), 7),
  Frequency =  as.vector(t(mean_data2)),
  Condition = factor(rep(c("83ms", "117ms", "150ms"), each = 14), 
                     levels = c("83ms", "117ms", "150ms"))
)

# 予測値のデータフレームを作成（仮の予測値）
predicted_data_stackedbar <- data.frame(
  Degree = factor(c("gray","9deg", "13deg", "17deg", "21deg", "25deg", "color"),
                  levels = c("gray", "9deg", "13deg", "17deg", "21deg", "25deg", "color")),
  Type = rep(c("1_full-color"), 7),
  Predicted = as.vector(t(mean_predicted2)),  # 仮の予測値
  Condition = factor(rep(c("83ms", "117ms", "150ms"), each = 7), 
                     levels = c("83ms", "117ms", "150ms"))
)

# グラフの作成
ggplot(data_stackedbar, aes(x = Degree, y = Frequency, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = NULL, y = "Response proportion(%)") +
  theme_minimal() +
  scale_fill_manual(values = c("1_full-color" = "gray50","0_others" = "gray80" )) +
  facet_grid(. ~ Condition) +  
  geom_point(data = predicted_data_stackedbar, aes(x = Degree, y = Predicted), 
             color = "red", size = 3) +  # モデルの予測値を示す点を追加
  theme(legend.position = "right")+  # テーマの設定
  theme(
    plot.title = element_text(size = 20 * 2),    # タイトルのサイズを5倍に
    axis.title.x = element_text(size = 14 * 2),  # x軸ラベルのサイズを5倍に
    axis.title.y = element_text(size = 14 * 2),  # y軸ラベルのサイズを5倍に
    axis.text.x = element_text(size = 8 * 2),   # x軸目盛りのサイズを5倍に
    axis.text.y = element_text(size = 11 * 2),    # y軸目盛りのサイズを5倍に
    legend.position = "right",
    strip.text = element_text(size = 18),  # 83ms、117ms、150msの文字サイズを大きくする
    legend.text = element_text(size = 18),  # 凡例の文字サイズを大きくする
    legend.title = element_text(size = 18)  # 凡例タイトルの文字サイズを大きくする
    )


################################################################################

data_dprimeC <- data.frame(
  presentationtime = c("0_83ms","1_117ms", "2_150ms", rep("0_83ms",5),rep("1_117ms",5),rep("2_150ms",5)),
  imagetype =  c(rep("0_gray",3),"1_9deg","2_13deg","3_17deg","4_21deg","5_25deg","1_9deg","2_13deg","3_17deg","4_21deg","5_25deg","1_9deg","2_13deg","3_17deg","4_21deg","5_25deg"),
  dprime = rowMeans(dprime_mat),
  C = rowMeans(C_mat)
)

g <- ggplot(data_dprimeC, aes(x = dprime, y = C, group = presentationtime, fill = imagetype))
g <- g + geom_point(aes(shape = presentationtime, colour = imagetype), size = 4)
g <- g + labs(x = "dprime", y = "C")
g <- g + coord_cartesian(xlim = c(0,3.5))
g <- g + coord_cartesian(ylim = c(-0.6,0.8))
g <- g + theme_classic()
g <- g+ theme(
  axis.title.x = element_text(size = 14 * 2),  # x軸ラベルのサイズを5倍に
  axis.title.y = element_text(size = 14 * 2),  # y軸ラベルのサイズを5倍に
  axis.text.x = element_text(size = 7* 2),   # x軸目盛りのサイズを5倍に
  axis.text.y = element_text(size = 7 * 2),    # y軸目盛りのサイズを5倍に
  legend.position = "right",
  strip.text = element_text(size = 18),  # 83ms、117ms、150msの文字サイズを大きくする
  legend.text = element_text(size = 18),  # 凡例の文字サイズを大きくする
  legend.title = element_text(size = 18)  # 凡例タイトルの文字サイズを大きくする
)


plot(g)
	