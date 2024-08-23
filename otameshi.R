library(tidyverse)#前処理に必要
library(magrittr) #パイプ演算子(%>%)を使うために必要

dat <- read.table("subdata.csv", header=T,sep=",") # 生データ

result = dat %>%             
  filter(colorcond == 1 | colorcond == 2) %>%  # 中心課題正答試行のみ抽出
  select(subnum, exposure, maskwidth, colorcond, TorF1, colorres, TorF2, mix_mix, mix_color, color_mix, color_color) %>%    # 分散分析に必要な変数の列だけを抽出して
  group_by(subnum, exposure, colorcond, maskwidth) %>%                # グループ化で並び替えて
  summarize(CR= sum(mix_mix), FA= sum(mix_color), miss= sum(color_mix), hit= sum(color_color),) %>%       # 要約統計（ここでは平均）を計算
  print()  

#write.csv(result, "otameshi.csv")

sub <- 1
dat1 <- read.table("hit.csv", header=F,sep=",") 
dat2 <- read.table("miss.csv", header=F,sep=",") 
dat3 <- read.table("correctrejection.csv", header=F,sep=",") 
dat4 <- read.table("fa.csv", header=F,sep=",") 
hit<-dat1[sub,2:16]
miss<-dat2[sub,2:16]
correctrejection<-dat3[sub+1,2:16]
fa<-dat4[sub+1,2:16]
kk

# add_constant = TRUE adds a small value to the response frequency vectors.
### Function for model fitting
fit_uvsdt_mle <- function(hit, miss, correctrejection, fa, add_constant = TRUE) {
  
  # correction against extreme estimates 
  if (add_constant) {
    hit <- hit + 0.001
    miss <- miss + 0.001
    correctrejection <- correctrejection + 0.001
    fa <- fa + 0.001
  }

  sigma83ms <- 1

  sigma117ms <- sigma150ms * a
  sigma150ms <- sigma150ms * b
  
  mu83ms_9deg <- 1
  
  mu83ms_13deg <- mu83ms_9deg * c
  mu83ms_17deg <- mu83ms_9deg * d   
  mu83ms_21deg <- mu83ms_9deg * e
  mu83ms_25deg <- mu83ms_9deg * f
  
  mu117ms_9deg <- mu83ms_9deg * g 
  mu117ms_13deg <- mu83ms_9deg * h
  mu117ms_17deg <- mu83ms_9deg * i 
  mu117ms_21deg <- mu83ms_9deg * j 
  mu117ms_25deg <- mu83ms_9deg * k
  
  mu150ms_9deg <- mu83ms_9deg * l
  mu150ms_13deg <- mu83ms_9deg * m
  mu150ms_17deg <- mu83ms_9deg * n
  mu150ms_21deg <- mu83ms_9deg * o
  mu150ms_25deg <- mu83ms_9deg * p
  
  # initial guess for parameter values
  

  a <- 1
  b <- 1
  
  c <- 1
  d <- 1
  e <- 1
  f <- 1
  g <- 1
  h <- 1
  i <- 1
  j <- 1
  k <- 1
  l <- 1
  m <- 1
  n <- 1
  o <- 1
  p <- 1
  
  cri <- 0 #β
  
  guess <- c(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, cri)
  
  # model fit
  fit <- suppressWarnings(optim(uvsdt_logL, 
                                par = guess, 
                                inputs = list("hit" = hit, "miss" = miss, "correctrejection" = correctrejection, "fa" = fa), 
                                gr = NULL, method = "BFGS", control = list("maxit" = 10000)))
  
  # outputs
  sigma117ms <- fit$par[1] * sigma83ms
  sigma150ms <- fit$par[2] * sigma83ms
  
  mu83ms_13deg <- fit$par[3] * mu83ms_9deg
  mu83ms_17deg <- fit$par[4] * mu83ms_9deg
  mu83ms_21deg <- fit$par[5] * mu83ms_9deg
  mu83ms_25deg <- fit$par[6] * mu83ms_9deg
  
  mu117ms_9deg <-  fit$par[7] * mu83ms_9deg
  mu117ms_13deg <- fit$par[8] * mu83ms_9deg
  mu117ms_17deg <- fit$par[9] * mu83ms_9deg
  mu117ms_21deg <- fit$par[10] * mu83ms_9deg
  mu117ms_25deg <- fit$par[11] * mu83ms_9deg
  
  mu150ms_9deg <-  fit$par[12] * mu83ms_9deg
  mu150ms_13deg <- fit$par[13] * mu83ms_9deg
  mu150ms_17deg <- fit$par[14] * mu83ms_9deg
  mu150ms_21deg <- fit$par[15] * mu83ms_9deg
  mu150ms_25deg <- fit$par[16] * mu83ms_9deg
  
    #da <-    mu / sqrt((1 + sigma^2) / 2)
  logL <- -fit$value
  
  #########?###############################################
  cri <- data.frame(matrix(vector(), 0, 2 * n_ratings - 1))
  for (i in 1:(2 * n_ratings - 1)) {
    cri[1, i] <- fit$par[2 * n_ratings + 2 - i]
  }
  #########################################################
  
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
                    cri = cri, 
                    logL = logL)
  return(est)
  
}

######################################################
### Likelihood function
uvsdt_logL <- function(x, inputs) {
  
  # target parameters
  mu <-    x[1]　
  sigma <- x[2]
  cri <-   x[3:(2 * inputs$n_ratings + 1)]　
  #【シミュレーション】今回は確信度をとっていないのでROC曲線はかけない
  
  # empirical data
  nr_s1 <- inputs$nr_s1
  nr_s2 <- inputs$nr_s2
  n_ratings <- inputs$n_ratings
  
  # model predictions
  #反応率　criterionより上側の面積を出す
  pred_far <- c(0, pnorm(0 - cri, 0, 1),            1)#s1分布
  pred_hr <-  c(0, pnorm((mu - cri) / sigma, 0, 1), 1)#S2分布なのでsigmaで割る
  #qnormは引数に確率を受け、確率密度を返す。例）累積確率が95%になる確率密度の値(1.6448)を返す
  #pnormは引数に確率密度を受て、確率を返す。例）確率密度がマイナス無限大から1.6448までの累積確率(94.9%)を返す
  #これだとなぜhit とfalse alarmで良い？
  
  #反応数
  pred_nr_s1 <- sum(nr_s1) * diff(pred_far)
  pred_nr_s2 <- sum(nr_s2) * diff(pred_hr)
  #diffは連続する要素間の差を計算。例 x <- c(5,3,4) diff(x) = (-2  1)
  
  # log likelihood
  logL <- sum(nr_s1 * log(pred_nr_s1 / sum(nr_s1)) + nr_s2 * log(pred_nr_s2 / sum(nr_s2)))#対数にしたので、S1刺激に対する対数尤度とS2刺激に対する対数尤度の単純な和をとればＯＫ
  if (is.nan(logL)) {
    logL <- -Inf
  }
  logL <- -logL　#対数尤度を対数損失にする
  return(logL)
  
}


#尤度により、1000回中の700回なのか、100回中の70回なのかを区別できる


#非等分散の場合は、small cとしてしか出ない。
#等分散なら、S1（左側分布）とS2（右側分布）の交点は、正答率を最大化する点という意味がある。それより左／右だったら、「右選択しやすいバイアス」などとして有意味に解釈できる。
#非等分散の場合は、densityが違うので、何の意味も持たない。交点も複数存在する。
#なので、非等分散の場合、S1のピークを0とし、それを基準とした計量をする。


### Fitting
fit <- fit_uvsdt_mle(nr_s1, nr_s2, add_constant = FALSE)
fit

###################################
#mean+variance,mean,variance
###################################
if (par[2] < par[3] & par[3] < par[4]) {
  return(mse)
} else {
  return(999999)
}