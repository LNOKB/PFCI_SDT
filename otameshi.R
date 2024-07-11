sub <- 1
dat1 <- read.table("hit.csv", header=F,sep=",") 
dat2 <- read.table("miss.csv", header=F,sep=",") 
dat3 <- read.table("correctrejection.csv", header=F,sep=",") 
dat4 <- read.table("fa.csv", header=F,sep=",") 
hit<-dat1[sub+1,2:16]
miss<-dat2[sub+1,2:16]
correctrejection<-dat3[sub+1,2:16]
fa<-dat4[sub+1,2:16]
# \response frequency vectors 


# add_constant = TRUE adds a small value to the response frequency vectors.
### Function for model fitting
fit_uvsdt_mle <- function(hit, miss, correctrejection, fa, add_constant = FALSE) {
  
  # correction against extreme estimates 
  if (add_constant) {
    hit <- hit + (1 / length(hit))
    miss <- miss + (1 / length(miss))
    correctrejection <- miss + (1 / length(correctrejection))
    fa <- fa + (1 / length(fa))
  }

  
  # initial guess for parameter values
  mu <- 0
  sigma <- 1.5
  cri <-  0
  cri <-   0
  guess <- c(mu, sigma, cri)
  
  # model fit
  fit <- suppressWarnings(optim(uvsdt_logL, 
                                par = guess, 
                                inputs = list("hit" = hit, "miss" = miss, "correctrejection" = correctrejection, "fa" = fa), 
                                gr = NULL, method = "BFGS", control = list("maxit" = 10000)))
  
  # outputs
  mu <-    fit$par[1]
  sigma <- fit$par[2]
  da <-    mu / sqrt((1 + sigma^2) / 2)
  logL <- -fit$value
  #########?############
  cri <- data.frame(matrix(vector(), 0, 2 * n_ratings - 1))
  for (i in 1:(2 * n_ratings - 1)) {
    cri[1, i] <- fit$par[2 * n_ratings + 2 - i]
  }
  est <- data.frame(mu = mu, sigma = sigma, da = da, cri = cri, logL = logL)
  return(est)
  
}

######################################################
### Likelihood function
uvsdt_logL <- function(x, inputs) {
  
  # target parameters
  mu <-    x[1]　#【シミュレーション】これが刺激に応じた値になるように
  sigma <- x[2]
  cri <-   x[3:(2 * inputs$n_ratings + 1)]　#【シミュレーション】ここをいくつかの値（5水準とか？）に変えて推定の際のパラメータのあたりをつける
  #【シミュレーション】今回は分散が条件間で違うようにする。たとえば、最も短い呈示時間のときのS1の平均μが０、標準偏差σが1とおく。そのほかの条件のμやσが刺激の関数で決まるようにする。
  #【シミュレーション】今回は確信度をとっていないのでROC曲線はかけない
  
  # empirical data
  nr_s1 <- inputs$nr_s1
  nr_s2 <- inputs$nr_s2
  n_ratings <- inputs$n_ratings
  
  # model predictions
  #反応率　criterionより上側の面積を出す
  pred_far <- c(0, pnorm(0 - cri, 0, 1),            1)
  pred_hr <-  c(0, pnorm((mu - cri) / sigma, 0, 1), 1)#S2分布なのでsigmaで割る
  #反応数
  pred_nr_s1 <- sum(nr_s1) * diff(pred_far)
  pred_nr_s2 <- sum(nr_s2) * diff(pred_hr)
  
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
