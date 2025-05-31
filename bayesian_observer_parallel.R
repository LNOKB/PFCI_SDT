#Bayesian observer

library(tidyverse)
library(ggplot2)
library(effectsize)
library(plotly)
library(viridis)

#Parallel 
library(doParallel)
library(foreach)
num_cores <- parallel::detectCores(logical = FALSE)
cl <- makeCluster(num_cores - 1) #keep 1 core for system stability
registerDoParallel(cl)


dat <- read_csv("subdata.csv")
result = dat %>%             
  filter(TorF1 %in% c(0, 1)) %>%
  select(
    subnum, exposure, maskwidth, colorcond, TorF1, colorres, TorF2, 
    gray_gray, gray_mix, gray_color, mix_gray, mix_mix, mix_color, 
    color_gray, color_mix, color_color
  ) %>%    
  group_by(subnum, colorcond, exposure, maskwidth) %>%                
  summarize(
    gray_gray = sum(gray_gray), gray_mix = sum(gray_mix), gray_color = sum(gray_color), 
    mix_gray = sum(mix_gray), mix_mix = sum(mix_mix), mix_color = sum(mix_color), 
    color_gray = sum(color_gray), color_mix = sum(color_mix), color_color = sum(color_color))     

# fixed parameters
sigma <- 1
mu83ms_gray <- 0
prior_gray <- 2/12
x_seq <- seq(-2, 6, length.out = 4000)

col_names <- c("mu83ms_9deg", "mu83ms_13deg", "mu83ms_17deg", "mu83ms_21deg", "mu83ms_25deg",
               "mu83ms_color", "lambda117ms", "lambda150ms", "theta83ms", "theta117ms", "theta150ms", 
               "prior_chimera_levels","prior_fullcolor", "logL")


posterior_probs <- function(x, mu83ms_chimeras, mu83ms_color, lambdas, prior_fullcolor) {
  criteria <- numeric(3)  
  prior_chimera_levels <- (1 - prior_gray - prior_fullcolor)/5
  
  for (cond2 in 1:3) {
    
    # Posterior probability calculation
    numerator_chimera_levels <- sapply(mu83ms_chimeras, function(mu) {
      dnorm(x, mean = mu * lambdas[cond2], sd = sigma) * prior_chimera_levels
    })
    numerator_chimera <- rowSums(numerator_chimera_levels)
    numerator_gray <- dnorm(x, mean = mu83ms_gray * lambdas[cond2], sd = sigma) * prior_gray
    numerator_full <- dnorm(x, mean = mu83ms_color * lambdas[cond2], sd = sigma) * prior_fullcolor
    marginal_likelihood <- numerator_chimera + numerator_gray + numerator_full
    posterior_gray <- numerator_gray / marginal_likelihood
    posterior_chimera <- numerator_chimera / marginal_likelihood
    posterior_full <- numerator_full / marginal_likelihood
    
    # Finding the minimum x in which the largest posterior probability is full-color
    posterior_mat <- cbind(posterior_gray, posterior_chimera, posterior_full)
    max_class <- c("Gray", "Chimera", "FullColor")[apply(posterior_mat, 1, which.max)]
    tmp <- data.frame(
      x = x,
      cond = cond2,
      numerator_chimera = numerator_chimera,
      numerator_gray = numerator_gray,
      numerator_full = numerator_full,
      posterior_chimera = posterior_chimera,
      posterior_gray = posterior_gray,
      posterior_full = posterior_full,
      max_class = max_class
    )
    color_tmp <- subset(tmp, max_class == "FullColor")
    min_x_row <- color_tmp[which.min(color_tmp$x), ]
    # result_ideal <- rbind(result_ideal, tmp)
    criteria[cond2] <- min_x_row$x
  }
  return(criteria)
}


### Fitting function
fit_PFCI_mle <- function(data, add_constant = TRUE) {
  if (add_constant) {
    data <- data + 0.5
  }
  
  initial_mu <- function(data, each_index) {
    each_color_rate <- data[each_index, 2]/(data[each_index, 1] + data[each_index, 2])
    gray_color_rate <- data[1, 2]/(data[1, 1] + data[1, 2])
    return(qnorm(each_color_rate) - qnorm(gray_color_rate))
  }
  
  color_83ms <- sum(data[c(1, 4:8, 19), 2])
  others_83ms <- sum(data[c(1, 4:8, 19), 1])
  color_rate_83ms <- color_83ms / (color_83ms + others_83ms)
  
  color_117ms <- sum(data[c(2, 9:13, 20), 2])
  others_117ms <- sum(data[c(2, 9:13, 20), 1])
  color_rate_117ms <- color_117ms / (color_117ms + others_117ms)
  
  color_150ms <- sum(data[c(3, 14:18, 21), 2])
  others_150ms <- sum(data[c(3, 14:18, 21), 1])
  color_rate_150ms <- color_150ms / (color_150ms + others_150ms)
  
  initial_lambda117ms <- color_rate_117ms/color_rate_83ms
  initial_lambda150ms <- color_rate_150ms/color_rate_83ms
  
  # setting initial values
  guess <- c(
    initial_mu(data, 4),  # mu83ms_9deg
    initial_mu(data, 5),  # mu83ms_13deg
    initial_mu(data, 6),  # mu83ms_17deg
    initial_mu(data, 7),  # mu83ms_21deg
    initial_mu(data, 8),  # mu83ms_25deg
    initial_mu(data, 19), # mu83ms_color
    initial_lambda117ms,
    initial_lambda150ms ,
    5/12
  )
  #   1.2,1.4
  # )

  lambda_parscale <- 1
  prior_parscale <- 0.001
  
  if (i %in% c(11, 13, 19, 20, 22, 26, 31, 36, 37, 39, 44)) {
    lambda_parscale <- 0.001
  } else {
  }

  if (i %in% c(11, 100)) {
    prior_parscale <- 1
  } else {
  }
  
  #origin: 11 19 20 22 31 36 39 add:13 26 37 44
  
  # fitting specifications
  lower_bounds <- c(0, 0, 0, 0, 0, 0, 0.5, 0.5, 0)
  upper_bounds <- c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 2.0, 2.0, 10/12)
  control_params <- list(
    "maxit" = 10000,
    "parscale" = c(1, 1, 1, 1, 1, 1, lambda_parscale, lambda_parscale, prior_parscale)
  )
  

  PFCI_logL_wrapper <- function(par) {
    PFCI_logL(par)$logL  
  }
  
  
  fit <- suppressWarnings(
    optim(
      PFCI_logL_wrapper,
      par = guess,
      lower = lower_bounds,
      upper = upper_bounds,
      gr = NULL,
      method = "BFGS",
      control = control_params
    )
  )
  
  # outputs
  fitresult <- PFCI_logL(fit$par)  
  criteria <- fitresult$criteria
  predicted_data <- fitresult$predicted_data
  est <- as.data.frame(matrix(0, nrow = 1, ncol = length(col_names)))
  colnames(est) <- col_names
  est[, col_names] <- c(
    fit$par[1:8],
    criteria[1:3],
    (1 - prior_gray - fit$par[9])/5,
    fit$par[9],
    fit$value
  )
  
  return(list(est = est, predicted_data = predicted_data))
}


### Likelihood function
#PFCI_logL <- function(x, inputs) {
PFCI_logL <- function(x) {
    
  predicted_data <- matrix(NA, nrow = 21, ncol = 2)
  
  # target parameters
  mu83ms_9deg  <- x[1]
  mu83ms_13deg <- x[2] 
  mu83ms_17deg <- x[3] 
  mu83ms_21deg <- x[4] 
  mu83ms_25deg <- x[5] 
  mu83ms_color <- x[6] 
  lambda117ms <-  x[7] 
  lambda150ms <-  x[8] 
  prior_fullcolor <- x[9]
  
  mu83ms_chimeras <- c(mu83ms_9deg, mu83ms_13deg, mu83ms_17deg, mu83ms_21deg, mu83ms_25deg) 
  lambdas <- c(1, lambda117ms, lambda150ms) 
  criteria <- posterior_probs(x_seq, mu83ms_chimeras, mu83ms_color, lambdas, prior_fullcolor)
  
  #Prediction of full-color / other response 
  mean_gray_mat <- mu83ms_gray * lambdas
  for (cond in 1:3) {
    predicted_data[cond, 1] <- pnorm(criteria[cond], mean = mean_gray_mat[cond], sd = sigma) # pred_gray_others_rate
    predicted_data[cond, 2] <- 1 - predicted_data[cond, 1] # pred_gray_color_rate
  }
  
  mean_chimera_mat <- rep(mu83ms_chimeras, times = length(lambdas)) *
    rep(lambdas, each = length(mu83ms_chimeras))
  criteria_chimera <- rep(criteria, each = 5)
  for (cond in 4:18) {
    predicted_data[cond, 1] <- pnorm(criteria_chimera[cond-3], mean = mean_chimera_mat[cond-3], sd = sigma) # pred_chimera_others_rate
    predicted_data[cond, 2] <- 1 - predicted_data[cond, 1] # pred_chimera_color_rate
  }
  
  mean_color_mat <- mu83ms_color * lambdas
  for (cond in 19:21) {
    predicted_data[cond, 1] <- pnorm(criteria[cond-18], mean = mean_color_mat[cond-18], sd = sigma) # pred_color_others_rate
    predicted_data[cond, 2] <- 1 - predicted_data[cond, 1] # pred_color_color_rate
  }
  
  #Calculation of log-Likelihood
  logL <- sum(data * log(predicted_data))
  if (is.nan(logL)) {
    logL <- -Inf
  }
  logL <- -logL
  
  return(list(logL = logL, criteria = criteria, predicted_data = predicted_data))
}


### Conducting fitting on individual data
estimates <- data.frame(matrix(NA, nrow = 44, ncol = length(col_names)))
colnames(estimates) <- col_names

data_rate_array <- array(NA, dim = c(21, 2, 44))
predicted_array <- array(NA, dim = c(21, 2, 44))

clusterExport(cl, varlist = c(
  "fit_PFCI_mle", "result"
))

# i_vals <- setdiff(4:47, c(11, 19, 20, 22, 31, 36, 39))
# i_vals <- c(11, 19, 20, 22, 31, 36, 39)
fromsub <- 4
tosub <- 10
results <- foreach(i = fromsub:tosub, .combine = 'rbind') %dopar% {
  
# results <- foreach(i = c(5, 7, 10, 11, 14, 19, 21, 22, 26, 30, 33, 36, 41, 43, 44), .combine = 'rbind') %dopar% {
  subdata <- result[((i - 4) * 21 + 1):((i - 4) * 21 + 21), ]
  data <- matrix(NA, nrow = 21, ncol = 2)
  
  # gray image
  data[1:3, 1] <- unlist(subdata[1:3, 5] + subdata[1:3, 6])  # _others
  data[1:3, 2] <- unlist(subdata[1:3, 7])                    # _color
  
  # chimera image
  data[4:18, 1] <- unlist(subdata[4:18, 8] + subdata[4:18, 9])  # _others
  data[4:18, 2] <- unlist(subdata[4:18, 10])                    # _color
  
  # full-color image
  data[19:21, 1] <- unlist(subdata[19:21, 11] + subdata[19:21, 12])  # _others
  data[19:21, 2] <- unlist(subdata[19:21, 13])                       # _color
  
  data_rate <- matrix(NA, nrow = 21, ncol = 2)
  data_rate[, 1] <- data[, 1] / (data[, 1] + data[, 2])
  data_rate[, 2] <- data[, 2] / (data[, 1] + data[, 2])
  
  fit <- fit_PFCI_mle(data, add_constant = TRUE)
  # df <- fit$est
  # df$sub <- i
  # 
  # list(est = df, pred = fit$predicted_data, rate = data_rate)
  est_vec <- as.numeric(fit$est)  # 13個のパラメータをベクトルに
  pred_vec <- as.numeric(t(fit$predicted_data))  # 21*2 = 42個、列ベクトル化
  rate_vec <- as.numeric(t(data_rate))            # 同様に42個
  
  # 名前付け用のベクトル作成
  names_est <- names(fit$est) #paste0("est_", names(fit$est))
  names_pred <- paste0("pred_", rep(1:21, each=2), "_", rep(1:2, times=21))
  names_rate <- paste0("rate_", rep(1:21, each=2), "_", rep(1:2, times=21))
  
  df <- data.frame(t(c(est_vec, pred_vec, rate_vec)))
  colnames(df) <- c(names_est, names_pred, names_rate)
  df$sub <- i
  df
}

#8
for (i in (fromsub - 3):(tosub - 3)){
  estvec <- unlist(unname(results[i, 1:14])) #1:12
  estimates[i, ] <- estvec
  
  predvec <- unlist(unname(results[i, 15:56])) #13:54
  pred_data <- matrix(predvec, nrow = 21, ncol = 2, byrow = TRUE)
  predicted_array[, , i] <- pred_data
  
  datavec <- unlist(unname(results[i, 57:98])) #55:96
  data_rate <- matrix(datavec, nrow = 21, ncol = 2, byrow = TRUE)
  data_rate_array[, , i] <- data_rate
}

# クラスターの停止
stopCluster(cl)


### Each subject plot
for (i in fromsub : tosub){
  #47
  #Bar plot
  data_rate <- data_rate_array[, , (i - 3)]
  predicted_data <- predicted_array[, , (i - 3)]
  
  index_order <- c(
    1, 4:8,   19,
    2, 9:13,  20,
    3, 14:18, 21
  )
  predicted2 <- predicted_data[index_order, 2] * 100
  data_rate2 <- data_rate[index_order, 2] * 100
  
  data_bar_sub <- data.frame(
    Imagetype = factor(rep(c("gray","chimera 9 degree", "chimera 13 degree", "chimera 17 degree",
                             "chimera 21 degree", "chimera 25 degree", "full-color"), 3),
                       levels = c("gray","chimera 9 degree", "chimera 13 degree", "chimera 17 degree",
                                  "chimera 21 degree", "chimera 25 degree", "full-color")),
    Proportion =  as.vector(t(data_rate2)),
    Condition = factor(c(rep("83 ms", 7), rep("117 ms", 7), rep("150 ms", 7)),
                       levels = c("83 ms", "117 ms", "150 ms"))
  )
  
  predicted_data_bar_sub <- data_bar_sub
  predicted_data_bar_sub <- predicted_data_bar_sub %>% rename(Predicted = Proportion)
  predicted_data_bar_sub$Predicted <- as.vector(t(predicted2))
  ann_text <- data.frame(Imagetype = "chimera 17 degree",Proportion = 100, Condition = factor("150 ms",levels = c("83 ms", "117 ms", "150 ms")))
  
  bar_graph_sub <- ggplot(data_bar_sub, aes(x = Imagetype, y = Proportion)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = NULL, y = "Full-color response proportion (%)") +
    facet_grid(. ~ Condition) +
    geom_point(data = predicted_data_bar_sub, aes(x = Imagetype, y = Predicted),
               color = "red", size = 3) +
    theme(legend.position = "right", text = element_text(family = "Arial")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title =   element_text(size = 20 * 2),
      axis.title.x = element_text(size = 14 * 2),
      axis.title.y = element_text(size = 14 * 2),
      axis.text.x =  element_text(size = 24, angle = 45, hjust = 1, color = "black"),
      axis.text.y =  element_text(size = 11 * 2),
      legend.position = "right",
      strip.text =   element_text(size = 36),
      legend.text =  element_text(size = 18),
      legend.title = element_text(size = 18)  #
    ) +
    geom_text(data = ann_text, label = paste("Log-likelihood =",  -round(estimates[(i - 3), 14], 1)),
              size = 7,
              color = "red") #12
  
  ggsave(file = file.path("SubPlotBar2", paste0("bar_graph_bayes_", i, ".png")),
         plot = bar_graph_sub, dpi = 150, width = 14, height = 8)
  
  
  #Distribution
  
  plot_sdt_distributions_ideal_sub <- function(means, sds, attention_levels, image_types, colors) {
    data_SDT_plot_ideal_sub <- data.frame()
    
    for (ii in 1:length(attention_levels)) {
      for (jj in 1:length(image_types)) {
        x <- seq(means[ii, jj] - 3 * sds[ii, jj], means[ii, jj] + 3 * sds[ii, jj], length.out = 100)
        y <- dnorm(x, mean = means[ii, jj], sd = sds[ii, jj])
        
        data_SDT_plot_ideal_sub <- rbind(data_SDT_plot_ideal_sub, data.frame(
          x = x,
          y = y,
          Attention = attention_levels[ii],
          ImageType = image_types[jj],
          color = colors[jj]
        ))
      }
    }
    
    vlines_ideal_sub <- data.frame(
      Attention = attention_levels,
      xintercept = c(estimates[(i - 3), 9], estimates[(i - 3), 10], estimates[(i - 3), 11]))
    
    Distribution_ideal_sub <- ggplot(data_SDT_plot_ideal_sub, aes(x = x, y = y, color = ImageType)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = colors) +
      labs(x = "Strength of peripheral color signal",
           y = "Probability density") +
      scale_x_continuous(limits = c(-3, 8)) +
      scale_y_continuous(breaks = seq(0, 0.6, length = 4),limits = c(0, 0.6)) +
      geom_vline(data = vlines_ideal_sub, aes(xintercept = xintercept),
                 linetype = "dashed", color = "black") +
      facet_wrap(~ Attention, nrow = 3, scales = "free_y") +
      theme_minimal(base_size = 18) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.3, "cm"),
        legend.position =  c(0.1, 5),
      )
    
    ggsave(file = file.path("SubPlotDistribution2", paste0("Distribution_ideal_", i, ".png")),
           plot = Distribution_ideal_sub, dpi = 150, width = 8, height = 6)
  }
  
  means_83ms <- c(mu83ms_gray, estimates[(i - 3), 1], estimates[(i - 3), 2],
                  estimates[(i - 3), 3], estimates[(i - 3), 4],
                  estimates[(i - 3), 5], estimates[(i - 3), 6])
  means_117ms <- means_83ms * estimates[(i - 3), 7]
  means_150ms <- means_83ms * estimates[(i - 3), 8]
  means <- rbind(means_83ms, means_117ms, means_150ms)
  sd_83ms <-  c(rep(sigma, 7))
  sd_117ms <- sd_83ms
  sd_150ms <- sd_83ms
  sds <- rbind(sd_83ms, sd_117ms, sd_150ms)
  
  attention_levels <- factor(c("83 ms", "117 ms", "150 ms"), levels = c("83 ms", "117 ms", "150 ms"))
  image_types <- factor(
    c("gray", "9 deg", "13 deg", "17 deg", "21 deg", "25 deg", "color"),
    levels = c("gray", "9 deg", "13 deg", "17 deg", "21 deg", "25 deg", "color"))
  colors <- viridis(7, option = "plasma")
  
  plot_sdt_distributions_ideal_sub(means, sds, attention_levels, image_types, colors)
}

### End of Each subject plot
# 
# ### Preparation for data visualization
# mean_predicted <- apply(predicted_array, c(1, 2), mean)*100
# mean_data <-  apply(data_rate_array, c(1, 2), mean)*100
# se_data <- (apply(data_rate_array, c(1, 2), sd)*100)/sqrt(44)
# index_order <- c(
#   1, 4:8,   19,
#   2, 9:13,  20,
#   3, 14:18, 21
# )
# mean_predicted2 <- mean_predicted[index_order, 2]
# mean_data2 <- mean_data[index_order, 2]
# se_data2 <- se_data[index_order, 2]
# 
# 
# ## t-test
# print(t.test(estimates[, 9], estimates[, 10], paired = TRUE))
# print(t.test(estimates[, 9], estimates[, 11], paired = TRUE))
# print(t.test(estimates[, 10], estimates[, 11], paired = TRUE))
# 
# print(cohens_d(estimates[, 9] - estimates[, 10], mu = 0))
# print(cohens_d(estimates[, 9] - estimates[, 11], mu = 0))
# print(cohens_d(estimates[, 10] - estimates[, 11], mu = 0))
# 
# # violin plot
# data_parameter_plot <- data.frame(
#   Parameters = rep(c("1.θ83ms","2.θ117ms", "3.θ150ms"), each = 44),
#   Value = c(estimates[, 9], estimates[, 10], estimates[, 11])
# )
# parameters_graph <- ggplot(data_parameter_plot, aes(x = Parameters, y = Value)) +
#   geom_violin(fill = "skyblue", color = "black") +
#   geom_jitter(width = 0.1) +
#   labs(y = "Value") +
#   scale_y_continuous(breaks = seq(1, 4, length = 4), limits = c(1, 4)) +
#   scale_x_discrete("Parameters",
#                    labels = c(expression("θ"[Bayes_83*ms]), expression("θ"[Bayes_117*ms]),
#                               expression("θ"[Bayes_150*ms]))) +
#   stat_summary(fun = mean, geom = "point",
#                shape = 18, size = 6, color = "black") +
#   theme_classic() +
#   theme(
#     plot.title =   element_text(size = 20 * 2),
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(size = 14 * 2),
#     axis.text.x =  element_text(size = 13 * 2),
#     axis.text.y =  element_text(size = 11 * 2)
#   )
# plot(parameters_graph)
# ggsave(file = "parameters_graph_bayes.png", plot = parameters_graph, dpi = 150, width = 8, height = 6)
# #
# # violin plot
# data_parameter_plot_lambda <- data.frame(
#   Parameters = rep(c("1.λ117ms", "2.λ150ms"), each = 44),
#   Value = c(estimates[, 7], estimates[, 8])
# )
# parameters_graph_lambda <- ggplot(data_parameter_plot_lambda, aes(x = Parameters, y = Value)) +
#   geom_violin(fill = "skyblue", color = "black") +
#   geom_jitter(width = 0.1) +
#   labs(y = "Value") +
#   scale_y_continuous(breaks = seq(0.5, 2, length = 4), limits = c(0.5, 2)) +
#   scale_x_discrete("Parameters",
#                    labels = c(expression("λ"[Bayes_117*ms]),expression("λ"[Bayes_150*ms]))) +
#   stat_summary(fun = mean, geom = "point",
#                shape = 18, size = 6, color = "black") +
#   theme_classic() +
#   theme(
#     plot.title =   element_text(size = 20 * 2),
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(size = 14 * 2),
#     axis.text.x =  element_text(size = 13 * 2),
#     axis.text.y =  element_text(size = 11 * 2)
#   )
# plot(parameters_graph)
# ggsave(file = "parameters_graph_bayes_lambda.png", plot = parameters_graph_lambda, dpi = 150, width = 8, height = 6)
# #
# #
# 
# ### Estimated distribution
# plot_sdt_distributions <- function(means, sds, attention_levels, image_types, colors) {
#   data_SDT_plot <- data.frame()
# 
#   for (i in 1:length(attention_levels)) {
#     for (j in 1:length(image_types)) {
#       x <- seq(means[i, j] - 3 * sds[i, j], means[i, j] + 3 * sds[i, j], length.out = 100)
#       y <- dnorm(x, mean = means[i, j], sd = sds[i, j])
# 
#       data_SDT_plot <- rbind(data_SDT_plot, data.frame(
#         x = x,
#         y = y,
#         Attention = attention_levels[i],
#         ImageType = image_types[j],
#         color = colors[j]
#       ))
#     }
#   }
# 
#   vlines <- data.frame(
#     Attention = attention_levels,
#     xintercept = c(mean(estimates[, 9]), mean(estimates[, 10]), mean(estimates[, 11])))
# 
#   Distribution <- ggplot(data_SDT_plot, aes(x = x, y = y, color = ImageType)) +
#     geom_line(size = 1.2) +
#     scale_color_manual(values = colors) +
#     labs(x = "Strength of peripheral color signal",
#          y = "Probability density") +
#     scale_x_continuous(limits = c(-3, 6)) +
#     scale_y_continuous(breaks = seq(0, 0.6, length = 4),limits = c(0, 0.6)) +
#     geom_vline(data = vlines, aes(xintercept = xintercept),
#                linetype = "dashed", color = "black") +
#     facet_wrap(~ Attention, nrow = 3, scales = "free_y") +
#     theme_minimal(base_size = 18) +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       axis.line = element_line(size = 0.5, color = "black"),
#       axis.ticks = element_line(color = "black"),
#       axis.ticks.length = unit(0.3, "cm"),
#       legend.position =  c(0.1, 5),
#     )
# 
#   plot(Distribution)
#   ggsave(file = "Distribution_bayes.png", plot = Distribution, dpi = 150, width = 8, height = 6)
# 
# 
# }
# 
# means_83ms <- c(
#   mu83ms_gray, mean(estimates[, 1]), mean(estimates[, 2]), mean(estimates[, 3]),
#   mean(estimates[, 4]), mean(estimates[, 5]),mean(estimates[, 6]))
# means_117ms <- means_83ms * mean(estimates[, 7])
# means_150ms <- means_83ms * mean(estimates[, 8])
# means <- rbind(means_83ms, means_117ms, means_150ms)
# sd_83ms <-  c(rep(sigma, 7))
# sd_117ms <- sd_83ms
# sd_150ms <- sd_83ms
# sds <- rbind(sd_83ms, sd_117ms, sd_150ms)
# 
# attention_levels <- factor(c("83 ms", "117 ms", "150 ms"), levels = c("83 ms", "117 ms", "150 ms"))
# image_types <- factor(
#   c("gray", "9 deg", "13 deg", "17 deg", "21 deg", "25 deg", "color"),
#   levels = c("gray", "9 deg", "13 deg", "17 deg", "21 deg", "25 deg", "color"))
# colors <- viridis(7, option = "plasma")
# 
# plot_sdt_distributions(means, sds, attention_levels, image_types, colors)
# 
# ### Model prediction
# data_bar <- data.frame(
#   Imagetype = factor(rep(c("gray","chimera 9 degree", "chimera 13 degree", "chimera 17 degree",
#                            "chimera 21 degree", "chimera 25 degree", "full-color"), 3),
#                      levels = c("gray","chimera 9 degree", "chimera 13 degree", "chimera 17 degree",
#                                 "chimera 21 degree", "chimera 25 degree", "full-color")),
#   Proportion =  as.vector(t(mean_data2)),
#   SE =  as.vector(t(se_data2)),
#   Condition = factor(c(rep("83 ms", 7), rep("117 ms", 7), rep("150 ms", 7)),
#                      levels = c("83 ms", "117 ms", "150 ms"))
# )
# 
# predicted_data_bar <- data_bar
# predicted_data_bar <- predicted_data_bar %>% rename(Predicted = Proportion)
# predicted_data_bar$Predicted <- as.vector(t(mean_predicted2))
# 
# # rsquareds <- array(NA, dim = c(44, 1, 1))
# # for (ii in 1:44) {
# #
# #   yhat <- predicted_array_ideal[,2,ii]
# #   y <- data_rate_array[,2,ii]
# #   ymean <- mean(y)
# #   rss <- sum((yhat - y)^2)
# #   tss <- sum((y - ymean)^2)
# #
# #   rsquareds[ii] <- (1 - (rss/tss))
# #   meanR2 <- mean(rsquareds)
# #   sdR2 <- sd(rsquareds)
# # }
# 
# ann_text <- data.frame(Imagetype = "chimera 17 degree",Proportion = 100,
#                        Condition = factor("150 ms",levels = c("83 ms", "117 ms", "150 ms")))
# 
# bar_graph <- ggplot(data_bar, aes(x = Imagetype, y = Proportion)) +
#   geom_bar(stat = "identity", position = "stack") +
#   labs(x = NULL, y = "Full-color response proportion (%)") +
#   facet_grid(. ~ Condition) +
#   geom_point(data = predicted_data_bar, aes(x = Imagetype, y = Predicted),
#              color = "red", size = 3) +
#   geom_errorbar(data = data_bar, aes(ymin = Proportion - SE, ymax = Proportion + SE),
#                 width = 0.2) +
#   theme(legend.position = "right", text = element_text(family = "Arial")) +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     plot.title =   element_text(size = 20 * 2),
#     axis.title.x = element_text(size = 14 * 2),
#     axis.title.y = element_text(size = 14 * 2),
#     axis.text.x =  element_text(size = 24, angle = 45, hjust = 1, color = "black"),
#     axis.text.y =  element_text(size = 11 * 2),
#     legend.position = "right",
#     strip.text =   element_text(size = 36),
#     legend.text =  element_text(size = 18),
#     legend.title = element_text(size = 18)  #
#   ) +
#   geom_text(data = ann_text, label = paste("Summed log-likelihood =",  -round(sum(estimates[, 12]), 1)),
#           size = 7,
#           color = "red")
# # +
# #   scale_y_continuous(
# #   breaks = seq(0, 100, by = 20),
# #   limits = c(0, 100)
# # )
# 
# plot(bar_graph)
# ggsave(file = "bar_graph_bayes.png", plot = bar_graph, dpi = 150, width = 14, height = 8)
