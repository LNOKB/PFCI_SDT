#9 model -set λ as 1
library(tidyverse)
library(ggplot2)
library(effectsize)
library(plotly)
library(viridis)

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
sigma83ms <- 1
mu_gray <- 0


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
  
  # setting initial values
  guess <- c(
    initial_mu(data, 4),  # mu_9deg
    initial_mu(data, 5),  # mu_13deg
    initial_mu(data, 6),  # mu_17deg
    initial_mu(data, 7),  # mu_21deg
    initial_mu(data, 8),  # mu_25deg
    initial_mu(data, 19), # mu_color
    1, 1, 2         # sigma117ms, sigma150ms,  theta 
  )
  
  # fitting specifications
  lower_bounds <- c(0, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5)#0.001, 0.001, 
  upper_bounds <- c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 2.0, 2.0,  5.0)#)
  control_params <- list(
    "maxit" = 10000,
    "parscale" = c(1, 1, 1, 1, 1, 1, 0.001, 0.001,  1)
  )
  
  # fitting
  fit <- suppressWarnings(
    optim(
      PFCI_logL,
      par = guess,
      lower = lower_bounds,
      upper = upper_bounds,
      gr = NULL,
      method = "BFGS",
      control = control_params
    )
  )
  
  # outputs
  mu_9deg  <- fit$par[1]
  mu_13deg <- fit$par[2] 
  mu_17deg <- fit$par[3] 
  mu_21deg <- fit$par[4] 
  mu_25deg <- fit$par[5] 
  mu_color <- fit$par[6] 
  sigma117ms <-  fit$par[7] 
  sigma150ms <-  fit$par[8] 
  theta <-   fit$par[9] 
  logL <-         fit$value
  
  
  est <- data.frame(mu_9deg   = mu_9deg, 
                    mu_13deg  = mu_13deg, 
                    mu_17deg  = mu_17deg, 
                    mu_21deg  = mu_21deg, 
                    mu_25deg  = mu_25deg, 
                    mu_color = mu_color,
                    sigma117ms = sigma117ms,
                    sigma150ms = sigma150ms,
                    theta = theta,
                    logL = logL)
  return(list(est, global_predicted_data))
}


### Likelihood function
PFCI_logL <- function(x, inputs) {
  
  # target parameters
  mu_9deg  <- x[1]
  mu_13deg <- x[2] 
  mu_17deg <- x[3] 
  mu_21deg <- x[4] 
  mu_25deg <- x[5] 
  mu_color <- x[6] 
  sigma117ms <-  x[7] 
  sigma150ms <-  x[8] 
  theta <-   x[9] 
  
  # model predictions
  mu_chimeras <- c(mu_9deg, mu_13deg, mu_17deg,mu_21deg,mu_25deg)
  predicted_data <- matrix(NA, nrow = 21, ncol = 2)
  # gray image
  sigma_mat <- c(sigma83ms, sigma117ms, sigma150ms)
  for (cond in 1:3) {
    predicted_data[cond, 1] <- pnorm(theta, mean = mu_gray, sd = sigma_mat[cond]) # pred_gray_others_rate
    predicted_data[cond, 2] <- 1 - predicted_data[cond, 1] # pred_gray_color_rate
  }
  # chimera image
  mean_one_mat <- c(rep(mu_chimeras,3))
  sigma_mat <- c(rep(sigma83ms, 5), rep(sigma117ms, 5), rep(sigma150ms, 5))
  for (cond in 4:18) {
    predicted_data[cond, 1] <- pnorm(theta, mean = mean_one_mat[cond - 3], sd = sigma_mat[cond - 3]) # pred_chimera_others_rate
    predicted_data[cond, 2] <- 1 - predicted_data[cond, 1] # pred_chimera_color_rate
  }
  # full-color image
  sigma_mat <- c(sigma83ms, sigma117ms, sigma150ms)
  for (cond in 19:21) {
    predicted_data[cond, 1] <- pnorm(theta, mean = mu_color, sd = sigma_mat[cond - 18]) # pred_color_others_rate
    predicted_data[cond, 2] <- 1 -  predicted_data[cond, 1] # pred_color_color_rate
  }
  
  global_predicted_data <<- predicted_data 
  # log likelihood
  logL <- sum(data * log(predicted_data))
  if (is.nan(logL)) {
    logL <- -Inf
  }
  logL <- -logL
  return(logL)
  
}


### Conducting fitting on individual data
estimates <- c()
predicted_array <- array(NA, dim = c(21, 2, 44))
data_rate_array <- array(NA, dim = c(21, 2, 44))

for (i in 4:47) {
  
  subdata <- result[((i - 4)*21 + 1):((i - 4)*21 + 21), ]
  data <- matrix(NA, nrow = 21, ncol = 2)
  # gray image
  data[1:3, 1] <- unlist(subdata[1:3, 5] + subdata[1:3, 6])# _others
  data[1:3, 2] <- unlist(subdata[1:3, 7])# _color
  # chimera image
  data[4:18, 1] <- unlist(subdata[4:18, 8] + subdata[4:18, 9])# _others
  data[4:18, 2] <- unlist(subdata[4:18, 10])# _color
  # full-color image
  data[19:21, 1] <- unlist(subdata[19:21, 11] + subdata[19:21, 12])# _others
  data[19:21, 2] <- unlist(subdata[19:21, 13])# _color
  
  data_rate_array[, 1, (i - 3)] <- data[, 1]/(data[, 1] + data[, 2])
  data_rate_array[, 2, (i - 3)] <- data[, 2]/(data[, 1] + data[, 2])
  
  # fitting
  fit <- fit_PFCI_mle(data, add_constant = TRUE)
  df <- fit[[1]]
  df$sub <- i
  estimates <- rbind(estimates, df)
  predicted_array[, , (i - 3)] <-  fit[[2]]
}


### Preparation for data visualization
mean_predicted <- apply(predicted_array, c(1, 2), mean)*100
mean_data <-  apply(data_rate_array, c(1, 2), mean)*100
se_data <- (apply(data_rate_array, c(1, 2), sd)*100)/sqrt(44)
index_order <- c(
  1, 4:8,   19, 
  2, 9:13,  20, 
  3, 14:18, 21
)
mean_predicted2 <- mean_predicted[index_order, 2]
mean_data2 <- mean_data[index_order, 2]
se_data2 <- se_data[index_order, 2]


### t-test against sigma117ms, sigma150ms
for (i in 1:2) {
  now <- log(estimates[, i + 6])
  t_test_below1 <- t.test(now, mu = 0, alternative = "less")
  t_test_above1 <- t.test(now, mu = 0, alternative = "greater")
  effect <- cohens_d(now, mu = 0)
  print(t_test_below1)
  print(t_test_above1)
  print(effect)
}

# violin plot
data_parameter_plot <- data.frame(
  Parameters = rep(c("σ117ms", "σ150ms"), each = 44),
  Value = c(estimates[, 7], estimates[, 8])
)
parameters_graph <- ggplot(data_parameter_plot, aes(x = Parameters, y = Value)) +
  geom_violin(fill = "skyblue", color = "black", scale = "width") +  
  geom_jitter(width = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  labs(y = "Value") + 
  scale_y_continuous(breaks = seq(0.5, 1.5, length = 5), limits = c(0.5, 1.5)) +
  #scale_y_continuous(breaks = seq(0, 1.5, length = 4), limits = c(0, 1.5)) +
  scale_x_discrete("Parameters", labels = c(expression("σ"[117*ms]), expression("σ"[150*ms]))) +
  stat_summary(fun = mean, geom = "point", 
               shape = 18, size = 6, color = "black") +
  theme_classic() +  
  theme(
    plot.title =   element_text(size = 20 * 2),    
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14 * 2),  
    axis.text.x =  element_text(size = 13 * 2),  
    axis.text.y =  element_text(size = 11 * 2)    
  )
plot(parameters_graph)
ggsave(file = "parameters_graph1.png", plot = parameters_graph, dpi = 150, width = 8, height = 6)



### Estimated distribution
plot_sdt_distributions <- function(means, sds, attention_levels, image_types, colors) {
  data_SDT_plot <- data.frame()
  
  for (i in 1:length(attention_levels)) {
    for (j in 1:length(image_types)) {
      x <- seq(means[i, j] - 3 * sds[i, j], means[i, j] + 3 * sds[i, j], length.out = 100)
      y <- dnorm(x, mean = means[i, j], sd = sds[i, j]) 
      
      data_SDT_plot <- rbind(data_SDT_plot, data.frame(
        x = x,
        y = y,
        Attention = attention_levels[i],
        ImageType = image_types[j],
        color = colors[j]
      ))
    }
  }
  
  Distribution <- ggplot(data_SDT_plot, aes(x = x, y = y, color = ImageType)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = colors) +
    labs(x = "Strength of peripheral color signal",
         y = "Probability density") +
    scale_x_continuous(limits = c(-3, 6)) +
    scale_y_continuous(breaks = seq(0, 0.6, length = 4),limits = c(0, 0.6)) +
    geom_vline(xintercept = mean(estimates[, 9]), linetype = "dashed", color = "black") + 
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
  
  
  plot(Distribution)
  ggsave(file = "Distribution1.png", plot = Distribution, dpi = 150, width = 8, height = 6)
  
  
}

means_83ms <- c(
  mu_gray, mean(estimates[, 1]), mean(estimates[, 2]), mean(estimates[, 3]), 
  mean(estimates[, 4]), mean(estimates[, 5]),mean(estimates[, 6]))  
means_117ms <- means_83ms 
means_150ms <- means_83ms 
means <- rbind(means_83ms, means_117ms, means_150ms)
sd_83ms <-  c(rep(sigma83ms, 7))  
sd_117ms <- c(rep(mean(estimates[, 7]), 7))  
sd_150ms <- c(rep(mean(estimates[, 8]), 7))  
sds <- rbind(sd_83ms, sd_117ms, sd_150ms)

attention_levels <- factor(c("83 ms", "117 ms", "150 ms"), levels = c("83 ms", "117 ms", "150 ms"))
image_types <- factor(
  c("gray", "9 deg", "13 deg", "17 deg", "21 deg", "25 deg", "color"),
  levels = c("gray", "9 deg", "13 deg", "17 deg", "21 deg", "25 deg", "color"))
colors <- viridis(7, option = "plasma") 

plot_sdt_distributions(means, sds, attention_levels, image_types, colors)


### Model prediction
data_bar <- data.frame(
  Imagetype = factor(rep(c("gray","chimera 9 degree", "chimera 13 degree", "chimera 17 degree", 
                           "chimera 21 degree", "chimera 25 degree", "full-color"), 3),
                     levels = c("gray","chimera 9 degree", "chimera 13 degree", "chimera 17 degree", 
                                "chimera 21 degree", "chimera 25 degree", "full-color")),
  Proportion =  as.vector(t(mean_data2)),
  SE =  as.vector(t(se_data2)),
  Condition = factor(c(rep("83 ms", 7), rep("117 ms", 7), rep("150 ms", 7)), 
                     levels = c("83 ms", "117 ms", "150 ms"))
)

predicted_data_bar <- data.frame(
  Imagetype = factor(rep(c("gray","chimera 9 degree", "chimera 13 degree", "chimera 17 degree",
                           "chimera 21 degree", "chimera 25 degree", "full-color"), 3),
                     levels = c("gray","chimera 9 degree", "chimera 13 degree", "chimera 17 degree",
                                "chimera 21 degree", "chimera 25 degree", "full-color")),
  Predicted = as.vector(t(mean_predicted2)),  
  Condition = factor(c(rep("83 ms", 7),rep("117 ms", 7),rep("150 ms", 7)), 
                     levels = c("83 ms", "117 ms", "150 ms"))
)



rsquareds <- array(NA, dim = c(44, 1, 1))
for (ii in 1:44) {
  
  yhat <- predicted_array[,2,ii]
  y <- data_rate_array[,2,ii]
  ymean <- mean(y)
  rss <- sum((yhat - y)^2)
  tss <- sum((y - ymean)^2)
  
  rsquareds[ii] <- (1 - (rss/tss))
  meanR2 <- mean(rsquareds)
  sdR2 <- sd(rsquareds)
}

ann_text <- data.frame(Imagetype = "chimera 17 degree",Proportion = 100, Condition = factor("150 ms",levels = c("83 ms", "117 ms", "150 ms")))

bar_graph <- ggplot(data_bar, aes(x = Imagetype, y = Proportion)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = NULL, y = "Full-color response proportion (%)") +
  facet_grid(. ~ Condition) +  
  geom_point(data = predicted_data_bar, aes(x = Imagetype, y = Predicted), 
             color = "red", size = 3) +  
  geom_errorbar(data = data_bar, aes(ymin = Proportion - SE, ymax = Proportion + SE), 
                width = 0.2) + 
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
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),    
    limits = c(0, 100)                  
  ) + 
  geom_text(data = ann_text, label = paste("Summed logL =",  -round(sum(estimates[, 10]), 1)), 
            size = 10,
            color = "red")


plot(bar_graph)
ggsave(file = "bar_graph1.png", plot = bar_graph, dpi = 150, width = 14, height = 8)
