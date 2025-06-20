library(tidyverse)
library(ggplot2)
library(effectsize)
library(plotly)
library(viridis)

# ---- Fitting based on mean+variance model ----
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
mu83ms_gray <- 0

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
    initial_mu(data, 4),  # mu83ms_9deg
    initial_mu(data, 5),  # mu83ms_13deg
    initial_mu(data, 6),  # mu83ms_17deg
    initial_mu(data, 7),  # mu83ms_21deg
    initial_mu(data, 8),  # mu83ms_25deg
    initial_mu(data, 19), # mu83ms_color
    1, 1, 1, 1, 2         # lambda117ms, lambda150ms, sigma117ms, sigma150ms, theta 
  )
  
  # fitting specifications
  lower_bounds <- c(0, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5)
  upper_bounds <- c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 2.0, 2.0, 2.0, 2.0, 5.0)
  control_params <- list(
    "maxit" = 10000,
    "parscale" = c(1, 1, 1, 1, 1, 1, 0.001, 0.001, 0.001, 0.001, 1)
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
  mu83ms_9deg  <- fit$par[1]
  mu83ms_13deg <- fit$par[2] 
  mu83ms_17deg <- fit$par[3] 
  mu83ms_21deg <- fit$par[4] 
  mu83ms_25deg <- fit$par[5] 
  mu83ms_color <- fit$par[6] 
  lambda117ms <-  fit$par[7] 
  lambda150ms <-  fit$par[8] 
  sigma117ms <-   fit$par[9] 
  sigma150ms <-   fit$par[10]
  theta <-        fit$par[11] 
  logL <-         fit$value
  
  est <- data.frame(mu83ms_9deg   = mu83ms_9deg, 
                    mu83ms_13deg  = mu83ms_13deg, 
                    mu83ms_17deg  = mu83ms_17deg, 
                    mu83ms_21deg  = mu83ms_21deg, 
                    mu83ms_25deg  = mu83ms_25deg, 
                    mu83ms_color = mu83ms_color,
                    lambda117ms = lambda117ms,
                    lambda150ms = lambda150ms,
                    sigma117ms = sigma117ms,
                    sigma150ms = sigma150ms,
                    theta = theta,
                    logL = logL)
  return(list(est, global_predicted_data))
}


### Likelihood function
PFCI_logL <- function(x, inputs) {
  
  # target parameters
  mu83ms_9deg  <- x[1]
  mu83ms_13deg <- x[2] 
  mu83ms_17deg <- x[3] 
  mu83ms_21deg <- x[4] 
  mu83ms_25deg <- x[5] 
  mu83ms_color <- x[6] 
  lambda117ms <-  x[7] 
  lambda150ms <-  x[8] 
  sigma117ms <-   x[9] 
  sigma150ms <-   x[10]
  theta <-        x[11] 
  
  # model predictions
  mu83ms_chimeras <- c(mu83ms_9deg, mu83ms_13deg, mu83ms_17deg,mu83ms_21deg,mu83ms_25deg)
  predicted_data <- matrix(NA, nrow = 21, ncol = 2)
  # gray image
  mean_one_mat <- c(mu83ms_gray, mu83ms_gray*lambda117ms, mu83ms_gray*lambda150ms)
  sigmamat <- c(sigma83ms,sigma117ms,sigma150ms)
  for (cond in 1:3) {
    predicted_data[cond, 1] <- pnorm(theta, mean = mean_one_mat[cond], sd = sigmamat[cond]) # pred_gray_others_rate
    predicted_data[cond, 2] <- 1 - predicted_data[cond, 1] # pred_gray_color_rate
  }
  # chimera image
  mean_one_mat <- c(mu83ms_chimeras, mu83ms_chimeras*lambda117ms, mu83ms_chimeras*lambda150ms)
  sigmamat <- c(rep(sigma83ms,5),rep(sigma117ms,5),rep(sigma150ms,5))
  for (cond in 4:18) {
    predicted_data[cond, 1] <- pnorm(theta, mean = mean_one_mat[cond - 3], sd = sigmamat[cond - 3]) # pred_chimera_others_rate
    predicted_data[cond, 2] <- 1 - predicted_data[cond, 1] # pred_chimera_color_rate
  }
  # full-color image
  mean_two_mat <- c(mu83ms_color,mu83ms_color*lambda117ms,mu83ms_color*lambda150ms)
  sigmamat <- c(sigma83ms,sigma117ms,sigma150ms)
  for (cond in 19:21) {
    predicted_data[cond, 1] <- pnorm(theta, mean = mean_two_mat[cond - 18], sd = sigmamat[cond - 18]) # pred_color_others_rate
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

# ---- Prediction of confidence based on mean+variance model ----

mu83ms_17deg <- mean(estimates$mu83ms_17deg)
mu83ms_color <- mean(estimates$mu83ms_color)
lambda150ms <-  mean(estimates$lambda150ms)
sigma150ms <-   mean(estimates$sigma150ms)
theta <-  mean(estimates$theta) 
theta_conf_no <- theta - 1
theta_conf_yes <- theta + 1

predicted_conf <- matrix(NA, nrow = 4, ncol = 4) #Column1:Color_high, 2:Color_low, 3:No_low, 4:No_high
mean_no_mat <- c(mu83ms_17deg, mu83ms_17deg*lambda150ms)
sigmamat <- c(sigma83ms,sigma150ms)
for (cond in 1:2) {
  predicted_conf[cond, 4] <- pnorm(theta_conf_no, mean = mean_no_mat[cond], sd = sigmamat[cond]) # N_high_rate
  A <- pnorm(theta, mean = mean_no_mat[cond], sd = sigmamat[cond])
  predicted_conf[cond, 3] <- A - predicted_conf[cond, 4] # N_low_rate
  B <- pnorm(theta_conf_yes, mean = mean_no_mat[cond], sd = sigmamat[cond])
  predicted_conf[cond, 2] <- B - A # S_low_rate
  predicted_conf[cond, 1] <- 1 - (predicted_conf[cond, 4] + predicted_conf[cond, 3] + predicted_conf[cond, 2]) # S_high_rate
}

mean_yes_mat <- c(mu83ms_color,mu83ms_color*lambda150ms)
for (cond in 3:4) {
  predicted_conf[cond, 4] <- pnorm(theta_conf_no, mean = mean_yes_mat[cond-2], sd = sigmamat[cond-2]) # N_high_rate
  C <- pnorm(theta, mean = mean_yes_mat[cond-2], sd = sigmamat[cond-2])
  predicted_conf[cond, 3] <- C - predicted_conf[cond, 4] # N_low_rate
  D <- pnorm(theta_conf_yes, mean = mean_yes_mat[cond-2], sd = sigmamat[cond-2])
  predicted_conf[cond, 2] <- D - C # S_low_rate
  predicted_conf[cond, 1] <- 1 - (predicted_conf[cond, 4] + predicted_conf[cond, 3] + predicted_conf[cond, 2]) # S_high_rate
  
}

index_order <- c(
  1, 3, 
  2, 4 
)
predicted_conf <- predicted_conf[index_order,] * 100

### Estimated distribution
plot_sdt_distributions_conf <- function(means, sds, attention_levels, image_types, colors) {
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

  Distribution_conf <- ggplot(data_SDT_plot, aes(x = x, y = y, color = ImageType)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = colors) +
    labs(x = "Strength of peripheral color signal",
         y = "Probability density") +
    scale_x_continuous(breaks = seq(-2, 8, length = 6),limits = c(-2, 8)) +
    scale_y_continuous(breaks = seq(0, 0.6, length = 4),limits = c(0, 0.6)) +
    geom_vline(xintercept = theta, linetype = "dashed", color = "black") +
    geom_vline(xintercept = theta_conf_no, linetype = "dashed", color = "#1f78b4") +
    geom_vline(xintercept = theta_conf_yes, linetype = "dashed", color = "#ffcc00") +
    facet_wrap(~ Attention, nrow = 3, scales = "free_y") +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(size = 0.5, color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.3, "cm"),
      legend.position = "none"
    )

    ggsave(file = "Distribution_conf.png", plot = Distribution_conf, dpi = 150, width = 8, height = 6)

}

means_83ms <- c(mu83ms_17deg, mu83ms_color)
means_150ms <- means_83ms * lambda150ms
means <- rbind(means_83ms, means_150ms)
sd_83ms <-  c(rep(sigma83ms, 2))
sd_150ms <- c(rep(sigma150ms, 2))
sds <- rbind(sd_83ms, sd_150ms)

attention_levels <- factor(c("83 ms", "150 ms"), levels = c("83 ms", "150 ms"))
image_types <- factor(
  c("chimera 17 degree",  "full-color"),
  levels = c("chimera 17 degree",  "full-color"))
colors <- viridis(7, option = "plasma")
colors <- colors[c(4,7)]

plot_sdt_distributions_conf(means, sds, attention_levels, image_types, colors)


response_labels <- c(
  "Yes response with high confidence",
  "Yes response with low confidence",
  "No response with low confidence",
  "No response with high confidence"
)

data_bar_conf <- data.frame(
  Imagetype = rep(c("chimera 17 degree","full-color"), each = 4, times = 2),
  Condition = rep(c("83 ms", "150 ms"), each = 8),
  ResponseType = rep(response_labels, times = 2),
  Proportion = as.vector(t(predicted_conf)) 
)

data_bar_conf$Imagetype <- factor(data_bar_conf$Imagetype,
                                  levels = c("chimera 17 degree", "full-color"))
data_bar_conf$Condition <- factor(data_bar_conf$Condition,
                                  levels = c("83 ms", "150 ms"))
data_bar_conf$ResponseType <- factor(data_bar_conf$ResponseType,
                                     levels = response_labels)


bar_graph_conf <- ggplot(data_bar_conf, aes(x = Imagetype, y = Proportion, fill = ResponseType)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = NULL, y = "Proportion (%)", fill = NULL) +
  facet_grid(. ~ Condition) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    strip.text = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) +
  scale_fill_manual(
    values = c(
      "Yes response with high confidence" = "#ffcc00",
      "Yes response with low confidence" = "#ffe680",
      "No response with low confidence" = "#a6cee3",
      "No response with high confidence" = "#1f78b4"
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 25),
    limits = c(0, 100)
  )

ggsave(file = "bar_graph_conf.png", plot = bar_graph_conf, dpi = 180, width = 10, height = 6)



# ---- Prediction of confidence based on Bayesian observer ----

mu_S_83ms <- 2.95#mu83ms_color 
mu_N_83ms <- 1.73#mu83ms_17deg 
mu_S_150ms <- mu_S_83ms * 1.379
mu_N_150ms <- mu_N_83ms * 1.379
x_seq <- seq(0, 8, length.out = 4000)

# Posterior Probability_S
calc_posterior <- function(x, mu_S, mu_N, sigma) {
  likelihood_S <- dnorm(x, mean = mu_S, sd = sigma)
  likelihood_N <- dnorm(x, mean = mu_N, sd = sigma)
  posterior <- likelihood_S / (likelihood_S + likelihood_N)
  return(data.frame(x = x, likelihood_S = likelihood_S, likelihood_N = likelihood_N, posterior = posterior))
}

posterior_83ms  <- calc_posterior(x_seq, mu_S_83ms,  mu_N_83ms,  sigma83ms)
posterior_150ms <- calc_posterior(x_seq, mu_S_150ms, mu_N_150ms, sigma83ms)

# Find corresponding signal strength
x_at_post75_83ms <- posterior_83ms$x[which.min(abs(posterior_83ms$posterior - 0.75))]
x_at_post75_150ms <- posterior_150ms$x[which.min(abs(posterior_150ms$posterior - 0.75))]
x_at_post50_83ms <- (mu_S_83ms + mu_N_83ms) / 2
x_at_post50_150ms <- (mu_S_150ms + mu_N_150ms) / 2
x_at_post25_83ms <- posterior_83ms$x[which.min(abs(posterior_83ms$posterior - 0.25))]
x_at_post25_150ms <- posterior_150ms$x[which.min(abs(posterior_150ms$posterior - 0.25))]

# Prediction
posterior_conf <- matrix(NA, nrow = 4, ncol = 4) 
xpost75mat <- c(x_at_post75_83ms,x_at_post75_150ms)
xpost50mat <- c(x_at_post50_83ms,x_at_post50_150ms)
xpost25mat <- c(x_at_post25_83ms,x_at_post25_150ms)

#chimera
mean_N_mat <- c(mu_N_83ms, mu_N_150ms)
for (cond in 1:2) {
  posterior_conf[cond, 4] <- pnorm(xpost25mat[cond], mean = mean_N_mat[cond], sd = sigma83ms) # N_high_rate
  O <- pnorm(xpost50mat[cond], mean = mean_N_mat[cond], sd = sigma83ms)
  posterior_conf[cond, 3] <- O - posterior_conf[cond, 4] # N_low_rate
  P <- pnorm(xpost75mat[cond], mean = mean_N_mat[cond], sd = sigma83ms)
  posterior_conf[cond, 2] <- P - O # S_low_rate
  posterior_conf[cond, 1] <- 1 - (posterior_conf[cond, 4] + posterior_conf[cond, 3] + posterior_conf[cond, 2]) # S_high_rate
}

#full-color
mean_S_mat <- c(mu_S_83ms,mu_S_150ms)
for (cond in 3:4) {
  posterior_conf[cond, 4] <- pnorm(xpost25mat[cond-2], mean = mean_S_mat[cond-2], sd = sigma83ms) # N_high_rate
  Q <- pnorm(xpost50mat[cond-2], mean = mean_S_mat[cond-2], sd = sigma83ms)
  posterior_conf[cond, 3] <- Q - posterior_conf[cond, 4] # N_low_rate
  R <- pnorm(xpost75mat[cond-2], mean = mean_S_mat[cond-2], sd = sigma83ms)
  posterior_conf[cond, 2] <- R - Q # S_low_rate
  posterior_conf[cond, 1] <- 1 - (posterior_conf[cond, 4] + posterior_conf[cond, 3] + posterior_conf[cond, 2]) # S_high_rate

}

posterior_conf <- posterior_conf[index_order,] * 100

# Creating graphs
# Distributions
plot_sdt_distributions_post_conf <- function(means, sds, attention_levels, image_types, colors) {
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
  
  theta_values <- data.frame(
    Attention = factor(c("83 ms", "150 ms"), levels = c("83 ms", "150 ms")),
    theta = xpost50mat,
    theta_conf_no = xpost25mat,
    theta_conf_yes = xpost75mat
  )
  
  Distribution_post_conf <- ggplot(data_SDT_plot, aes(x = x, y = y, color = ImageType)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = colors) +
    labs(x = "Strength of peripheral color signal",
         y = "Probability density") +
    scale_x_continuous(breaks = seq(-2, 8, length = 6),limits = c(-2, 8)) +
    scale_y_continuous(breaks = seq(0, 0.6, length = 4),limits = c(0, 0.6)) +
    geom_vline(data = theta_values, aes(xintercept = theta), linetype = "dashed", color = "black") +
    geom_vline(data = theta_values, aes(xintercept = theta_conf_no), linetype = "dashed", color = "#1f78b4") +
    geom_vline(data = theta_values, aes(xintercept = theta_conf_yes), linetype = "dashed", color = "#ffcc00") +
    facet_wrap(~ Attention, nrow = 3, scales = "free_y") +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(size = 0.5, color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.3, "cm"),
      legend.position = "none"
    )

  ggsave(file = "Distribution_post_conf.png", plot = Distribution_post_conf, dpi = 150, width = 8, height = 6)
  
}

means_83ms <- c(mu_N_83ms, mu_S_83ms)
means_150ms <- c(mu_N_150ms, mu_S_150ms)
means <- rbind(means_83ms, means_150ms)
sd_83ms <-  c(rep(sigma83ms, 2))
sd_150ms <- sd_83ms
sds <- rbind(sd_83ms, sd_150ms)
attention_levels <- factor(c("83 ms", "150 ms"), levels = c("83 ms", "150 ms"))
image_types <- factor(
  c("chimera 17 degree",  "full-color"),
  levels = c("chimera 17 degree",  "full-color"))
colors <- viridis(7, option = "plasma")
colors <- colors[c(4,7)]
plot_sdt_distributions_post_conf(means, sds, attention_levels, image_types, colors)

# Bar graphs

data_bar_post_conf <- data_bar_conf
data_bar_post_conf$Proportion <- as.vector(t(posterior_conf))

bar_graph_conf <- ggplot(data_bar_post_conf, aes(x = Imagetype, y = Proportion, fill = ResponseType)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = NULL, y = "Proportion (%)", fill = NULL) +
  facet_grid(. ~ Condition) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    strip.text = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) +
  scale_fill_manual(
    values = c(
      "Yes response with high confidence" = "#ffcc00",
      "Yes response with low confidence" = "#ffe680",
      "No response with low confidence" = "#a6cee3",
      "No response with high confidence" = "#1f78b4"
    )
  ) 
ggsave(file = "bar_graph_posterior_conf.png", plot = bar_graph_conf, dpi = 180, width = 10, height = 6)

# #Detailed bar graph
# 
# data_summary <- data_bar_post_conf %>%
#   mutate(
#     Confidence = ifelse(grepl("high confidence", ResponseType), "High", "Low"),
#     Correctness = case_when(
#       Imagetype == "chimera 17 degree" & grepl("^No response", ResponseType) ~ "Correct",
#       Imagetype == "full-color" & grepl("^Yes response", ResponseType) ~ "Correct",
#       TRUE ~ "Incorrect"
#     )
#   ) %>%
#   group_by(Imagetype, Condition, Correctness, Confidence) %>%
#   summarise(Proportion = sum(Proportion), .groups = "drop") %>%
#   pivot_wider(names_from = Confidence, values_from = Proportion, values_fill = 0) %>%
#   mutate(
#     Total = High + Low,
#     HighConfidenceRate = High / Total * 100
#   )
# 
# data_summary %>%
#   select(Condition, Imagetype, Correctness, HighConfidenceRate) %>%
#   arrange(Condition, Imagetype, Correctness)
# 
# 
# ggplot(data_summary, aes(x = Imagetype, y = HighConfidenceRate, fill = Correctness)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_grid(. ~ Condition) +
#   labs(
#     x = NULL,
#     y = "High Confidence Rate (%)",
#     fill = "Response Accuracy"
#   ) +
#   scale_fill_manual(
#     values = c("Correct" = "#4CAF50", "Incorrect" = "#F44336")  # 緑：正答、赤：誤答
#   ) +
#   theme_minimal(base_size = 18) +
#   coord_cartesian(ylim = c(0, 100)) +  # ここでy軸の範囲を指定
#   theme(
#     axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
#     axis.text.y = element_text(size = 14),
#     axis.title.y = element_text(size = 16),
#     strip.text = element_text(size = 16),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 13),
#     axis.ticks.length = unit(0.25, "cm"),              # 目盛の長さ
#     axis.ticks = element_line(color = "black"),        # 目盛の線を描く
#     axis.ticks.y.right = element_line(),               # y軸右側も表示したい場合
#     axis.line = element_line(color = "black"),         # 軸線を描く
#     panel.grid.major = element_blank(),                # メジャー罫線を消す
#     panel.grid.minor = element_blank()                 # マイナー罫線を消す
#   )
