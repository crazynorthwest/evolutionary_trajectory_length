library(ggplot2)
library(car)
library(lme4)
library(patchwork)
library(dplyr)
library(rstatix)
library(ggpubr)

plot_mean_errorbar <- function(data, x_col,  y_col,  x_label, y_label) {
  x_levels <- sort(unique(data[[x_col]]))
  summary_df <- data %>%
    filter(
      !is.na(.data[[x_col]]),
      !is.na(.data[[y_col]]),
      .data[[x_col]] %in% x_levels
    ) %>%
    mutate(
      x_factor = factor(.data[[x_col]], levels = x_levels)
    ) %>%
    group_by(x_factor) %>%
    summarise(
      n = n(),
      mean_value = mean(.data[[y_col]], na.rm = TRUE),
      sd_value = sd(.data[[y_col]], na.rm = TRUE),
      se_value = sd_value / sqrt(n),
      .groups = "drop"
    )
  
  ggplot(summary_df, aes(x = x_factor, y = mean_value, group = 1)) +
    geom_point(size = 2.8) +
    geom_line(linewidth = 0.8) +
    geom_errorbar(
      aes(
        ymin = mean_value - 1.96 *se_value,
        ymax = mean_value + 1.96 *se_value
      ),
      width = 0.12,
      linewidth = 0.6
    ) +
    labs(
      x = x_label,
      y = y_label
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.6
      )
    )
}


generation_record <-read.csv(paste0(path_figures,"/generation_record.dat"), header = TRUE, sep = " ")
mut_fix_sum <-read.csv(paste0(path_figures,"/mut_fix_sum.dat"), header = TRUE, sep = " ")

generation_record_sub <- subset(generation_record, par == 1)
mut_fix_sum_sub <- subset(mut_fix_sum, par == 1)

data_clean <- generation_record_sub$generation_record[!is.na(generation_record_sub$generation_record)]
source(path_selection)
source(path_mu_selection)

figure_hist <- ggplot(data.frame(x = data_clean), aes(x = x)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black") +
  labs(title = "", x = "Generation", y = "Frequency") +
  theme_minimal()
figure_hist

figure_ccdf <- ggplot() +
  geom_point(
    data = emp_df,
    aes(x = x, y = ccdf, color = model),
    size = 1.6,
    alpha = 0.85
  ) +
  geom_line(
    data = fit_df,
    aes(
      x = x,
      y = ccdf,
      color = model,
      linetype = model,
      linewidth = model
    )
  ) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(
    xlim = c(xmin, xmax),
    ylim = c(y_min, 1)
  ) +
  scale_color_manual(
    name = NULL,
    values = plot_cols,
    breaks = legend_order
  ) +
  scale_linetype_manual(
    name = NULL,
    values = plot_ltys,
    breaks = legend_order
  ) +
  scale_linewidth_manual(
    name = NULL,
    values = plot_lwds,
    breaks = legend_order
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = c(16, rep(NA, length(fit_models))),
        linetype = unname(plot_ltys),
        linewidth = unname(plot_lwds),
        size = 1.2
      ),
      keyheight = unit(0.35, "cm"),
      keywidth = unit(0.7, "cm"),
      label.theme = element_text(size = 5)
    ),
    linetype = "none",
    linewidth = "none"
  ) +
  labs(
    x = "Generation",
    y = "P(generation >= x)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = c(0.05, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.75),
      color = "black",
      linewidth = 0.2
    ),
    legend.key = element_blank(),
    legend.key.width = unit(1.2, "cm"),
    legend.text = element_text(size = 6),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.2
    )
  )

print(figure_ccdf)


figure_generation_number_mean <- plot_mean_errorbar(
  data = generation_record,
  x_col = "fix_mut_number",
  y_col = "generation_record",
  x_label = "Number of Fixed Mutation",
  y_label = "Generation"
)

print(figure_generation_number_mean)

figure_generation_selection_mean <- plot_mean_errorbar(
  data = generation_record,
  x_col = "par",
  y_col = "generation_record",
  x_label = "Selection Strength",
  y_label = "Generation"
)

print(figure_generation_selection_mean)

figure_mut_size_selection_mean <- plot_mean_errorbar(
  data = mut_fix_sum,
  x_col = "par",
  y_col = "mut_size",
  x_label = "Selection Strength",
  y_label = "Effect Size"
)

print(figure_mut_size_selection_mean)

figure_mu <- ggplot(mu_table, aes(x = parameter, y = mu)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mu_lower, ymax = mu_upper),
    width = 0.04,
    linewidth = 0.7
  ) +
  scale_x_continuous(
    breaks = sort(unique(mu_table$parameter)),
    labels = function(x) sub("\\.?0+$", "", as.character(x))
  ) +
  labs(
    x = "Selection Strength",
    y = expression(mu)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    )
  )
figure_mu

mut_generation_mean <- aggregate(
  mut_size ~ generation_record,
  data = mut_fix_sum_sub,
  FUN = mean
)

names(mut_generation_mean)[2] <- "mut_size_mean"

fit_log_generation_mutSize <- lm(
  log(generation_record) ~ mut_size_mean,
  data = mut_generation_mean
)
summary(fit_log_generation_mutSize)
pred_log_generation_mutSize <- data.frame(
  mut_size_mean = seq(
    min(mut_generation_mean$mut_size_mean),
    max(mut_generation_mean$mut_size_mean),
    length.out = 200
  )
)

pred_log <- predict(
  fit_log_generation_mutSize,
  newdata = pred_log_generation_mutSize,
  interval = "confidence"
)

pred_log_generation_mutSize$fit_generation <- exp(pred_log[, "fit"])
pred_log_generation_mutSize$lower_generation <- exp(pred_log[, "lwr"])
pred_log_generation_mutSize$upper_generation <- exp(pred_log[, "upr"])

eq_label <- paste(
  sprintf(
    "log(y) = %.3f %+.3f·x",
    coef(fit_log_generation_mutSize)[1],
    coef(fit_log_generation_mutSize)[2]
  ),
  sprintf(
    "R² = %.3f",
    summary(fit_log_generation_mutSize)$r.squared
  ),
  sep = "\n"
)

figure_generation_mutSize <- ggplot(
  mut_generation_mean,
  aes(x = mut_size_mean, y = generation_record)
) +
  geom_point(
    alpha = 0.75,
    color = "#2C7FB8",
    size = 2.4
  ) +
  geom_ribbon(
    data = pred_log_generation_mutSize,
    aes(
      x = mut_size_mean,
      ymin = lower_generation,
      ymax = upper_generation
    ),
    inherit.aes = FALSE,
    fill = "#FDB863",
    alpha = 0.28
  ) +
  geom_line(
    data = pred_log_generation_mutSize,
    aes(x = mut_size_mean, y = fit_generation),
    inherit.aes = FALSE,
    color = "#D7301F",
    linewidth = 1.2
  ) +
  scale_y_log10(
    breaks = c(100, 1000, 5000),
    labels = c("100", "1,000", "5,000")
  ) +
  annotate(
    "label",
    x = max(mut_generation_mean$mut_size_mean, na.rm = TRUE) * 0.65,
    y = max(mut_generation_mean$generation_record, na.rm = TRUE) * 0.9,
    label = eq_label,
    hjust = 0,
    vjust = 1,
    size = 3.3,
    color = "white",
    fill = "#333333",
    alpha = 0.85,
    label.size = 0,
    fontface = "bold"
  ) +
  labs(
    x = "Mean Effect Size of Fixed Mutations",
    y = "Generation"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    )
  )

figure_generation_mutSize

fit_log_generationFix_mutSize <- lm(
  log(fix_time) ~ mut_size + I(mut_size^2),
  data = mut_fix_sum_sub
)

summary(fit_log_generationFix_mutSize)

pred_log_generationFix_mutSize <- data.frame(
  mut_size = seq(
    min(mut_fix_sum_sub$mut_size, na.rm = TRUE),
    max(mut_fix_sum_sub$mut_size, na.rm = TRUE),
    length.out = 200
  )
)

pred_log <- predict(
  fit_log_generationFix_mutSize,
  newdata = pred_log_generationFix_mutSize,
  interval = "confidence"
)

pred_log_generationFix_mutSize$fit <- exp(pred_log[, "fit"])
pred_log_generationFix_mutSize$lwr <- exp(pred_log[, "lwr"])
pred_log_generationFix_mutSize$upr <- exp(pred_log[, "upr"])

eq_label <- paste(
  sprintf(
    "log(y) = %.3f %+.3f·x %+.3f·x²",
    coef(fit_log_generationFix_mutSize)[1],
    coef(fit_log_generationFix_mutSize)[2],
    coef(fit_log_generationFix_mutSize)[3]
  ),
  sprintf(
    "R² = %.3f",
    summary(fit_log_generationFix_mutSize)$r.squared
  ),
  sep = "\n"
)

figure_fix_time_mutSize <- ggplot(
  mut_fix_sum_sub,
  aes(x = mut_size, y = fix_time)
) +
  geom_point(
    color = "#2C7FB8",
    size = 2.5,
    alpha = 0.75
  ) +
  geom_ribbon(
    data = pred_log_generationFix_mutSize,
    aes(
      x = mut_size,
      ymin = lwr,
      ymax = upr
    ),
    inherit.aes = FALSE,
    fill = "#FDB863",
    alpha = 0.30
  ) +
  geom_line(
    data = pred_log_generationFix_mutSize,
    aes(
      x = mut_size,
      y = fit
    ),
    inherit.aes = FALSE,
    color = "#D7301F",
    linewidth = 1.2
  ) +
  scale_y_log10(
    breaks = c(100, 1000, 5000),
    labels = c("100", "1,000", "5,000")
  ) +
  annotate(
    "label",
    x = max(mut_fix_sum_sub$mut_size, na.rm = TRUE) * 0.5,
    y = max(mut_fix_sum_sub$fix_time, na.rm = TRUE) * 0.9,
    label = eq_label,
    hjust = 0,
    vjust = 1,
    size = 3.3,
    color = "white",
    fill = "#333333",
    alpha = 0.85,
    label.size = 0,
    fontface = "bold"
  ) +
  labs(
    x = "Effect Size",
    y = "Generation"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    )
  )

figure_fix_time_mutSize

fitness_mean <- aggregate(
  cbind(fitness_final, fitness_f1, fitness_f2) ~ repeat_number,
  data = generation_record_sub,
  FUN = function(x) mean(x, na.rm = TRUE)
)

generation_sum <- aggregate(
  generation_record ~ repeat_number,
  data = generation_record_sub,
  FUN = function(x) sum(x, na.rm = TRUE)
)

names(generation_sum)[2] <- "generation_sum"

result_final_f1_f2 <- merge(
  fitness_mean,
  generation_sum,
  by = "repeat_number"
)

result_final_f1_f2$final_f1 <- result_final_f1_f2$fitness_final - result_final_f1_f2$fitness_f1
result_final_f1_f2$final_f2 <- result_final_f1_f2$fitness_final - result_final_f1_f2$fitness_f2

fit_generation_sum_final_f2_1 <- lm(
  final_f2 ~ generation_sum,
  data = result_final_f1_f2
)

summary(fit_generation_sum_final_f2_1)
Anova(fit_generation_sum_final_f2_1)

pred <- data.frame(
  generation_sum = seq(
    min(result_final_f1_f2$generation_sum, na.rm = TRUE),
    max(result_final_f1_f2$generation_sum, na.rm = TRUE),
    length.out = 200
  )
)

pred_ci <- predict(
  fit_generation_sum_final_f2_1,
  newdata = pred,
  interval = "confidence"
)

pred$fit <- pred_ci[, "fit"]
pred$lwr <- pred_ci[, "lwr"]
pred$upr <- pred_ci[, "upr"]

eq_label <- paste(
  sprintf(
    "y = %.3f %+.7f·x",
    coef(fit_generation_sum_final_f2_1)[1],
    coef(fit_generation_sum_final_f2_1)[2]
  ),
  sprintf(
    "R² = %.3f",
    summary(fit_generation_sum_final_f2_1)$r.squared
  ),
  sep = "\n"
)

figure_generation_sum_final_f2 <- ggplot(
  result_final_f1_f2,
  aes(x = generation_sum, y = final_f2)
) +
  geom_point(
    color = "#2C7FB8",
    size = 2.5,
    alpha = 0.75
  ) +
  geom_ribbon(
    data = pred,
    aes(
      x = generation_sum,
      ymin = lwr,
      ymax = upr
    ),
    inherit.aes = FALSE,
    fill = "#FDB863",
    alpha = 0.30
  ) +
  geom_line(
    data = pred,
    aes(
      x = generation_sum,
      y = fit
    ),
    inherit.aes = FALSE,
    color = "#D7301F",
    linewidth = 1.2
  ) +
  annotate(
    "label",
    x = max(result_final_f1_f2$generation_sum, na.rm = TRUE) * 0.5,
    y = max(result_final_f1_f2$final_f2, na.rm = TRUE) * 0.8,
    label = eq_label,
    hjust = 0,
    vjust = 1,
    size = 3.3,
    color = "white",
    fill = "#333333",
    alpha = 0.85,
    label.size = 0,
    fontface = "bold"
  ) +
  labs(
    x = "Sum of Generation",
    y = "Fitness Difference"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    )
  )

figure_generation_sum_final_f2


(figure_hist + figure_ccdf) + plot_annotation(tag_levels = "A")
ggsave(filename = paste0(path_figures,"/figure1.jpg"),dpi = 600, width=8, height=4)
figure_generation_number_mean / figure_generation_mutSize / figure_fix_time_mutSize + plot_annotation(tag_levels = "A")
ggsave(filename = paste0(path_figures,"/figure2.jpg"),dpi = 600, width=7, height=9)
figure_mu / figure_generation_selection_mean / figure_mut_size_selection_mean + plot_annotation(tag_levels = "A")
ggsave(filename = paste0(path_figures,"/figure3.jpg"),dpi = 600, width=7, height=9)
print(figure_generation_sum_final_f2)
ggsave(filename = paste0(path_figures,"/figure4.jpg"),dpi = 600, width=5, height=5)


 figure_generation_mutSize / figure_fix_time_mutSize / 
   (figure_generation_number_mean + figure_mu) / 
   (figure_generation_selection_mean + figure_mut_size_selection_mean)/
  figure_generation_sum_final_f2 +
  plot_annotation(tag_levels = "A")
ggsave(filename = paste0(path_figures,"/figureS.jpg"),dpi = 600, width=7, height=10)


summary(fit_log_generationFix_mutSize)
summary(fit_log_generation_mutSize)
summary(fit_generation_sum_final_f2_1)
Anova(fit_log_generationFix_mutSize)
Anova(fit_generation_sum_final_f2_1)
Anova(fit_log_generation_mutSize)
