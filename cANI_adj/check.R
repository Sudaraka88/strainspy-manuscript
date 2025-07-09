library(ggplot2)
library(patchwork)
library(glmmTMB)

see_boxplot = function(df){
  p1 = ggplot2::ggplot(df, ggplot2::aes(x = spiked, y = value)) + 
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    ggbeeswarm::geom_quasirandom() +
    ggplot2::labs(x = "Spiked", y = "ANI") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::scale_color_brewer(palette = "Set2") +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  p2 = ggplot2::ggplot(df[df$value>0, ], ggplot2::aes(x = spiked, y = value)) + 
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    ggbeeswarm::geom_quasirandom() +
    ggplot2::labs(x = "Spiked", y = "ANI") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::scale_color_brewer(palette = "Set2") +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  patchwork::wrap_plots(p1/p2)
  
}
nbeta = 2
## PRIORS
p1 = data.frame(
  prior = rep("normal(0,5)", 2*nbeta),
  class = rep(c("fixef", "fixef_zi"), each=nbeta),
  coef  = rep(as.character(seq(1,nbeta)), 2))

p2 = data.frame(
  prior = rep("normal(0,0.5)", 2*nbeta),
  class = rep(c("fixef", "fixef_zi"), each=nbeta),
  coef  = rep(as.character(seq(1,nbeta)), 2))

## DATA and prelim
df_fp = readRDS("some_data.rds")
df_tp = readRDS("signal_data.rds")
df_null = df_fp; df_null$spiked = sample(df_fp$spiked)

df_tp$value_ = df_tp$value/100
df_fp$value_ = df_fp$value
df_null$value_ = df_null$value
df_fp$spiked = ifelse(df_fp$spiked == 1, "No", "Yes")
df_null$spiked = ifelse(df_null$spiked == 1, "No", "Yes")



# prelim setup
# df_tp$value = base::pmin(as.vector(df_tp$value) / 100, 0.99999)


# df = rbind(df, data.frame(value = 0.99, spiked = 'Yes'))
# res = function(x){(x-min(x))/(max(x)-min(x))*0.5}
# # 
# df$value[df$value != 0] = res(df$value[df$value != 0])

see_boxplot(df_fp)
see_boxplot(df_tp)
see_boxplot(df_null)

reset_1 = function(x, v = 0.99999){
  base::pmin(x, 0.99999)
}

squash = function(x, eps = 1e-2) {
  x[x > 0] <- x[x > 0] * (1 - 2 * eps) + eps
  # x[x>0] <- x[x>0] - eps # offset
  x
}

offset = function(x, eps = 1e-2){
  x[x>0] <- x[x>0] - eps
  x
}

rescale_p01p99 <- function(x) {
  out <- x
  is_nonzero <- x != 0
  x_nz <- x[is_nonzero]
  
  rng <- range(x_nz, na.rm = TRUE)
  rescaled <- (x_nz - rng[1]) / (rng[2] - rng[1])
  out[is_nonzero] <- rescaled * 0.98 + 0.01
  
  return(out)
}


map_range <- function(x, min_val = 0.9) {
  out <- x
  is_nonzero <- x > 0
  out[is_nonzero] <- (x[is_nonzero] - min_val) / (1-min_val) * 0.98 + 0.01
  return(out)
}

compare_methods = function(x){
  plt_df = data.frame(val_org = rep(x,4),
                      method = c(rep('squash',length(x)), rep('offset',length(x)), rep('rescale',length(x)), rep('map',length(x))),
                      val_adj = c(squash(x), offset(x), rescale_p01p99(x), map_range(x)))
  
  ggplot(plt_df[plt_df$val_org!=0, ]) + 
    geom_point(aes(x = val_org, y = val_adj, col = method)) + 
    facet_wrap(~method, scales = "free_y", nrow = 1) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 20)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5))
  
}

p1 = compare_methods(c(0, 9000:10000/1e4))
p2 = compare_methods(c(0, 925:1000/1000))
p3 = compare_methods(c(0, 0.95, 0.96, 0.97, 0.98, 0.99, 1))
patchwork::wrap_plots(p1/p2/p3)

# Let's try adjusting and check outcomes

combined_formula = as.formula('value~spiked')
design = as.formula('~spiked')


fit_and_summary_zib = function(df, combined_formula = as.formula('value~spiked'), design = as.formula('~spiked')){
  x = summary(glmmTMB(formula = combined_formula, ziformula = design, data = df, family = glmmTMB::beta_family(link = "logit")))
  tmp = c(x$coefficients$cond[2, c(1,4)], x$coefficients$zi[2, c(1,4)])
  names(tmp) = c("beta", "p_beta", "z_beta", "z_p_beta")
  tmp
}

fit_and_summary_ob = function(df, combined_formula = as.formula('value~spiked')){
  x = summary(glmmTMB(formula = combined_formula, data = df, family = glmmTMB::ordbeta()))
  tmp = x$coefficients$cond[2, c(1,4)]
  names(tmp) = c("beta", "p_beta")
  tmp
}

fit_and_summary_qb = function(df, combined_formula = as.formula('value~spiked')){
  x = summary(stats::glm(formula = combined_formula, data = df, family = 'quasibinomial'))
  tmp = x$coefficients[2,c(1,4)]
  names(tmp) = c("beta", "p_beta")
  tmp
}

fit_all <- function(df, dpoint){
  op <- data.frame()
  
  adj_methods <- list(
    reset_1 = reset_1,
    squash = squash,
    offset = offset,
    rescale = rescale_p01p99,
    map = map_range
  )
  
  models <- list(
    zib = fit_and_summary_zib,
    ob  = fit_and_summary_ob,
    qb  = fit_and_summary_qb
  )
  for (model_name in names(models)) {
    fit_fun <- models[[model_name]]
    
    for (adj_name in names(adj_methods)) {
      adj_fun <- adj_methods[[adj_name]]
      
      df$value <- adj_fun(df$value_)
      # see_boxplot(df_tp)
      
      tmp <- fit_fun(df)
      op <- rbind(op, data.frame(
        dpoint = dpoint,
        model = model_name,
        adj = adj_name,
        stat = names(tmp),
        vals = unname(tmp)
      ))
    }
  }
  
  op
}

op = rbind(fit_all(df_tp, 'tp'),
           fit_all(df_fp, 'fp'),
           fit_all(df_null, 'null'))


op$log_vals <- with(op, ifelse(stat %in% c("p_beta", "z_p_beta"), -log10(vals), vals))
op$dpoint = factor(op$dpoint, levels = c("tp", "fp","null"))
op$model = factor(op$model, levels = c("zib", "ob","qb"))
op$adj = factor(op$adj, levels = c("reset_1" , "squash", "offset", "map", "rescale"))

ggplot(op, aes(x = dpoint, y = ifelse(stat %in% c("p_beta", "z_p_beta"), -log10(vals), vals), colour = dpoint)) +
  geom_point() +
  facet_grid(stat ~ adj + model, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()

# False positive
# df_fp$value = squash(df_fp$value_); see_boxplot(df_fp); fit_and_summary_zib(df_fp)
# df_fp$value = offset(df_fp$value_); see_boxplot(df_fp); fit_and_summary_zib(df_fp)
# df_fp$value = rescale_p01p99(df_fp$value_); see_boxplot(df_fp); fit_and_summary_zib(df_fp)
# df_fp$value = map_range(df_fp$value_); see_boxplot(df_fp); fit_and_summary_zib(df_fp)
# 
# # True Negative
# df_null$value = squash(df_null$value_); see_boxplot(df_null); fit_and_summary_zib(df_null)
# df_null$value = offset(df_null$value_); see_boxplot(df_null); fit_and_summary_zib(df_null)
# df_null$value = rescale_p01p99(df_null$value_); see_boxplot(df_null); fit_and_summary_zib(df_null)
# df_null$value = map_range(df_null$value_); see_boxplot(df_null); fit_and_summary_zib(df_null)
# 
# 
# df_tp_offset = df_tp; df_tp_offset$value[df_tp_offset$value > 0] = df_tp_offset$value[df_tp_offset$value > 0] - 0.02
# see_boxplot(df_tp_offset)
# summary(glmmTMB(formula = combined_formula, ziformula = design, data = df_tp_offset, family = glmmTMB::beta_family(link = "logit")))
# 
# # Shifting the false positive
# see_boxplot(df_fp)
# summary(glmmTMB(formula = combined_formula, ziformula = design, data = df_fp, family = glmmTMB::beta_family(link = "logit")))
# 
# df_fp_offset = df_fp; df_fp_offset$value[df_fp_offset$value > 0] = df_fp_offset$value[df_fp_offset$value > 0] - 0.02
# see_boxplot(df_fp_offset)
# summary(glmmTMB(formula = combined_formula, ziformula = design, data = df_fp_offset, family = glmmTMB::beta_family(link = "logit")))

# Shifting the null

# 
# 
# summary(glmmTMB(formula = combined_formula, ziformula = design, data = df, family = glmmTMB::beta_family(link = "logit")))
# summary(glmmTMB(formula = combined_formula, ziformula = design, data = df, family = glmmTMB::beta_family(link = "logit"), priors = p1))
# summary(glmmTMB(formula = combined_formula, ziformula = design, data = df, family = glmmTMB::beta_family(link = "logit"), priors = p2))
# 
# 
# ## NULL DATA
# 
# see_boxplot(df_null)
# summary(glmmTMB(formula = combined_formula, ziformula = design, data = df_null, family = glmmTMB::beta_family(link = "logit")))
# summary(glmmTMB(formula = combined_formula, ziformula = design, data = df_null, family = glmmTMB::beta_family(link = "logit"), priors = p1))
# summary(glmmTMB(formula = combined_formula, ziformula = design, data = df_null, family = glmmTMB::beta_family(link = "logit"), priors = p2))
# 
# # quasi binomial
# summary(glm(formula = combined_formula, data = df,family = 'quasibinomial'))
# summary(glm(formula = combined_formula, data = df_null,family = 'quasibinomial'))
