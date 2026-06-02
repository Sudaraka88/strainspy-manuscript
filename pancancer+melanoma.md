Prediction of immune therapy outcomes in pan cancer: pooling Ash’s rare
cancer and Kshitij’s melanoma data
================
2025-12-18

## Load dependencies

``` r
library(strainspy)
library(SummarizedExperiment)
library(caret)
library(glmnet) # elastic net
library(ranger) # RF
library(pROC)
library(doParallel)
library(foreach)
library(parallel)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggtree)
library(ggsci)
```

## Load metadata

``` r
meta_path <- "./data/ash_pancancer/metadata_full.tsv"
meta_pan <- read.csv(meta_path, sep = '\t') 
meta_pan = cbind(run_acc = meta_pan$run_accession, meta_pan)
# From paper:
# RvsP = CR or PR VS. PD or cPD - excluded patients with a BOR of stable disease (SD) (n=29)

# Outcome
meta_pan$RvsP = "R"
meta_pan$RvsP[which(meta_pan$BOR == "PD" | meta_pan$BOR == "cPD")] = "NR" 
meta_pan$RvsP = factor(meta_pan$RvsP, levels = c("NR", "R"))

# Melanoma meta
meta_m = read.csv("data/melanoma_pooled/Clean_metadata_Aug5_Allyson_match_2025.txt", sep = '\t')


meta = rbind(data.frame(X = meta_pan$run_acc, c_type = "RARE", type = meta_pan$histology_cohort.x, RvsP = meta_pan$RvsP), 
             data.frame(X = meta_m$X, c_type = "MEL", type = meta_m$Study_simplified, RvsP = meta_m$ORR))
```

## Load sylph outputs

``` r
sy <- read_sylph("./data/ash_pancancer/merged_pan_melanoma.tsv.gz") # q99)
```

    ## Detected Sylph query output file.

``` r
# annoying renames to match meta V sylph file
colnames(sy) <- gsub("_1", "", colnames(sy))
colData(sy)$Sample_file <- gsub("_1", "", basename(colData(sy)$Sample_file))

# Reorder meta 
meta = meta[match(colnames(sy), meta$X), ]

sy <- filter_by_presence(sy, min_nonzero = 53) # filter at 10%
```

    ## Retained 19691 rows after filtering

``` r
dim(sy)
```

    ## [1] 19691   526

``` r
# Checks before merging metadata
all(colnames(sy) %in% meta$X)
```

    ## [1] TRUE

``` r
all(meta$X %in% colnames(sy))
```

    ## [1] TRUE

``` r
sy = modify_metadata(sy, meta, replace = T)
```

## Fit the univariate model

``` r
design <- as.formula("~ RvsP")

save_path <- "output_rds/ASH+mela_zib_q_99_ebp.rds"

if(file.exists(save_path)){
  ZB_fit <- readRDS(save_path)
} else {
  # Run with weak prior for prediction
  ebp = compute_eb_priors(sy, design, nthreads = 10L, low_cutoff = 0, high_cutoff = Inf)
  ZB_fit <- glmZiBFit(sy, design, nthreads = parallel::detectCores(), MAP_prior = ebp)
  saveRDS(ZB_fit, save_path)
}
```

# Prediction

## Using Elastic net

### Set up pipeline

``` r
run_enet <- function(train_idx, test_idx, meta, sy, feature_names, return_fit = F) {
  
  train_mx <- strainspy::prep_for_prediction(
    sy[, train_idx],
    'RvsP',
    feature_names
  )
  
  test_mx <- strainspy::prep_for_prediction(
    sy[, test_idx],
    'RvsP',
    feature_names
  )
  
  # Skip if test labels are all one class
  if (length(unique(test_mx$RvsP)) < 2) {
    warning("Test class only has one type of label")
    return(NA)
  }
  
  set.seed(1988)
  fit <- caret::train(
    RvsP ~ .,
    data = train_mx,
    method = "glmnet",
    preProcess = c("center", "scale"),
    metric = "ROC",
    trControl = caret::trainControl(
      method = "cv",
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      savePredictions = "final"
    )
  )
  
  if(return_fit) {
    # If return_fit, we only use the test_idx to see if it has both classes
    return(fit)
  }
  
  preds <- predict(fit, test_mx, type = "prob")$R
  roc_obj <- roc(test_mx$RvsP, preds, levels = c("NR","R"))
  as.numeric(auc(roc_obj))
}

print_auc <- function(boot_auc, ret = F){
  point_est <- median(boot_auc)
  ci_lower  <- quantile(boot_auc, 0.025)
  ci_upper  <- quantile(boot_auc, 0.975)
  
  cat("AUC =", round(point_est, 3), "[95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), "]\n")
  
  if(ret){
    return(c(point_est, ci_lower, ci_upper))
  }
}

# Ordered top hits

th = top_hits(ZB_fit, alpha = 1)
```

    ## Found 19691 tophits for RvsPR at alpha = 1 using holm

``` r
th$min_p <- pmin(th$p_value, th$zi_p_value)
th <- th[order(th$min_p), ]
```

## Train and test on random subsets

    ## AUC = 0.635 [95% CI: 0.43 - 0.797 ]

## Train on rare cancers

### Predict all melanoma

    ## AUC = 0.634 [95% CI: 0.584 - 0.69 ]

### Predict each melanoma dataset

    ## Frankel :

    ## AUC = 0.778 [95% CI: 0.595 - 0.918 ]
    ## Gopalakrishnan :

    ## AUC = 0.66 [95% CI: 0.419 - 0.896 ]
    ## Lee :

    ## AUC = 0.605 [95% CI: 0.47 - 0.726 ]
    ## Matson :

    ## AUC = 0.468 [95% CI: 0.269 - 0.659 ]
    ## McCulloch :

    ## AUC = 0.642 [95% CI: 0.499 - 0.764 ]
    ## Spencer :

    ## AUC = 0.612 [95% CI: 0.521 - 0.701 ]

## Train on melanoma

### Predict all rare cancers

    ## AUC = 0.722 [95% CI: 0.614 - 0.818 ]

### Predict each rare cancer

    ## GYN :

    ## AUC = 0.685 [95% CI: 0.481 - 0.851 ]
    ## NEN :

    ## AUC = 0.777 [95% CI: 0.596 - 0.922 ]
    ## UGB :

    ## AUC = 0.761 [95% CI: 0.571 - 0.92 ]

## Viasualise

``` r
df_all <- bind_rows(
  data.frame(AUC = boot_auc_tRandpRand, scenario = "T-Rand,P-Rand"),
  data.frame(AUC = boot_auc_tRCpMel, scenario = "T-Rare,P-MEL"),
  data.frame(AUC = boot_auc_tMelpRC, scenario = "T-MEL,P-RARE")
)

ggplot(df_all, aes(x = scenario, y = AUC, fill = scenario)) +
  # geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", color = "red", size = 2) +
  ylab("AUC") +
  xlab("") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
```

![](pancancer+melanoma_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# Train and Test on each

## Logic

    ## Training set: RARE - GYN

    ## Testing set: RARE - GYN

    ## AUC = 1 [95% CI: 1 - 1 ]
    ## Testing set: RARE - NEN

    ## AUC = 0.555 [95% CI: 0.299 - 0.783 ]
    ## Testing set: RARE - UGB

    ## AUC = 0.771 [95% CI: 0.607 - 0.901 ]
    ## Testing set: MEL - Frankel

    ## AUC = 0.66 [95% CI: 0.444 - 0.839 ]
    ## Testing set: MEL - Gopalakrishnan

    ## AUC = 0.675 [95% CI: 0.429 - 0.891 ]
    ## Testing set: MEL - Lee

    ## AUC = 0.552 [95% CI: 0.431 - 0.675 ]
    ## Testing set: MEL - Matson

    ## AUC = 0.512 [95% CI: 0.328 - 0.695 ]
    ## Testing set: MEL - McCulloch

    ## AUC = 0.583 [95% CI: 0.434 - 0.721 ]
    ## Testing set: MEL - Spencer

    ## AUC = 0.573 [95% CI: 0.486 - 0.666 ]
    ## Testing set: RARE

    ## AUC = 0.813 [95% CI: 0.726 - 0.891 ]
    ## Testing set: MEL

    ## AUC = 0.599 [95% CI: 0.542 - 0.653 ]
    ## Training set: RARE - NEN

    ## Testing set: RARE - GYN

    ## AUC = 0.513 [95% CI: 0.309 - 0.7 ]
    ## Testing set: RARE - NEN

    ## AUC = 1 [95% CI: 1 - 1 ]
    ## Testing set: RARE - UGB

    ## AUC = 0.637 [95% CI: 0.446 - 0.812 ]
    ## Testing set: MEL - Frankel

    ## AUC = 0.751 [95% CI: 0.552 - 0.901 ]
    ## Testing set: MEL - Gopalakrishnan

    ## AUC = 0.726 [95% CI: 0.487 - 0.921 ]
    ## Testing set: MEL - Lee

    ## AUC = 0.584 [95% CI: 0.454 - 0.701 ]
    ## Testing set: MEL - Matson

    ## AUC = 0.543 [95% CI: 0.329 - 0.76 ]
    ## Testing set: MEL - McCulloch

    ## AUC = 0.604 [95% CI: 0.465 - 0.748 ]
    ## Testing set: MEL - Spencer

    ## AUC = 0.602 [95% CI: 0.507 - 0.691 ]
    ## Testing set: RARE

    ## AUC = 0.744 [95% CI: 0.641 - 0.831 ]
    ## Testing set: MEL

    ## AUC = 0.638 [95% CI: 0.586 - 0.696 ]
    ## Training set: RARE - UGB

    ## Testing set: RARE - GYN

    ## AUC = 0.766 [95% CI: 0.591 - 0.911 ]
    ## Testing set: RARE - NEN

    ## AUC = 0.552 [95% CI: 0.295 - 0.775 ]
    ## Testing set: RARE - UGB

    ## AUC = 1 [95% CI: 1 - 1 ]
    ## Testing set: MEL - Frankel

    ## AUC = 0.608 [95% CI: 0.395 - 0.804 ]
    ## Testing set: MEL - Gopalakrishnan

    ## AUC = 0.594 [95% CI: 0.372 - 0.831 ]
    ## Testing set: MEL - Lee

    ## AUC = 0.599 [95% CI: 0.479 - 0.715 ]
    ## Testing set: MEL - Matson

    ## AUC = 0.4 [95% CI: 0.222 - 0.59 ]
    ## Testing set: MEL - McCulloch

    ## AUC = 0.61 [95% CI: 0.47 - 0.749 ]
    ## Testing set: MEL - Spencer

    ## AUC = 0.609 [95% CI: 0.52 - 0.693 ]
    ## Testing set: RARE

    ## AUC = 0.836 [95% CI: 0.758 - 0.911 ]
    ## Testing set: MEL

    ## AUC = 0.604 [95% CI: 0.548 - 0.659 ]
    ## Training set: MEL - Frankel

    ## Testing set: RARE - GYN

    ## AUC = 0.505 [95% CI: 0.385 - 0.628 ]
    ## Testing set: RARE - NEN

    ## AUC = 0.688 [95% CI: 0.524 - 0.851 ]
    ## Testing set: RARE - UGB

    ## AUC = 0.563 [95% CI: 0.391 - 0.732 ]
    ## Testing set: MEL - Frankel

    ## AUC = 0.881 [95% CI: 0.729 - 0.975 ]
    ## Testing set: MEL - Gopalakrishnan

    ## AUC = 0.542 [95% CI: 0.5 - 0.65 ]
    ## Testing set: MEL - Lee

    ## AUC = 0.469 [95% CI: 0.35 - 0.59 ]
    ## Testing set: MEL - Matson

    ## AUC = 0.571 [95% CI: 0.417 - 0.717 ]
    ## Testing set: MEL - McCulloch

    ## AUC = 0.591 [95% CI: 0.449 - 0.727 ]
    ## Testing set: MEL - Spencer

    ## AUC = 0.505 [95% CI: 0.45 - 0.56 ]
    ## Testing set: RARE

    ## AUC = 0.6 [95% CI: 0.514 - 0.681 ]
    ## Testing set: MEL

    ## AUC = 0.597 [95% CI: 0.549 - 0.648 ]
    ## Training set: MEL - Gopalakrishnan

    ## Testing set: RARE - GYN

    ## AUC = 0.559 [95% CI: 0.35 - 0.746 ]
    ## Testing set: RARE - NEN

    ## AUC = 0.713 [95% CI: 0.464 - 0.9 ]
    ## Testing set: RARE - UGB

    ## AUC = 0.65 [95% CI: 0.461 - 0.844 ]
    ## Testing set: MEL - Frankel

    ## AUC = 0.803 [95% CI: 0.628 - 0.936 ]
    ## Testing set: MEL - Gopalakrishnan

    ## AUC = 1 [95% CI: 1 - 1 ]
    ## Testing set: MEL - Lee

    ## AUC = 0.402 [95% CI: 0.284 - 0.523 ]
    ## Testing set: MEL - Matson

    ## AUC = 0.517 [95% CI: 0.325 - 0.707 ]
    ## Testing set: MEL - McCulloch

    ## AUC = 0.68 [95% CI: 0.532 - 0.812 ]
    ## Testing set: MEL - Spencer

    ## AUC = 0.613 [95% CI: 0.519 - 0.698 ]
    ## Testing set: RARE

    ## AUC = 0.659 [95% CI: 0.55 - 0.764 ]
    ## Testing set: MEL

    ## AUC = 0.631 [95% CI: 0.579 - 0.682 ]
    ## Training set: MEL - Lee

    ## Testing set: RARE - GYN

    ## AUC = 0.712 [95% CI: 0.519 - 0.873 ]
    ## Testing set: RARE - NEN

    ## AUC = 0.345 [95% CI: 0.142 - 0.577 ]
    ## Testing set: RARE - UGB

    ## AUC = 0.536 [95% CI: 0.346 - 0.722 ]
    ## Testing set: MEL - Frankel

    ## AUC = 0.404 [95% CI: 0.223 - 0.597 ]
    ## Testing set: MEL - Gopalakrishnan

    ## AUC = 0.442 [95% CI: 0.218 - 0.674 ]
    ## Testing set: MEL - Lee

    ## AUC = 1 [95% CI: 1 - 1 ]
    ## Testing set: MEL - Matson

    ## AUC = 0.411 [95% CI: 0.225 - 0.592 ]
    ## Testing set: MEL - McCulloch

    ## AUC = 0.494 [95% CI: 0.358 - 0.64 ]
    ## Testing set: MEL - Spencer

    ## AUC = 0.489 [95% CI: 0.408 - 0.583 ]
    ## Testing set: RARE

    ## AUC = 0.509 [95% CI: 0.397 - 0.621 ]
    ## Testing set: MEL

    ## AUC = 0.624 [95% CI: 0.567 - 0.681 ]
    ## Training set: MEL - Matson

    ## Testing set: RARE - GYN

    ## AUC = 0.67 [95% CI: 0.48 - 0.84 ]
    ## Testing set: RARE - NEN

    ## AUC = 0.672 [95% CI: 0.447 - 0.864 ]
    ## Testing set: RARE - UGB

    ## AUC = 0.455 [95% CI: 0.276 - 0.662 ]
    ## Testing set: MEL - Frankel

    ## AUC = 0.568 [95% CI: 0.378 - 0.768 ]
    ## Testing set: MEL - Gopalakrishnan

    ## AUC = 0.707 [95% CI: 0.478 - 0.889 ]
    ## Testing set: MEL - Lee

    ## AUC = 0.451 [95% CI: 0.328 - 0.573 ]
    ## Testing set: MEL - Matson

    ## AUC = 1 [95% CI: 1 - 1 ]
    ## Testing set: MEL - McCulloch

    ## AUC = 0.489 [95% CI: 0.345 - 0.638 ]
    ## Testing set: MEL - Spencer

    ## AUC = 0.548 [95% CI: 0.456 - 0.64 ]
    ## Testing set: RARE

    ## AUC = 0.593 [95% CI: 0.476 - 0.689 ]
    ## Testing set: MEL

    ## AUC = 0.587 [95% CI: 0.527 - 0.637 ]
    ## Training set: MEL - McCulloch

    ## Testing set: RARE - GYN

    ## AUC = 0.342 [95% CI: 0.176 - 0.539 ]
    ## Testing set: RARE - NEN

    ## AUC = 0.686 [95% CI: 0.453 - 0.891 ]
    ## Testing set: RARE - UGB

    ## AUC = 0.439 [95% CI: 0.244 - 0.625 ]
    ## Testing set: MEL - Frankel

    ## AUC = 0.498 [95% CI: 0.282 - 0.692 ]
    ## Testing set: MEL - Gopalakrishnan

    ## AUC = 0.533 [95% CI: 0.292 - 0.772 ]
    ## Testing set: MEL - Lee

    ## AUC = 0.505 [95% CI: 0.374 - 0.635 ]
    ## Testing set: MEL - Matson

    ## AUC = 0.414 [95% CI: 0.229 - 0.599 ]
    ## Testing set: MEL - McCulloch

    ## AUC = 1 [95% CI: 1 - 1 ]
    ## Testing set: MEL - Spencer

    ## AUC = 0.421 [95% CI: 0.326 - 0.517 ]
    ## Testing set: RARE

    ## AUC = 0.5 [95% CI: 0.387 - 0.609 ]
    ## Testing set: MEL

    ## AUC = 0.515 [95% CI: 0.459 - 0.571 ]
    ## Training set: MEL - Spencer

    ## Testing set: RARE - GYN

    ## AUC = 0.77 [95% CI: 0.588 - 0.901 ]
    ## Testing set: RARE - NEN

    ## AUC = 0.752 [95% CI: 0.519 - 0.938 ]
    ## Testing set: RARE - UGB

    ## AUC = 0.775 [95% CI: 0.585 - 0.92 ]
    ## Testing set: MEL - Frankel

    ## AUC = 0.639 [95% CI: 0.439 - 0.824 ]
    ## Testing set: MEL - Gopalakrishnan

    ## AUC = 0.793 [95% CI: 0.574 - 0.935 ]
    ## Testing set: MEL - Lee

    ## AUC = 0.541 [95% CI: 0.423 - 0.672 ]
    ## Testing set: MEL - Matson

    ## AUC = 0.459 [95% CI: 0.247 - 0.66 ]
    ## Testing set: MEL - McCulloch

    ## AUC = 0.679 [95% CI: 0.527 - 0.807 ]
    ## Testing set: MEL - Spencer

    ## AUC = 0.867 [95% CI: 0.806 - 0.912 ]
    ## Testing set: RARE

    ## AUC = 0.778 [95% CI: 0.691 - 0.859 ]
    ## Testing set: MEL

    ## AUC = 0.704 [95% CI: 0.656 - 0.755 ]
    ## Training set: RARE

    ## Testing set: RARE - GYN

    ## AUC = 1 [95% CI: 0.981 - 1 ]
    ## Testing set: RARE - NEN

    ## AUC = 1 [95% CI: 1 - 1 ]
    ## Testing set: RARE - UGB

    ## AUC = 1 [95% CI: 1 - 1 ]
    ## Testing set: MEL - Frankel

    ## AUC = 0.778 [95% CI: 0.595 - 0.918 ]
    ## Testing set: MEL - Gopalakrishnan

    ## AUC = 0.66 [95% CI: 0.419 - 0.896 ]
    ## Testing set: MEL - Lee

    ## AUC = 0.605 [95% CI: 0.47 - 0.726 ]
    ## Testing set: MEL - Matson

    ## AUC = 0.468 [95% CI: 0.269 - 0.659 ]
    ## Testing set: MEL - McCulloch

    ## AUC = 0.642 [95% CI: 0.499 - 0.764 ]
    ## Testing set: MEL - Spencer

    ## AUC = 0.612 [95% CI: 0.521 - 0.701 ]
    ## Testing set: RARE

    ## AUC = 0.999 [95% CI: 0.994 - 1 ]
    ## Testing set: MEL

    ## AUC = 0.634 [95% CI: 0.584 - 0.69 ]
    ## Training set: MEL

    ## Testing set: RARE - GYN

    ## AUC = 0.685 [95% CI: 0.481 - 0.851 ]
    ## Testing set: RARE - NEN

    ## AUC = 0.777 [95% CI: 0.596 - 0.922 ]
    ## Testing set: RARE - UGB

    ## AUC = 0.761 [95% CI: 0.571 - 0.92 ]
    ## Testing set: MEL - Frankel

    ## AUC = 0.969 [95% CI: 0.902 - 1 ]
    ## Testing set: MEL - Gopalakrishnan

    ## AUC = 0.895 [95% CI: 0.736 - 1 ]
    ## Testing set: MEL - Lee

    ## AUC = 0.904 [95% CI: 0.833 - 0.96 ]
    ## Testing set: MEL - Matson

    ## AUC = 0.812 [95% CI: 0.65 - 0.949 ]
    ## Testing set: MEL - McCulloch

    ## AUC = 0.879 [95% CI: 0.77 - 0.949 ]
    ## Testing set: MEL - Spencer

    ## AUC = 0.823 [95% CI: 0.752 - 0.882 ]
    ## Testing set: RARE

    ## AUC = 0.722 [95% CI: 0.614 - 0.818 ]
    ## Testing set: MEL

    ## AUC = 0.864 [95% CI: 0.83 - 0.9 ]

## Visualise

``` r
plot_df <- Op %>%
  mutate(
    Train_lab = ifelse(Train == "ALL", Ca_train, Train),
    Test_lab  = ifelse(Test  == "ALL", Ca_test,  Test)
  )

train_levels <- c("MEL", "Frankel", "Gopalakrishnan", "Lee", "Matson", "McCulloch", "Spencer", "RARE", "GYN", "NEN", "UGB")
test_levels  <- rev(c("MEL", "Frankel", "Gopalakrishnan", "Lee", "Matson", "McCulloch", "Spencer", "RARE", "GYN", "NEN", "UGB"))

plot_df <- plot_df %>%
  mutate(
    Train_lab = factor(Train_lab, levels = train_levels),
    Test_lab  = factor(Test_lab,  levels = test_levels)
  )

ggplot(plot_df, aes(x = Test_lab, y = Train_lab, fill = bmed)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", bmed)), size = 7) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = "Test set", y = "Train set", fill = "Median AUC") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
```

![](pancancer+melanoma_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

# Train on all except Test (leave out)

## Logic

    ## Testing set: RARE - GYN

    ## AUC = 0.696 [95% CI: 0.505 - 0.857 ]
    ## Testing set: RARE - NEN

    ## AUC = 0.71 [95% CI: 0.493 - 0.879 ]
    ## Testing set: RARE - UGB

    ## AUC = 0.697 [95% CI: 0.51 - 0.85 ]
    ## Testing set: MEL - Frankel

    ## AUC = 0.644 [95% CI: 0.44 - 0.843 ]
    ## Testing set: MEL - Gopalakrishnan

    ## AUC = 0.74 [95% CI: 0.519 - 0.917 ]
    ## Testing set: MEL - Lee

    ## AUC = 0.466 [95% CI: 0.346 - 0.588 ]
    ## Testing set: MEL - Matson

    ## AUC = 0.554 [95% CI: 0.356 - 0.741 ]
    ## Testing set: MEL - McCulloch

    ## AUC = 0.644 [95% CI: 0.496 - 0.781 ]
    ## Testing set: MEL - Spencer

    ## AUC = 0.643 [95% CI: 0.554 - 0.729 ]
    ## Testing set: RARE

    ## AUC = 0.722 [95% CI: 0.614 - 0.818 ]
    ## Testing set: MEL

    ## AUC = 0.634 [95% CI: 0.584 - 0.69 ]

## Visualise

### Boxplot

``` r
boot_df1 = boot_df1 %>%
  mutate(
    Test_lab = ifelse(Test == "ALL", Ca_test, Test),
    Test_lab = factor(Test_lab,
                      levels = unique(Test_lab))  # preserve ordering
  )

# boot_df1$Ca_test[which(boot_df1$Ca_test == "RARE" & boot_df1$Test == "ALL")] = "ALL_RARE"
# boot_df1$Ca_test[which(boot_df1$Ca_test == "MEL" & boot_df1$Test == "ALL")] = "ALL_MEL"

boot_df1$Test_lab = factor(boot_df1$Test_lab, levels = c("ALL_RARE", "ALL_MEL", "GYN", "NEN", "UGB",  "Frankel", "Gopalakrishnan", "Lee", "Matson", "McCulloch", "Spencer"))

ggplot(boot_df1, aes(x = Test_lab, y = boot, fill = Ca_test)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(
    values = c(
      "MEL"  = "#1f78b4",   # blue
      "ALL_MEL" = "#0b4f8a",
      "RARE" = "#33a02c",    # green
      "ALL_RARE" = "#1f6f1a"
    )
  ) +
  labs(x = "Test set", y = "Bootstrap AUC", fill = "Cohort") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](pancancer+melanoma_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

### Forest plot - like

``` r
Op2$Label <- ifelse(Op2$Test == "ALL", Op2$Ca_test, Op2$Test)
ggplot(Op2, aes(x = bmed, y = reorder(Label, bmed))) +
  geom_segment(aes(x = lci, xend = uci,
                   y = Label, yend = Label),
               linewidth = 1.2, color = "grey40") +
  geom_point(size = 4, color = "steelblue") +
  
  # Optional: put the numeric value next to the point
  geom_text(aes(label = sprintf("%.2f", bmed)),
            hjust = -0.3, size = 4.5) +
  
  scale_x_continuous(limits = c(min(plot_df$lci) - 0.05,
                                max(plot_df$uci) + 0.05),
                     breaks = seq(0.3, 1, by = 0.1)) +
  
  labs(x = "AUC (median, 95% bootstrap CI)",
       y = "",
       title = "") +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
```

![](pancancer+melanoma_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

## Build a Predictive Species list

### Plot this on a tree

``` r
tax = read_taxonomy("data/TAXONOMY/sylph_DB_taxonomy_99.tsv")
th2 = strainspy:::add_tax2tophits(th, tax, columns = c("Phylum","Genus","Species"))

# Fit a model using all data and take coefficients
f = 500
fit = run_enet(train_idx = 1:dim(sy)[2], test_idx = dummy_test_idx, return_fit = T, meta = meta, sy = sy, feature_names = th$Contig_name[1:f])
```

    ## Prepared data: 526 samples and 500 predictors.

    ## Prepared data: 38 samples and 500 predictors.

    ## Warning: from glmnet C++ code (error code -4); Convergence for 4th lambda value
    ## not reached after maxit=100000 iterations; solutions for larger lambdas
    ## returned

``` r
coefs <- coef(fit$finalModel, s = fit$bestTune$lambda); coefs = coefs[coefs[,1] != 0,]; coefs = coefs[-1]

# write.table(names(coefs), quote = F, row.names = F, col.names = F, "pancancer_tree_accessions.txt") 

# tree

tree = read.tree("data/ash_pancancer/tree.nwk")
tree$tip.label = unname(sapply(tree$tip.label, function(x) paste(strsplit(x, "_")[[1]][1:2], collapse = "_")))
tree = phytools::midpoint.root(tree)

df <- th2 %>%
  filter(Genome_file %in% names(coefs)) %>%
  transmute(
    tip     = Genome_file,
    Phylum  = Phylum,
    Genus   = Genus,
    Species = Species,
    Coef    = coefs[Genome_file]   # attach coefficients directly
  ) %>%
  arrange(match(tip, tree$tip.label))  # order to match tree

df = as.data.frame(df)

rownames(df) = df$tip

com_gen <- c(
  "Ruminococcus_B",
  "Faecalibacterium",
  "Phascolarctobacterium",
  "Bacteroides",
  "Actinomyces",
  "Alistipes",
  "Gemmiger",
  "Clostridium_A",
  "Collinsella",
  "Veillonella",
  "Streptococcus"
)

genus_cols <- c(
  "#000000", "#f781bf", "#ff7f00",
  "#009E73", "#984ea3", "#4daf4a",
  "#ffff33", "#a65628", "#005999", 
  "#377eb8", "#E41A1C", "#DDDDDD"
)

names(genus_cols) = c(com_gen, "Other")

df$Genus[!df$Genus %in% com_gen] = "Other"

df$Genus = factor(df$Genus, levels = c(com_gen, "Other"))


df = df[rownames(df) %in% tree$tip.label,]
tree = ape::keep.tip(tree,df$tip)

p <- ggtree(tree, layout = "circular")
p$data <- merge(
  p$data,
  df,
  by.x = "label",
  by.y = "tip",
  all.x = TRUE
)

p <- p +
  geom_tippoint(
    aes(color = Genus),
    size = 2
  )  +
  scale_color_manual(
    values = genus_cols,
    name = "Genus"
  )




df_tmp <- data.frame(Coef = df$Coef)
rownames(df_tmp) <- df$tip

p <- gheatmap(
  p,
  df_tmp,
  offset = 0,
  width = 0.05,
  colnames = FALSE
) +
  scale_fill_gradient2(
    high = "#084594",
    mid = "white",
    low = "#99000D",
    midpoint = 0,
    transform = "pseudo_log",
    name = "ENet coef"
  ) +
  theme(
    text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12)
  )
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

``` r
p
```

![](pancancer+melanoma_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

### Scatter plot with line to see the variation

``` r
df2 = df
med_df <- df2 %>%
  filter(Genus != "Other") %>%
  group_by(Genus) %>%
  summarise(
    med_coef = median(Coef, na.rm = TRUE),
    .groups = "drop"
  )


scfa_gen <- c(
  # beneficial / SCFA
  "Faecalibacterium",
  "Ruminococcus_B",
  "Phascolarctobacterium",
  "Gemmiger",
  "Eubacterium_F",
  "Bifidobacterium")

genus_order <- df2 %>%
  filter(Genus != "Other") %>%
  group_by(Genus) %>%
  summarise(
    med_coef = median(Coef, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    is_scfa = Genus %in% scfa_gen
  ) %>%
  arrange(
    desc(is_scfa),   # SCFA first
    desc(med_coef)  # then strictly decreasing coef
  ) %>%
  pull(Genus)



med_df$Genus <- factor(med_df$Genus, levels = genus_order)
df2$Genus <- factor(df2$Genus, levels = med_df$Genus[order(med_df$med_coef, decreasing = T)])

fit <- lm(med_coef ~ as.numeric(factor(Genus, levels = med_df$Genus)), data = med_df)

# Coefficients
coef(fit)
```

    ##                                      (Intercept) 
    ##                                      0.018040411 
    ## as.numeric(factor(Genus, levels = med_df$Genus)) 
    ##                                     -0.004202066

``` r
ggplot(df2, aes(x = Genus, y = Coef, color = Genus, colour = genus_cols)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.25
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_smooth(
    data = med_df,
    aes(y = med_coef, group = 1),
    method = "lm",
    se = FALSE,
    colour = "black",
    linewidth = 0.9
  ) +
  scale_color_manual(
    values = genus_cols,
    guide = "none"   # legend already in the tree
  ) +
  labs(
    x = "Ordered Genus",
    y = "Elastic-net coefficient"
  ) +
  ggthemes::theme_clean(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
```

    ## Warning: Duplicated aesthetics after name standardisation: colour

    ## `geom_smooth()` using formula = 'y ~ x'

![](pancancer+melanoma_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->
