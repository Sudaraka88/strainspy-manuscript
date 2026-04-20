Simulations with spiked P. dis strain analysed using Sylph Query
================
2026-04-15

### Load data and prepare for testing

``` r
library(strainspy)
# Queries
d_q99 <- read_sylph("p_distasonis_sims/pdis_q99.tsv")
```

    ## Detected Sylph query output file.

``` r
# spiked samples
ss = readLines("p_distasonis_sims/spike_genome_ids.txt") # These are the baseline samples, not spiked
ss_spikes = c(paste(ss, "_sx1", sep = ""), paste(ss, "_sx3", sep = ""), paste(ss, "_sx5", sep = "")) # These are the spikes

# Remove _1 annoyance
colnames(d_q99) = sub("_1", "", colnames(d_q99))
d_q99@colData$Sample_file = sub("_1", "", colnames(d_q99))

# add the metadata
d_q99 = strainspy:::modify_metadata(se = d_q99, meta_data = data.frame(run_accession = colnames(d_q99),
                                                                       spiked = as.factor(sapply(colnames(d_q99), function(x) x%in% c(ss, ss_spikes)))))

# I am marking the baseline unspiked ones as spiked to show the null doesn't detect anything

d_q99 <- filter_by_presence(d_q99, min_nonzero = 25)
```

    ## Retained 16081 rows after filtering

``` r
# taxonomies
tax_99 = read_taxonomy("data/TAXONOMY/sylph_DB_taxonomy_99.tsv")

design <- as.formula(" ~ spiked")
```

## Fit ZIB model

### query 99, baseline (no signal)

``` r
# This simulation does not have any signal, we shouldn't detect anything!

# Subset
dim(d_q99)
```

    ## [1] 16081   160

``` r
rm_idx = grep("_sx", colnames(d_q99))
d_q99_test = d_q99[,-rm_idx]
dim(d_q99_test)
```

    ## [1] 16081   100

``` r
# These TRUE = 20 are NOT spiked in, this is just to show there is no detectable signal
table(d_q99_test@colData$spiked) 
```

    ## 
    ## FALSE  TRUE 
    ##    80    20

``` r
save_path = "output_rds/spike_pdis_zib_q_99_bl.rds"
if(file.exists(save_path)){
  fit_zib_99_bl <- readRDS(save_path)
} else {
  system.time({
    ebp = compute_eb_priors(d_q99_test, design, nthreads = parallel::detectCores())
    fit_zib_99_bl <- glmZiBFit(d_q99_test, design, nthreads = parallel::detectCores(), MAP_prior = ebp)
  })
  saveRDS(fit_zib_99_bl, save_path)
}

top_hits(fit_zib_99_bl) # no hits as expected
```

    ## Warning in top_hits(fit_zib_99_bl): Multiple testing correction using `holm`:
    ## No significant associations detected for coef = 2 at alpha = 0.050000

    ## # A tibble: 0 × 10
    ## # ℹ 10 variables: Contig_name <chr>, Genome_file <chr>, coefficient <dbl>,
    ## #   std_error <dbl>, p_value <dbl>, p_adjust <dbl>, zi_coefficient <dbl>,
    ## #   zi_std_error <dbl>, zi_p_value <dbl>, zi_p_adjust <dbl>

### query 99, cov1x

``` r
# baseline false samples and higher coverage ones
rm_idx = c(sapply(ss, function(x) which(colnames(d_q99) %in% x)), grep("_sx3", colnames(d_q99)), grep("_sx5", colnames(d_q99)))
d_q99_test = d_q99[,-rm_idx]
dim(d_q99_test)
```

    ## [1] 16081   100

``` r
# These 20 are spiked ata coverage of 1x
table(d_q99_test@colData$spiked) 
```

    ## 
    ## FALSE  TRUE 
    ##    80    20

``` r
save_path = "output_rds/spike_pdis_zib_q_99_cov1x.rds"
if(file.exists(save_path)){
  fit_zib_99_c1 <- readRDS(save_path)
} else {
  system.time({
    ebp = compute_eb_priors(d_q99_test, design, nthreads = parallel::detectCores())
    fit_zib_99_c1 <- glmZiBFit(d_q99_test, design, nthreads = parallel::detectCores(), MAP_prior = ebp)
  })
  saveRDS(fit_zib_99_c1, save_path)
}

th_c1 = strainspy:::add_tax2tophits(top_hits(fit_zib_99_c1), tax_99)
```

    ## Found 76 tophits for spikedTRUE at alpha = 0.05 using holm

``` r
# All 76 hits from the same species
table(th_c1$Species)
```

    ## 
    ##   Parabacteroides distasonis Parabacteroides distasonis_A 
    ##                           55                           21

``` r
plot_manhattan(fit_zib_99_c1, taxonomy = tax_99, method = "BH", tax_levels = c("Order", "Phylum", "Genus", "Species"))
```

    ## Warning: `aes_()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`
    ## ℹ The deprecated feature was likely used in the ggtree package.
    ##   Please report the issue at <https://github.com/YuLab-SMU/ggtree/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## ℹ The deprecated feature was likely used in the strainspy package.
    ##   Please report the issue to the authors.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## ℹ The deprecated feature was likely used in the ggthemes package.
    ##   Please report the issue at <https://github.com/jrnold/ggthemes/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

``` r
plot_manhattan(fit_zib_99_c1, taxonomy = tax_99, aggregate_by_taxa = F)
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-3-2.png" width="100%" />

``` r
plot_volcano(fit_zib_99_c1, label = T)
```

    ## Found 16081 tophits for spikedTRUE at alpha = 1 using holm

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-3-3.png" width="100%" />

``` r
plot_ani_dist(d_q99, phenotype = 'spiked', contigs = th_c1$Contig_name[1:5], contig_names = th_c1$Genome_file[1:5])
```

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-3-4.png" width="100%" />

``` r
# Parabacteroides distasonis is in all samples, so, this can only be a beta hit!
# GCF_024791025.1 is the top hit - it is the spiked strain. And for the spiked samples, the ANI is 100%
```

### query 99, cov3x

``` r
# baseline false samples and higher coverage ones
rm_idx = c(sapply(ss, function(x) which(colnames(d_q99) %in% x)), grep("_sx1", colnames(d_q99)), grep("_sx5", colnames(d_q99)))
d_q99_test = d_q99[,-rm_idx]
dim(d_q99_test)
```

    ## [1] 16081   100

``` r
# These 20 are spiked ata coverage of 3x
table(d_q99_test@colData$spiked) 
```

    ## 
    ## FALSE  TRUE 
    ##    80    20

``` r
save_path = "output_rds/spike_pdis_zib_q_99_cov3x.rds"
if(file.exists(save_path)){
  fit_zib_99_c3 <- readRDS(save_path)
} else {
  system.time({
    ebp = compute_eb_priors(d_q99_test, design, nthreads = parallel::detectCores())
    fit_zib_99_c3 <- glmZiBFit(d_q99_test, design, nthreads = parallel::detectCores(), MAP_prior = ebp)
  })
  saveRDS(fit_zib_99_c3, save_path)
}

th_c3 = strainspy:::add_tax2tophits(top_hits(fit_zib_99_c3), tax_99)
```

    ## Found 91 tophits for spikedTRUE at alpha = 0.05 using holm

``` r
# All 91 hits from the same species
table(th_c3$Species)
```

    ## 
    ##   Parabacteroides distasonis Parabacteroides distasonis_A 
    ##                           70                           21

``` r
plot_manhattan(fit_zib_99_c3, taxonomy = tax_99, method = "BH", tax_levels = c("Order", "Phylum", "Genus", "Species"))
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-4-1.png" width="100%" />

``` r
plot_manhattan(fit_zib_99_c3, taxonomy = tax_99, aggregate_by_taxa = F)
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?
    ## Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-4-2.png" width="100%" />

``` r
plot_volcano(fit_zib_99_c3, label = T)
```

    ## Found 16081 tophits for spikedTRUE at alpha = 1 using holm

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-4-3.png" width="100%" />

``` r
plot_ani_dist(d_q99, phenotype = 'spiked', contigs = th_c3$Contig_name[1:5], contig_names = th_c3$Genome_file[1:5])
```

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-4-4.png" width="100%" />

``` r
# Parabacteroides distasonis is in all samples, so, this can only be a beta hit!
# GCF_024791025.1 is the top hit - it is the spiked strain. And for the spiked samples, the ANI is 100%
```

### query 99, cov5x

``` r
# baseline false samples and higher coverage ones
rm_idx = c(sapply(ss, function(x) which(colnames(d_q99) %in% x)), grep("_sx1", colnames(d_q99)), grep("_sx3", colnames(d_q99)))
d_q99_test = d_q99[,-rm_idx]
dim(d_q99_test)
```

    ## [1] 16081   100

``` r
# These 20 are spiked ata coverage of 3x
table(d_q99_test@colData$spiked) 
```

    ## 
    ## FALSE  TRUE 
    ##    80    20

``` r
save_path = "output_rds/spike_pdis_zib_q_99_cov5x.rds"
if(file.exists(save_path)){
  fit_zib_99_c5 <- readRDS(save_path)
} else {
  system.time({
    ebp = compute_eb_priors(d_q99_test, design, nthreads = parallel::detectCores())
    fit_zib_99_c5 <- glmZiBFit(d_q99_test, design, nthreads = parallel::detectCores(), MAP_prior = ebp)
  })
  saveRDS(fit_zib_99_c5, save_path)
}

th_c5 = strainspy:::add_tax2tophits(top_hits(fit_zib_99_c5), tax_99)
```

    ## Found 78 tophits for spikedTRUE at alpha = 0.05 using holm

``` r
# All 91 hits from the same species
table(th_c5$Species)
```

    ## 
    ##   Parabacteroides distasonis Parabacteroides distasonis_A 
    ##                           57                           21

``` r
plot_manhattan(fit_zib_99_c5, taxonomy = tax_99, method = "BH", tax_levels = c("Order", "Phylum", "Genus", "Species"))
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-5-1.png" width="100%" />

``` r
plot_manhattan(fit_zib_99_c5, taxonomy = tax_99, aggregate_by_taxa = F)
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?
    ## Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    ## ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • branch.length = branch.length
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-5-2.png" width="100%" />

``` r
plot_volcano(fit_zib_99_c5, label = T)
```

    ## Found 16081 tophits for spikedTRUE at alpha = 1 using holm

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-5-3.png" width="100%" />

``` r
plot_ani_dist(d_q99, phenotype = 'spiked', contigs = th_c5$Contig_name[1:5], contig_names = th_c5$Genome_file[1:5])
```

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-5-4.png" width="100%" />

``` r
# Parabacteroides distasonis is in all samples, so, this can only be a beta hit!
# GCF_024791025.1 is the top hit - it is the spiked strain. And for the spiked samples, the ANI is 100%
```

StrainSpy works regardless of coverage, now to test AnPan

### Load data and prepare for testing

``` r
library(anpan)
```

    ## • This is anpan version 0.3.0
    ## • Read the guide: run anpan::anpan_vignette()
    ## • Get help: Visit the biobakery help forum at <https://forum.biobakery.org/>
    ## • Parallelize: Run `future::plan()` as appropriate for your system.
    ## • Activate progress bars: `library(progressr); handlers(global=TRUE)`

``` r
library(ape)
library(tibble)
```

    ## Warning: package 'tibble' was built under R version 4.4.3

``` r
tree = read.tree("p_distasonis_sims/spiked_tree.tre")

metadata = tibble(sample_id = tree$tip.label,
                  spiked   = as.factor(sapply(tree$tip.label, function(x) x%in% c(ss, sub("_", "_1_", paste(ss_spikes, ".fastq", sep = "")) ))) )

table(metadata$spiked)
```

    ## 
    ## FALSE  TRUE 
    ##    64    72

``` r
plot_outcome_tree(tree,
                  metadata, 
                  covariates = c(),
                  outcome    = "spiked")
```

    ## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`.
    ## ℹ See also `vignette("ggplot2-in-packages")` for more information.
    ## ℹ The deprecated feature was likely used in the anpan package.
    ##   Please report the issue to the authors.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-6-1.png" width="100%" />

Looks like the left-most branch has a bunch of spiked ones, it’ll be
interesting to see variation with coverage

## anPan PGLMM

### baseline (no signal)

``` r
# This simulation does not have any signal, we shouldn't detect anything!

# Subset
rm_idx = grep("_sx", metadata$sample_id)

metadata_test = metadata[-rm_idx, ]
table(metadata_test$spiked)
```

    ## 
    ## FALSE  TRUE 
    ##    64    12

``` r
# Only 76/100 samples had P. dist at a detectable level for strainphlan
# These TRUE = 12 are NOT spiked in, this is just to show there is no detectable signal

# Drop the same leaves from the tree
tree_test = drop.tip(tree, metadata$sample_id[rm_idx])
plot_outcome_tree(tree_test,
                  metadata_test, 
                  covariates = c(),
                  outcome    = "spiked")
```

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-7-1.png" width="100%" />

``` r
# Clearly no signal

save_path = "output_rds/spike_pdis_anpan_bl.rds"
if(file.exists(save_path)){
  anpan_res_bl <- readRDS(save_path)
} else {
  system.time({
    anpan_res_bl = anpan_pglmm(meta_file       = metadata_test,
                               tree_file       = tree_test,
                               outcome         = "spiked",
                               family          = "binomial",
                               bug_name        = "Pdis",
                               reg_noise       = TRUE,
                               loo_comparison  = TRUE,
                               run_diagnostics = FALSE,
                               refresh         = 500,
                               show_plot_tree  = FALSE,
                               show_post       = FALSE)
  })
  saveRDS(anpan_res_bl, save_path)
}

anpan_res_bl$loo$comparison
```

    ##             elpd_diff   se_diff  elpd_loo se_elpd_loo    p_loo  se_p_loo
    ## base_fit   0.00000000 0.0000000 -34.10157    5.094281 0.855544 0.1756572
    ## pglmm_fit -0.09731679 0.1441665 -34.19889    5.065588 1.205773 0.2441358
    ##              looic se_looic
    ## base_fit  68.20314 10.18856
    ## pglmm_fit 68.39777 10.13118
    ## attr(,"class")
    ## [1] "compare.loo" "matrix"      "array"

Clearly no signal - confirmed by PGLMM

### cov1x

``` r
# Subset
rm_idx = unlist(c(sapply(ss, function(x) which(metadata$sample_id %in% x)), grep("_sx3", metadata$sample_id), grep("_sx5", metadata$sample_id)))

metadata_test = metadata[-rm_idx, ]
table(metadata_test$spiked)
```

    ## 
    ## FALSE  TRUE 
    ##    64    20

``` r
# Only 20/20 spiked samples are now detected, but only non spiked are present - should be plenty 

# Drop the same leaves from the tree
tree_test = drop.tip(tree, metadata$sample_id[rm_idx])
plot_outcome_tree(tree_test,
                  metadata_test, 
                  covariates = c(),
                  outcome    = "spiked")
```

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-8-1.png" width="100%" />

``` r
# Signal is a bit weak

save_path = "output_rds/spike_pdis_anpan_cov1x.rds"
if(file.exists(save_path)){
  anpan_res_c1 <- readRDS(save_path)
} else {
  system.time({
    anpan_res_c1 = anpan_pglmm(meta_file       = metadata_test,
                               tree_file       = tree_test,
                               outcome         = "spiked",
                               family          = "binomial",
                               bug_name        = "Pdis",
                               reg_noise       = TRUE,
                               loo_comparison  = TRUE,
                               run_diagnostics = FALSE,
                               refresh         = 500,
                               show_plot_tree  = FALSE,
                               show_post       = FALSE)
  })
  saveRDS(anpan_res_c1, save_path)
}

anpan_res_c1$loo$comparison
```

    ##           elpd_diff   se_diff  elpd_loo se_elpd_loo     p_loo  se_p_loo
    ## pglmm_fit  0.000000 0.0000000 -45.52922    4.278753 2.1446260 0.2959159
    ## base_fit  -1.612409 0.8160137 -47.14163    4.425549 0.9993287 0.1292906
    ##              looic se_looic
    ## pglmm_fit 91.05844 8.557506
    ## base_fit  94.28326 8.851099
    ## attr(,"class")
    ## [1] "compare.loo" "matrix"      "array"

The phylogenetic model seems to fit better, but the difference doesn’t
seem clear (less than 2 standard errors difference in ELPD).

### cov3x

``` r
# Subset
rm_idx = unlist(c(sapply(ss, function(x) which(metadata$sample_id %in% x)), grep("_sx1", metadata$sample_id), grep("_sx5", metadata$sample_id)))

metadata_test = metadata[-rm_idx, ]
table(metadata_test$spiked)
```

    ## 
    ## FALSE  TRUE 
    ##    64    20

``` r
# Only 20/20 spiked samples are now detected, but only non spiked are present - should be plenty 

# Drop the same leaves from the tree
tree_test = drop.tip(tree, metadata$sample_id[rm_idx])
plot_outcome_tree(tree_test,
                  metadata_test, 
                  covariates = c(),
                  outcome    = "spiked")
```

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-9-1.png" width="100%" />

``` r
# Clearly signal is present

save_path = "output_rds/spike_pdis_anpan_cov3x.rds"
if(file.exists(save_path)){
  anpan_res_c3 <- readRDS(save_path)
} else {
  system.time({
    anpan_res_c3 = anpan_pglmm(meta_file       = metadata_test,
                               tree_file       = tree_test,
                               outcome         = "spiked",
                               family          = "binomial",
                               bug_name        = "Pdis",
                               reg_noise       = TRUE,
                               loo_comparison  = TRUE,
                               run_diagnostics = FALSE,
                               refresh         = 500,
                               show_plot_tree  = FALSE,
                               show_post       = FALSE)
  })
  saveRDS(anpan_res_c3, save_path)
}

anpan_res_c3$loo$comparison
```

    ##           elpd_diff  se_diff  elpd_loo se_elpd_loo    p_loo  se_p_loo    looic
    ## pglmm_fit  0.000000 0.000000 -38.39092    3.918328 2.993786 0.4993756 76.78183
    ## base_fit  -8.657758 2.351852 -47.04867    4.399150 0.902703 0.1160073 94.09735
    ##           se_looic
    ## pglmm_fit 7.836656
    ## base_fit  8.798300
    ## attr(,"class")
    ## [1] "compare.loo" "matrix"      "array"

The phylogenetic model seems to fit better, and the difference seems
clear (more than 2 standard errors difference in ELPD).

### cov5x

``` r
# Subset
rm_idx = unlist(c(sapply(ss, function(x) which(metadata$sample_id %in% x)), grep("_sx1", metadata$sample_id), grep("_sx3", metadata$sample_id)))

metadata_test = metadata[-rm_idx, ]
table(metadata_test$spiked)
```

    ## 
    ## FALSE  TRUE 
    ##    64    20

``` r
# Only 20/20 spiked samples are now detected, but only non spiked are present - should be plenty 

# Drop the same leaves from the tree
tree_test = drop.tip(tree, metadata$sample_id[rm_idx])
plot_outcome_tree(tree_test,
                  metadata_test, 
                  covariates = c(),
                  outcome    = "spiked")
```

<img src="sims_spike_query_pdis_files/figure-gfm/unnamed-chunk-10-1.png" width="100%" />

``` r
# Clearly dominant signal

save_path = "output_rds/spike_pdis_anpan_cov5x.rds"
if(file.exists(save_path)){
  anpan_res_c5 <- readRDS(save_path)
} else {
  system.time({
    anpan_res_c5 = anpan_pglmm(meta_file       = metadata_test,
                               tree_file       = tree_test,
                               outcome         = "spiked",
                               family          = "binomial",
                               bug_name        = "Pdis",
                               reg_noise       = TRUE,
                               loo_comparison  = TRUE,
                               run_diagnostics = FALSE,
                               refresh         = 500,
                               show_plot_tree  = FALSE,
                               show_post       = FALSE)
  })
  saveRDS(anpan_res_c5, save_path)
}

anpan_res_c5$loo$comparison
```

    ##           elpd_diff  se_diff  elpd_loo se_elpd_loo    p_loo  se_p_loo    looic
    ## pglmm_fit   0.00000 0.000000 -36.09019    3.860277 2.839414 0.5063496 72.18038
    ## base_fit  -11.05142 2.735972 -47.14161    4.429170 1.000522 0.1284311 94.28322
    ##           se_looic
    ## pglmm_fit 7.720554
    ## base_fit  8.858339
    ## attr(,"class")
    ## [1] "compare.loo" "matrix"      "array"

The phylogenetic model seems to fit better, and the difference seems
clear (more than 2 standard errors difference in ELPD).

AnPan works as well - when there is signal in the phylogeny, it can
detect it accurately. However, building the tree can be tricky if the
strain coverage is low. StrainSpy can work equally well in the lower
coverage simulations as well, without requiring a phylogeny.
