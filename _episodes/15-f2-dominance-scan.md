---
title: "QTL scan with genetically dominant models in CC0001 by CC042 F2 mice"
teaching: 30
exercises: 30
questions:
- "How do I encode genotypes to reflect dominance models?"
objectives:
- Recode genotypes to reflect dominance for one allele or the other.
- Perform a basic QTL analysis with those recoded genotypes.
- Identify QTL with a genome scan.
keypoints:
- "."
- "."
source: Rmd
---





First, we read in the saved data files. These contain the genome map, genotype probabilities, and trait values.


~~~
library(magrittr)
map <- readRDS("../data/derived_data/map.rds")
probs <- readRDS("../data/derived_data/probs.rds")
pheno <- readRDS("../data/derived_data/pheno.rds")
~~~
{: .r}

We tweak the rownames so that subjects' phenotype values can be matched with their genotypes. We also create a set of vectors to indicate batch number.


~~~
rownames(pheno) <- rownames(probs$`1`)
batch <- tibble::tibble(batch1 = pheno$batch == 1, batch2 = pheno$batch == 2, batch3 = pheno$batch == 3) %>%
  purrr::map_df(.f = as.numeric) %>%
  as.matrix %>%
  (function(x) {rownames(x) <- rownames(probs$`1`); return(x)})
~~~
{: .r}

We want to incorporate polygenic effects into our models, so we calculate a relatedness matrix. 


~~~
kinship <- qtl2::calc_kinship(probs)
~~~
{: .r}

Now, we can perform the QTL scans.



~~~
aprobs <- qtl2::genoprob_to_alleleprob(probs)
s1out_ap <- qtl2::scan1(genoprobs = aprobs, pheno = pheno[ , 1:3], kinship = kinship, addcovar = cbind(batch, pheno$sex), reml = TRUE, cores = 0)
~~~
{: .r}



## Function to convert genotype probabilities to binary

Our goal is to recode the genotypes according to dominance models and then do univariate QTL scans with those recoded genotypes.

So, we'll do a single QTL scan for each of the two encodings:

1. AA -> 0; AB -> 1; BB -> 1 (B dominant)

2. AA -> 1; AB -> 1; BB -> 0 (A dominant)

So, inputs to `qtl2::scan1` will be a genoprobs object that is a list of arrays, each n by 1 by n.markers. The one column for each array is the probability of having the dominant allele.


~~~
length(probs)
~~~
{: .r}



~~~
[1] 20
~~~
{: .output}



~~~
colnames(probs$`1`)
~~~
{: .r}



~~~
[1] "AA" "AB" "BB"
~~~
{: .output}



~~~
# A dominant
probs$`1` -> pp
sum_columns <- function(x, dominant_allele = "A"){
  if (dominant_allele == "A") {da_numeric <- 1} else {da_numeric <- 3}
  out <- x[, da_numeric, ] + x[, 2, ]
  return(out)
}
#sum_columns(pp) %>% array(dim = c(nrow(.), 1, ncol(.)))

probs_A_dominant <- probs %>%
  purrr::map(.f = function(x){
    sum_columns(x) -> ss
  }) %>%
  purrr::map(.f = function(x){array(data = x, dim = c(nrow(x), 1, ncol(x)))}) %>%
  purrr::map(.f = function(x){abind::abind(x, 1-x, along = 2)}) %>%
  purrr::map2(.y = probs, .f = function(x, y) {dimnames(x)[[1]] <- dimnames(y)[[1]]; dimnames(x)[[2]] <- c("dominant_allele", "recessive_allele"); dimnames(x)[[3]] <- dimnames(y)[[3]]; return(x)}) %>%
  (function(x){attributes(x) <- attributes(probs); return(x)})
~~~
{: .r}



~~~
Error in loadNamespace(name): there is no package called 'abind'
~~~
{: .error}



~~~
probs_B_dominant <- probs %>%
  purrr::map(.f = sum_columns, dominant_allele = "B") %>%
  purrr::map(.f = function(x){array(data = x, dim = c(nrow(x), 1, ncol(x)))}) %>%
  purrr::map(.f = function(x){abind::abind(x, 1-x, along = 2)}) %>%
  purrr::map2(.y = probs, .f = function(x, y) {dimnames(x)[[1]] <- dimnames(y)[[1]]; dimnames(x)[[2]] <- c("dominant_allele", "recessive_allele"); dimnames(x)[[3]] <- dimnames(y)[[3]]; return(x)}) %>%
  (function(x){attributes(x) <- attributes(probs); return(x)})
~~~
{: .r}



~~~
Error in loadNamespace(name): there is no package called 'abind'
~~~
{: .error}


~~~
A_dom<- qtl2::scan1(genoprobs = probs_A_dominant, pheno = pheno[, 1:3], addcovar = cbind(batch, pheno$sex), kinship = kinship, reml = TRUE, cores = 0)
~~~
{: .r}



~~~
Error in qtl2::scan1(genoprobs = probs_A_dominant, pheno = pheno[, 1:3], : object 'probs_A_dominant' not found
~~~
{: .error}


~~~
B_dom<- qtl2::scan1(genoprobs = probs_B_dominant, pheno = pheno[, 1:3], addcovar = cbind(batch, pheno$sex), kinship = kinship, reml = TRUE, cores = 0)
~~~
{: .r}



~~~
Error in qtl2::scan1(genoprobs = probs_B_dominant, pheno = pheno[, 1:3], : object 'probs_B_dominant' not found
~~~
{: .error}


~~~
qtl2::find_peaks(A_dom, map, threshold = 2)
~~~
{: .r}



~~~
Error in align_scan1_map(scan1_output, map): object 'A_dom' not found
~~~
{: .error}


~~~
qtl2::find_peaks(B_dom, map, threshold = 2)
~~~
{: .r}



~~~
Error in align_scan1_map(scan1_output, map): object 'B_dom' not found
~~~
{: .error}


~~~
qtl2::plot_scan1(A_dom[, 2, drop = FALSE], map$`7`)
~~~
{: .r}



~~~
Error in align_scan1_map(x, map): object 'A_dom' not found
~~~
{: .error}



~~~
qtl2::plot_scan1(A_dom[, 3, drop = FALSE], map$`7`)
~~~
{: .r}



~~~
Error in align_scan1_map(x, map): object 'A_dom' not found
~~~
{: .error}


~~~
qtl2::plot_scan1(B_dom[, 2, drop = FALSE], map$`7`)
~~~
{: .r}



~~~
Error in align_scan1_map(x, map): object 'B_dom' not found
~~~
{: .error}



~~~
qtl2::plot_scan1(B_dom[, 3, drop = FALSE], map$`7`)
~~~
{: .r}



~~~
Error in align_scan1_map(x, map): object 'B_dom' not found
~~~
{: .error}




## Next step: get permutation pvalues with 1000 perms.


~~~
A_dom_perms<- qtl2::scan1perm(genoprobs = probs_A_dominant, pheno = pheno[, 1:3], addcovar = cbind(batch, pheno$sex), kinship = kinship, reml = TRUE, cores = 0, n_perm = 1000)
saveRDS(A_dom_perms, "../data/derived_data/perms_A_dominance_F2.rds")
~~~
{: .r}


~~~
B_dom_perms <- qtl2::scan1perm(genoprobs = probs_B_dominant, pheno = pheno[, 1:3], addcovar = cbind(batch, pheno$sex), kinship = kinship, reml = TRUE, cores = 0, n_perm = 1000)
saveRDS(B_dom_perms, "../data/derived_data/perms_B_dominance_F2.rds")
~~~
{: .r}


~~~
A_dom_perms <- readRDS("../data/derived_data/perms_A_dominance_F2.rds")
~~~
{: .r}



~~~
Warning in gzfile(file, "rb"): cannot open compressed file '../data/
derived_data/perms_A_dominance_F2.rds', probable reason 'No such file or
directory'
~~~
{: .error}



~~~
Error in gzfile(file, "rb"): cannot open the connection
~~~
{: .error}



~~~
B_dom_perms <- readRDS("../data/derived_data/perms_B_dominance_F2.rds")
~~~
{: .r}



~~~
Warning in gzfile(file, "rb"): cannot open compressed file '../data/
derived_data/perms_B_dominance_F2.rds', probable reason 'No such file or
directory'
~~~
{: .error}



~~~
Error in gzfile(file, "rb"): cannot open the connection
~~~
{: .error}


~~~
get_pvalue <- function(observed_lod, perm_lods){
  mean(perm_lods >= observed_lod)
}

qtl2::find_peaks(A_dom, map) %>%
  dplyr::mutate(pvalue = purrr::map2_dbl(.x = lod, .y = lodindex, .f = function(x, y) get_pvalue(x, A_dom_perms[ , y])))
~~~
{: .r}



~~~
Error in align_scan1_map(scan1_output, map): object 'A_dom' not found
~~~
{: .error}


~~~
# get quantiles of perms
quantile(A_dom_perms[, 1], c(0.95, 0.9))
~~~
{: .r}



~~~
Error in quantile(A_dom_perms[, 1], c(0.95, 0.9)): object 'A_dom_perms' not found
~~~
{: .error}



~~~
quantile(A_dom_perms[, 2], c(0.95, 0.9))
~~~
{: .r}



~~~
Error in quantile(A_dom_perms[, 2], c(0.95, 0.9)): object 'A_dom_perms' not found
~~~
{: .error}



~~~
quantile(A_dom_perms[, 3], c(0.95, 0.9))
~~~
{: .r}



~~~
Error in quantile(A_dom_perms[, 3], c(0.95, 0.9)): object 'A_dom_perms' not found
~~~
{: .error}



~~~
quantile(B_dom_perms[, 1], c(0.95, 0.9))
~~~
{: .r}



~~~
Error in quantile(B_dom_perms[, 1], c(0.95, 0.9)): object 'B_dom_perms' not found
~~~
{: .error}



~~~
quantile(B_dom_perms[, 2], c(0.95, 0.9))
~~~
{: .r}



~~~
Error in quantile(B_dom_perms[, 2], c(0.95, 0.9)): object 'B_dom_perms' not found
~~~
{: .error}



~~~
quantile(B_dom_perms[, 3], c(0.95, 0.9))
~~~
{: .r}



~~~
Error in quantile(B_dom_perms[, 3], c(0.95, 0.9)): object 'B_dom_perms' not found
~~~
{: .error}

## Analysis for hets v. non-hets


~~~
pr_hets <- probs %>%
  purrr::map(.f = function(x){
    x[, 1 , ] <- x[ , 3, ] + x[ , 1, ]
    return(x[ , 1:2, ])
  }) %>%
  (function(x){attributes(x) <- attributes(probs); return(x)})
s1out_het <- qtl2::scan1(genoprobs = pr_hets, 
                         pheno = pheno[, 1:3], 
                         addcovar = cbind(batch, pheno$sex), 
                         kinship = kinship, 
                         reml = TRUE, 
                         cores = 1
                         )
qtl2::plot_scan1(s1out_het, map = map)
~~~
{: .r}

<img src="../fig/rmd-15-unnamed-chunk-17-1.png" title="plot of chunk unnamed-chunk-17" alt="plot of chunk unnamed-chunk-17" style="display: block; margin: auto;" />



~~~
het_perms <- qtl2::scan1perm(genoprobs = pr_hets, 
                             pheno = pheno[, 1:3], 
                             addcovar = cbind(batch, pheno$sex), 
                             kinship = kinship, 
                             reml = TRUE, 
                             cores = parallel::detectCores() - 1, 
                             n_perm = 1000
                             )
saveRDS(het_perms, "../data/derived_data/het_perms.rds")
~~~
{: .r}


~~~
het_perms <- readRDS("../data/derived_data/het_perms.rds")
~~~
{: .r}



~~~
Warning in gzfile(file, "rb"): cannot open compressed file '../data/
derived_data/het_perms.rds', probable reason 'No such file or directory'
~~~
{: .error}



~~~
Error in gzfile(file, "rb"): cannot open the connection
~~~
{: .error}



~~~
summary(het_perms)
~~~
{: .r}



~~~
Error in summary(het_perms): object 'het_perms' not found
~~~
{: .error}


~~~
s1out_het %>%
  qtl2::find_peaks(map = map) %>%
  dplyr::mutate(pvalue = purrr::map2_dbl(.x = lod, .y = lodindex, .f = function(x, y) get_pvalue(x, het_perms[ , y])))
~~~
{: .r}



~~~
Error in mean(perm_lods >= observed_lod): object 'het_perms' not found
~~~
{: .error}




