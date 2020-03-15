---
title: "QTL scan with genetically dominant models in CC0001 by CC042 F2 mice"
teaching: 30
exercises: 30
questions:
- "How do I use 'contrasts' to reflect dominance models?"
objectives:
- Use contrasts with 'contrasts' argument to `qtl2::fit1`.
- Perform a basic QTL analysis.
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
out <- list()
for (j in 1:18){
  pr <- probs[[j]]
  f1out <- list()
  for (i in 1:dim(pr)[[3]]){
    f1out[[i]] <- qtl2::fit1(genoprobs = pr[ , , i], 
                         pheno = pheno[ , 1, drop = FALSE], #lung cfu trait
                         kinship = kinship, 
                         addcovar = cbind(batch, pheno$sex), 
                         contrasts = cbind(mu=c(1, 1, 1), 
                                           a=c(-1, 0, 1), 
                                           d=c(0, 1, 0)), 
                         reml = TRUE, 
                         cores = 0
                         )

  }
  out[[j]] <- f1out
}
saveRDS(out, "../data/out-fit1.rds")
~~~
{: .r}




~~~
s1out <- qtl2::scan1(genoprobs = probs, 
                     pheno = pheno[, 1:3], 
                     kinship = kinship, 
                     addcovar = cbind(batch, pheno$sex), 
                     cores = 0
                     )
qtl2::find_peaks(s1out, map = map, threshold = 2)
~~~
{: .r}



~~~
   lodindex lodcolumn chr    pos      lod
1         1   LungCFU   3 49.523 2.060117
2         1   LungCFU   6 78.325 2.353980
3         1   LungCFU   7 35.675 2.025513
4         1   LungCFU  14 65.779 2.378495
5         1   LungCFU   X 38.355 2.177884
6         2 SpleenCFU   7 68.790 6.849362
7         2 SpleenCFU   8 25.667 2.210979
8         3      IFNg   7 67.893 6.663188
9         3      IFNg   9  2.476 2.207654
10        3      IFNg  10 62.215 2.155271
11        3      IFNg  15 38.599 2.997543
~~~
{: .output}


~~~
aprobs <- qtl2::genoprob_to_alleleprob(probs)
s1outa <- qtl2::scan1(genoprobs = aprobs, 
                     pheno = pheno[, 1:3], 
                     kinship = kinship, 
                     addcovar = cbind(batch, pheno$sex), 
                     cores = 0
                     )
qtl2::find_peaks(s1outa, map = map, threshold = 2)
~~~
{: .r}



~~~
  lodindex lodcolumn chr    pos      lod
1        1   LungCFU   7 35.675 2.025201
2        2 SpleenCFU   7 68.790 3.792543
3        3      IFNg   7 68.401 2.974409
4        3      IFNg  15 38.849 2.729513
~~~
{: .output}


## Dominance effects for lung CFU

We'll extract lods from the object `out`, then look at the differences between these LODs and the LODs of the `scan1` output, `scan1out`.


~~~
lod_a_plus_d <- out %>%
  purrr::map(.f = function(x)purrr::map_dbl(.x = x, .f = function(x)x$lod)) %>%
  unlist()
~~~
{: .r}


~~~
domlod <- s1out %>%
  tibble::as_tibble() %>%
  dplyr::mutate(markernum = 1:7467) %>%
  dplyr::filter(markernum <= 6880) %>%
  dplyr::mutate(lungCFU_a_plus_d = lod_a_plus_d, 
                lod_dominance = lungCFU_a_plus_d - LungCFU) %>%
  dplyr::select(lod_dominance) %>%
  as.matrix() %>%
  (function(x){rownames(x)<- rownames(s1out)[1:6880]; return(x)}) 

qtl2::plot_scan1(domlod, map)
~~~
{: .r}

<img src="../fig/rmd-16-unnamed-chunk-9-1.png" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" style="display: block; margin: auto;" />

~~~
qtl2::find_peaks(domlod, map, threshold = 1, peakdrop = 1, drop = 0.5)
~~~
{: .r}



~~~
[1] lodindex  lodcolumn chr       pos       lod       ci_lo     ci_hi    
<0 rows> (or 0-length row.names)
~~~
{: .output}

## Dominance effects for spleen CFU


~~~
sp_out <- list()
for (j in 1:18){
  pr <- probs[[j]]
  f1out <- list()
  for (i in 1:dim(pr)[[3]]){
    f1out[[i]] <- qtl2::fit1(genoprobs = pr[ , , i], 
                         pheno = pheno[ , 2, drop = FALSE], #lung cfu trait
                         kinship = kinship, 
                         addcovar = cbind(batch, pheno$sex), 
                         contrasts = cbind(mu=c(1, 1, 1), 
                                           a=c(-1, 0, 1), 
                                           d=c(0, 1, 0)), 
                         reml = TRUE, 
                         cores = 0
                         )

  }
  sp_out[[j]] <- f1out
}
saveRDS(sp_out, "../data/out-fit1-spleenCFU.rds")
~~~
{: .r}


~~~
lod_a_plus_d <- sp_out %>%
  purrr::map(.f = function(x)purrr::map_dbl(.x = x, .f = function(x)x$lod)) %>%
  unlist()
~~~
{: .r}


~~~
domlod <- s1out %>%
  tibble::as_tibble() %>%
  dplyr::mutate(markernum = 1:7467) %>%
  dplyr::filter(markernum <= 6880) %>%
  dplyr::mutate(spleenCFU_a_plus_d = lod_a_plus_d, 
                lod_dominance = spleenCFU_a_plus_d - SpleenCFU) %>%
  dplyr::select(lod_dominance) %>%
  as.matrix() %>%
  (function(x){rownames(x)<- rownames(s1out)[1:6880]; return(x)})
qtl2::plot_scan1(domlod, map)                   
~~~
{: .r}

<img src="../fig/rmd-16-unnamed-chunk-12-1.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" style="display: block; margin: auto;" />

~~~
qtl2::find_peaks(domlod, map, threshold = 1, peakdrop = 1, drop = 0.5)
~~~
{: .r}



~~~
[1] lodindex  lodcolumn chr       pos       lod       ci_lo     ci_hi    
<0 rows> (or 0-length row.names)
~~~
{: .output}




## Dominance effects for IFNgamma


~~~
ifn_out <- list()
for (j in 1:18){
  pr <- probs[[j]]
  f1out <- list()
  for (i in 1:dim(pr)[[3]]){
    f1out[[i]] <- qtl2::fit1(genoprobs = pr[ , , i], 
                         pheno = pheno[ , 3, drop = FALSE], #lung cfu trait
                         kinship = kinship, 
                         addcovar = cbind(batch, pheno$sex), 
                         contrasts = cbind(mu=c(1, 1, 1), 
                                           a=c(-1, 0, 1), 
                                           d=c(0, 1, 0)), 
                         reml = TRUE, 
                         cores = 0
                         )

  }
  ifn_out[[j]] <- f1out
}
saveRDS(ifn_out, "../data/out-fit1-ifn.rds")
~~~
{: .r}



~~~
lod_a_plus_d <- ifn_out %>%
  purrr::map(.f = function(x)purrr::map_dbl(.x = x, .f = function(x)x$lod)) %>%
  unlist()
~~~
{: .r}


~~~
domlod <- s1out %>%
  tibble::as_tibble() %>%
  dplyr::mutate(markernum = 1:7467) %>%
  dplyr::filter(markernum <= 6880) %>%
  dplyr::mutate(ifng_a_plus_d = lod_a_plus_d, 
                lod_dominance = ifng_a_plus_d - IFNg) %>%
  dplyr::select(lod_dominance) %>%
  as.matrix() %>%
  (function(x){rownames(x)<- rownames(s1out)[1:6880]; return(x)})
qtl2::plot_scan1(domlod, map)                   
~~~
{: .r}

<img src="../fig/rmd-16-unnamed-chunk-15-1.png" title="plot of chunk unnamed-chunk-15" alt="plot of chunk unnamed-chunk-15" style="display: block; margin: auto;" />

~~~
qtl2::find_peaks(domlod, map, threshold = 1, peakdrop = 1, drop = 0.5)
~~~
{: .r}



~~~
[1] lodindex  lodcolumn chr       pos       lod       ci_lo     ci_hi    
<0 rows> (or 0-length row.names)
~~~
{: .output}

