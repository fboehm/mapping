---
title: "QTL scan with genetically additive models in CC0001 by CC042 F2 mice"
teaching: 30
exercises: 30
questions:
- "How do I perform a QTL analysis in F2 mice?"
objectives:
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
aprobs <- qtl2::genoprob_to_alleleprob(probs)
s1out_ap <- qtl2::scan1(genoprobs = aprobs, pheno = pheno[ , 1:3], kinship = kinship, addcovar = cbind(batch, pheno$sex), reml = TRUE, cores = 0)
~~~
{: .r}



~~~
qtl2::find_peaks(s1out_ap, map)
~~~
{: .r}



~~~
  lodindex lodcolumn chr   pos      lod
1        2 SpleenCFU   7 68.79 3.792543
~~~
{: .output}



~~~
qtl2::plot_scan1(s1out_ap[, 2, drop = FALSE], map = map, chr = 7)
~~~
{: .r}

<img src="../fig/rmd-14-unnamed-chunk-6-1.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" style="display: block; margin: auto;" />


~~~
qtl2::plot_scan1(s1out_ap[, 3, drop = FALSE], map = map, chr = 7)
~~~
{: .r}

<img src="../fig/rmd-14-unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" style="display: block; margin: auto;" />

