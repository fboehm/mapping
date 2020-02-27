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
map <- readRDS("../data/derived_data/map.rds")
~~~
{: .r}



~~~
Warning in gzfile(file, "rb"): cannot open compressed file '../data/
derived_data/map.rds', probable reason 'No such file or directory'
~~~
{: .error}



~~~
Error in gzfile(file, "rb"): cannot open the connection
~~~
{: .error}



~~~
probs <- readRDS("../data/derived_data/probs.rds")
~~~
{: .r}



~~~
Warning in gzfile(file, "rb"): cannot open compressed file '../data/
derived_data/probs.rds', probable reason 'No such file or directory'
~~~
{: .error}



~~~
Error in gzfile(file, "rb"): cannot open the connection
~~~
{: .error}



~~~
pheno <- readRDS("../data/derived_data/pheno.rds")
~~~
{: .r}



~~~
Warning in gzfile(file, "rb"): cannot open compressed file '../data/
derived_data/pheno.rds', probable reason 'No such file or directory'
~~~
{: .error}



~~~
Error in gzfile(file, "rb"): cannot open the connection
~~~
{: .error}

We tweak the rownames so that subjects' phenotype values can be matched with their genotypes. We also create a set of vectors to indicate batch number.


~~~
rownames(pheno) <- rownames(probs$`1`)
~~~
{: .r}



~~~
Error in rownames(probs$`1`): object 'probs' not found
~~~
{: .error}



~~~
batch <- tibble::tibble(batch1 = pheno$batch == 1, batch2 = pheno$batch == 2, batch3 = pheno$batch == 3) %>%
  purrr::map_df(.f = as.numeric) %>%
  as.matrix %>%
  (function(x) {rownames(x) <- rownames(probs$`1`); return(x)})
~~~
{: .r}



~~~
Error in eval_tidy(xs[[i]], unique_output): object 'pheno' not found
~~~
{: .error}

We want to incorporate polygenic effects into our models, so we calculate a relatedness matrix. 


~~~
kinship <- qtl2::calc_kinship(probs)
~~~
{: .r}



~~~
Error in qtl2::calc_kinship(probs): object 'probs' not found
~~~
{: .error}

Now, we can perform the QTL scans.



~~~
aprobs <- qtl2::genoprob_to_alleleprob(probs)
~~~
{: .r}



~~~
Error in qtl2::genoprob_to_alleleprob(probs): object 'probs' not found
~~~
{: .error}



~~~
s1out_ap <- qtl2::scan1(genoprobs = aprobs, pheno = pheno[ , 1:3], kinship = kinship, addcovar = cbind(batch, pheno$sex), reml = TRUE, cores = 0)
~~~
{: .r}



~~~
Error in qtl2::scan1(genoprobs = aprobs, pheno = pheno[, 1:3], kinship = kinship, : object 'aprobs' not found
~~~
{: .error}



~~~
qtl2::find_peaks(s1out_ap, map)
~~~
{: .r}



~~~
Error in is.data.frame(map): object 'map' not found
~~~
{: .error}



~~~
qtl2::plot_scan1(s1out_ap[, 2, drop = FALSE], map = map, chr = 7)
~~~
{: .r}



~~~
Error in qtl2::plot_scan1(s1out_ap[, 2, drop = FALSE], map = map, chr = 7): object 'map' not found
~~~
{: .error}


~~~
qtl2::plot_scan1(s1out_ap[, 3, drop = FALSE], map = map, chr = 7)
~~~
{: .r}



~~~
Error in qtl2::plot_scan1(s1out_ap[, 3, drop = FALSE], map = map, chr = 7): object 'map' not found
~~~
{: .error}

