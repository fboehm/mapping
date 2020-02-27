---
layout: page
title: Setup
permalink: /setup/
---
## Installation

R is a programming language that is especially powerful for data exploration, visualization, and statistical analysis. To interact with R, we use RStudio. 

1. Install the latest version of R from [CRAN](https://cran.r-project.org/).

2. Install the latest version of RStudio [here](https://www.rstudio.com/products/rstudio/download/). Choose the free RStudio Desktop version for Windows, Mac, or Linux. 

3. Start RStudio. The [qtl2](https://github.com/rqtl/qtl2) package contains code for haplotype reconstruction, QTL mapping and plotting. Install qtl2 by running the following code in the R console.

~~~
install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
~~~
{: .r}

Make sure that the installation was successful by loading the qtl2 library. You shouldn't get any error messages.

~~~
library(qtl2)
~~~
{: .r}

## Data files and project organization

1. Create a "new project" within Rstudio. To do this, click on the transparent cube icon in the upper right corner of Rstudio. From the drop-down menu, choose "New Project". In the resulting pop-up window, choose "New Directory". When prompted for "Project Type" choose "New Project". Enter "mapping" in the field for "Directory name". Choose an appropriate directory, perhaps your "Desktop", for the field "Create project as subdirectory of:". Then click the button to "Create Project".



2. Create  a `data` folder to hold the data, a `scripts` folder to house your scripts, and a `results` folder to hold results. Within the `data` folder, make a folder called `derived_data` to hold derived data objects (vs. raw data).

Alternatively, you can use the R console to run the following commands for steps 1 and 2.

~~~
dir.create("./data")
dir.create("./scripts")
dir.create("./results")
dir.create("./data/derived_data")
~~~
{: .r}


3. Please download the data files from my email. Save them in the "derived_data" subdirectory of your new Rstudio project.



