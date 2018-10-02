# forecASD
Forecasting autism gene discovery with machine learning and genome-scale data. 

Genes are one of the most powerful windows into the biology of autism. However, less than 100 genes are currently viewed as having robust enough evidence to be considered true "autism genes". 

forecASD approaches autism gene discovery as a machine learning problem, rather than a genetic association problem, and uses genome-scale data from [STRING](http://version10.string-db.org/) and [Brainspan](http://www.brainspan.org/) as predictors for identifying further genes that have similar properties in the feature space compared to established autism risk genes. 

Preprint available here: https://www.biorxiv.org/content/early/2018/07/16/370601

# Using this code

All scripts are written in the [R programming language](https://cran.r-project.org/), and should be excuted in order indicated by name:
+ `01_load_data.R` = extracts, loads and formats external data sets
+ `02_network_models.R` = trains STRING and BrainSpan models
+ `03_ensemble_model.R` = trains model with a combination of existing ASD-gene scores and the network models

## Package Requirements
+ Matrix
+ RBGL
+ snow
+ randomForest

**NOTE** 
forecASD was executed using randomForest version **4.6-12** -- newer versions of this packages will return slightly different results.  You can install the version used in forecASD with the command: `devtools::install_version("randomForest", version = "4.6-12", repos = "http://cran.us.r-project.org")`

## sessionInfo()
```r
R version 3.5.1 (2018-07-02)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.5 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] randomForest_4.6-12 snow_0.4-2          RBGL_1.56.0         graph_1.58.0        BiocGenerics_0.26.0
[6] Matrix_1.2-14      

```
