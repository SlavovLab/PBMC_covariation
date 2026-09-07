########################################################################
## Install and load packges using renv to manage library
########################################################################




############ Install renv if not already installed:

# install.packages("renv")





############ Use renv to store and save package versions:

# renv::init()





############ Install packages:

## Note:
# QuantQC R package is needed for mapping single-cell proteomic data from the sample prep used, termed nPOP, via the CellenONE platform
# QuantQC is typically being continuously updated via Github, so special package installation instructions are provided to help the user replicate this analysis
# Below are specific packages that needed versions required for QuantQC that were beyond the current CRAN release and for some, the developer version is installed
# Usually the solution is to download the current release/developer version of packages if the user replicating this analysis faces any install challenges

# install.packages(c("devtools", "remotes", "Seurat", "BiocManager", "dplyr",
#                  "reshape2", "ggplot2", "ggpubr", "patchwork",
#                  "gghighlight", "viridis", "ggpointdensity",
#                  "tidyr", "stringr", "ggridges", "rmarkdown",
#                  "seqinr", "dtplyr", "dataPreparation", "pheatmap",
#                  "ggplotify", "plotly", "Matrix", "ineq", "cluster",
#                  "tidyverse", "randomForest", "ggfortify", "factoextra",
#                  "caret", "pROC", "fastshap", "stats", "missMDA", "MASS", "rliger", "readr", "ggrepel", "ComplexHeatmap",
#                  "msigdbr", "ggokabeito", "tibble", "data.table", "limma", "circlize", "gprofiler2", "igraph", "STRINGdb", "uwot", "fastICA",
#                  "ppcor", "vioplot", "lmerTest", "ggh4x", "ggrastr"))
# 
# library(remotes)
# remotes::install_version("psych", version = "2.5.3", repos = "http://cran.us.r-project.org") # Update this to whatever version is required for QuantQC installed next
# remotes::install_github("HenrikBengtsson/matrixStats", ref="develop") # Gets current version for dependency. At the time of writing this script, it is version 1.5.0
# remotes::install_github("jeffreyevans/yaImpute") # Gets current version for dependency. At the time of writing this script, it is version 1.0.35
# remotes::install_github("https://github.com/vdemichev/diann-rpackage")
# 
# 
# library(BiocManager)
# BiocManager::install("genefilter", force = TRUE)
# BiocManager::install("sva", force = TRUE)
# BiocManager::install("mzR", force = TRUE)
# BiocManager::install("biomaRt", force = TRUE)
# BiocManager::install('org.Hs.eg.db')
# BiocManager::install('reactome.db')
# BiocManager::install('fgsea')


# library(devtools)
# devtools::install_github("https://github.com/SlavovLab/QuantQC")
# devtools::install_github("welch-lab/RcppPlanc")





############ Load libraries:

# library(data.table)
library(rmarkdown)
library(BiocManager)
library(org.Hs.eg.db)
library(reactome.db)
library(biomaRt)
library(genefilter)
library(fgsea)
library(limma)
library(sva)
library(mzR)
library(scales)
library(readr)
library(seqinr)
library(Matrix)
library(reshape2)
library(dtplyr)
library(dataPreparation)
library(tidyverse)
library(conflicted)
library(matrixStats)
library(psych)
library(missMDA)
library(stats)
library(ineq)
library(cluster)
library(ppcor)
library(vioplot)
library(randomForest)
library(caret)
library(pROC)
library(fastshap)
library(fastICA)
library(Seurat)
library(rliger)
library(RcppPlanc)
library(hdWGCNA)
library(uwot)
library(diann)
library(QuantQC)
library(yaImpute)
library(msigdbr)
library(gprofiler2)
library(STRINGdb)
library(igraph)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggridges)
library(ggpointdensity)
library(gghighlight)
library(ggfortify)
library(factoextra)
library(ggokabeito)
library(viridis)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(ggplotify)
library(plotly)
library(grid)
library(lmerTest)
library(ggh4x)
library(ggrastr)

conflict_prefer("select",     "dplyr")
conflict_prefer("filter",     "dplyr")
conflict_prefer("left_join",  "dplyr")
conflict_prefer("right_join", "dplyr")
conflict_prefer("full_join",  "dplyr")
conflict_prefer("mutate",     "dplyr")
conflict_prefer("summarise",  "dplyr")
conflict_prefer("summarize",  "dplyr")
conflict_prefer("arrange",    "dplyr")
conflict_prefer("rename",     "dplyr")
conflict_prefer("count",      "dplyr")
conflict_prefer("slice",      "dplyr")
conflict_prefer("desc",       "dplyr")
conflict_prefer("vars",       "dplyr")
conflict_prefer("first",      "dplyr")
conflict_prefer("cor",        "stats")
conflict_prefer("var",        "stats")
conflict_prefer("intersect",  "base")
conflict_prefer("setdiff",    "base")
conflict_prefer("union",      "base")
conflict_prefer("unname",     "base")
conflict_prefer("margin",     "ggplot2")
conflict_prefer("alpha",      "ggplot2")





############ Once packages are successfully installed, save the version via renv:snapshot()

# renv::snapshot()





############ If the user of this script ever wants to share or reproduce this environment elsewhere, run this renv function:

# renv::restore()




