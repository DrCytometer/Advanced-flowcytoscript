# installs the required packages

required.packages <- c(
  "digest", "dunn.test", "ggplot2", "ggridges", "ggrepel", "ggtext",
  "RColorBrewer", "Rtsne", "uwot", "dplyr", "tidyr", "RcppHNSW",
  "parallel", "data.table", "remotes", "BiocManager", "FastPG",
  "coda", "emmeans", "EmbedSOM", "ConsensusClusterPlus", "flowCore",
  "flowWorkspace", "readxl", "scattermore"
)

cran.packages <- c(
  "digest", "dunn.test", "ggplot2", "ggridges", "ggrepel", "ggtext",
  "RColorBrewer", "Rtsne", "uwot", "dplyr", "tidyr", "RcppHNSW",
  "parallel", "data.table", "remotes", "BiocManager", "coda", "emmeans",
  "EmbedSOM", "readxl", "scattermore"
)

bioconductor.packages <- setdiff(required.packages, cran.packages)

if ( length(setdiff(required.packages, rownames(installed.packages()))) !=0 ){
  
  install.packages(setdiff(cran.packages, rownames(installed.packages())))
  
  if ( length(setdiff(bioconductor.packages, rownames(installed.packages()))) !=0 ){
    BiocManager::install("sararselitsky/FastPG")
    BiocManager::install("ConsensusClusterPlus")
    BiocManager::install("flowCore")
    BiocManager::install("flowWorkspace")
  }
  
}

