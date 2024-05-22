## James Fulcher & Kyle Nadeau
## Pacific Northwest National Laboratory
## Single Cell Proteomics

## Imports
library(mclust)
library(tidyverse)
library(stats)
library(cluster)

#############################################################
#############Cluster comparison function #################
#############################################################
Clust_compare <-function(data, 
                         clusterNumbers= c(5,10),
                         nameAlgorithm = c('mclust'),
                         models = c('VEI','VII', 'EEI',
                                     'EII','EVI'),
                         shrinkage_values = c(0.01,0)) {
  
  dfList <- list()
  print('Initializing...')
  data <- t(scale(t(data)))
  for (algorithm in nameAlgorithm) {
  for (clusters in clusterNumbers) {
    ##MCLUST (REMEMBER SPECIFIC TYPES)
    if (algorithm == 'mclust'){
      print('Running mclust')
      for(model in models) {
        set.seed(420)
        mclust::mclust.options(subset = 4000)
        default_shrinkage <- shrinkage_values[1]
        tryCatch({
          result <- mclust::Mclust(data, 
                                   modelNames = model, 
                                   G = clusters,
                                   prior = priorControl(shrinkage = default_shrinkage))
      
          }, error = function(e) {
          cat("Error with shrinkage=", default_shrinkage,"\n")
            result <- NULL
        })
        if(is.null(result)) {
          for(shrink_val in shrinkage_values[-1]){
            tryCatch({
              result <- mclust::Mclust(data, 
                                       modelNames = model, 
                                       G = clusters,
                                       prior = priorControl(shrinkage =shrink_val))
            }, error = function(e) {
              cat("Error with shrinkage=", shrink_val,"\n")
            })
          }
        }
        
        df <- data.frame(
          Gene = names(result$classification),
          Cluster = result$classification,
          uncertainty = result$uncertainty)
        df_name <- paste( model, "_clusters_", clusters, sep = "")
        dfList[[df_name]] <- df
      }
    }
                      }
  }
  dfList <- dfList[sapply(dfList, nrow)>0]
  return(dfList)
}
