library(dplyr)
library(Seurat)
library(patchwork)
library(hash)
library(ggplot2)

Global_pbmc <- NULL
Global_datacells <- NULL
Global_cluserNames <- NULL
Global_hashClus <- NULL #hash map of cells in the clusters
Namesgenes <- NULL
Namescalls <- NULL



#private methods:

# Load the PBMC dataset
makePBMC <- function(link2seurat_objec){
  pbmc =readRDS(file = link2seurat_objec)
  return(pbmc)
}
#To create matrix of the data:
makeDataMatrix <- function(){
  #to make the list of names of the genes.
  Namesgenes<<-NULL
  tempnamesgenes= Global_pbmc@assays$RNA@data@Dimnames[1]
  for (abc in tempnamesgenes) {
    Namesgenes <<- abc
  }
  #to make the list of names of the cells.
  Namescalls<<- NULL
  tempnamescalls= Global_pbmc@assays$RNA@data@Dimnames[2]
  for (abc in tempnamescalls) {
    Namescalls <<- abc
  }
  
  
  #to make the data table
  datacells= matrix(Global_pbmc@assays$RNA@data[c(Namesgenes),] #the data
                    , nrow = length(Namesgenes) #the number of rows= genes
                    , ncol = length(Namescalls))#the number of colums= cells
  # Naming rows
  rownames(datacells) = Namesgenes
  # Naming columns
  colnames(datacells) = Namescalls
  
  return(datacells)
}

makeHashCluster <- function(){
  
  hashClusters <- hash() 
  # set values
  for (nameC in Global_cluserNames){
    hashClusters[[nameC]] <- list()
  }
  
  for (x in Namescalls){
    clus = NULL
    for (temp in Global_pbmc@active.ident[x]){
      clus = temp
    }
    hashClusters[[ clus ]] = append(hashClusters[[ clus ]], x)
  }
  
  return(hashClusters)
}





#public methods:

Constructor <- function(link2seurat_objec = "output/pbmc3k_final.rds", 
            clusNames = c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B",
                             "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")){
  #restart necessary variables
  Global_pbmc <<- makePBMC(link2seurat_objec)
  print("pbmc finish!")
  Global_clusNames <<- clusNames
  print("clusters names finish!")
  Global_datacells <<- makeDataMatrix()
  print("datacells finish!")
  Global_hashClus <<- makeHashCluster()
  print("hashcluster finish!")
}

InpoScript_finalProject <- function(){
  infoScript = "Explain of methods and objects to use from this file:
  
  object:
   * datacells[key_geneName, key_cellID]:
   matrix that given 2 keys (can be vectors). the first key it's for the name of
   the gene. The second is for the cell ID. the matrix return the expression of
   the gene in the cell.
  
  
  methods:
    * Cell_in_cluster:
    get - name of cluster (or vector of names).
    return - all id of the cells in the cluster.
    
    * FindMarkers_New:
    get - the data, 2 vectors of the names of cluster.
    return - [1]the result from the function FindMarkers.
             [2]cells from the cluster in id1.
             [3]cells from the cluster in id2.
    
    * ExpressionGenes:
    get- cluster name, gene name, min exp. The defulte of the min exp it is 0.1 (like in the func FindMarkers), but it can be change.
    return- 3 things.
            [1]id of cells with expression to the given gene (more then minExp).
            [2]mean of how many cells had expression.
            [3]mean of the expression in the cells.
  
  "
    
  return(infoScript)
}

#get name of cluster (or vectors of names) and return all the cells in the cluster.
Cells_in_clusters <- function(keys){
  
  valH <- list()
  for (key in keys){
    valH <- unique(c(valH, Global_hashClus[[key]]))
  }
  return(valH)
}

# update to FindMarkers
FindMarkers_New <- function(data, id1 = 5, id2 = c(0, 3), mipct = 0.25){
  
  cluster.markers <- FindMarkers(data, ident.1 = id1, ident.2 = id2, min.pct = mipct)
  
  vec1 = Cells_in_clusters(id1)
  
  vec2 = Cells_in_clusters(id2)
  resFunc = list(cluster.markers, vec1, vec2)
  return(resFunc)
}

# function:
# get- cluster name, gene name, min exp.
# the defulte of the min exp it is 0.1 (like in the func FindMarkers), but it can be change.
# return- 3 things.
#         [1]id of cells with expression to the given gene (more then minExp).
#         [2]mean of how many cells had expression.
#         [3]mean of the expression in the cells.
ExpressionGenes <- function(clusName, geneName, minExp =0.1){
  cellsList <- Cells_in_clusters(clusName)
  expCells = list()
  for (ce in cellsList){
    if(Global_datacells[geneName, ce]>=minExp){
      expCells = c(expCells, ce)
    }
  }
  
  #create the dataframe of Cells and there expression.
  df_exCELLS <-data.frame(nameCell = character(), exp = numeric())
  for (ce in expCells){
    df_exCELLS <- rbind(df_exCELLS, data.frame(nameCell = ce, exp = Global_datacells[geneName, ce]))
  }
  
  p <- ggplot(df_exCELLS, aes(x = exp)) + geom_histogram()
  p <- p + labs(x = "expretion", y = "number of cells with this expetion", title = "exp of the cells")
  
  anser = list( df_exCELLS, length(expCells)/length(cellsList),  mean(df_exCELLS$exp),  p)
  return(anser)
}









