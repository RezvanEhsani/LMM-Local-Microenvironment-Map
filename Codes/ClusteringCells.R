

#' Clustering cells of interest across of all samples

#' @param data        dataframe of features returend by ExtractFeatures.R function
#' @param NCellTypes  Number of cell types
#' @param Nseed       Seed to pass on to given clustering method
#' @param k           number of clusters

#' There is a side function to determine optimal number of clusters (k) using 
#'                      "ConsensusClusterPlus" algorithm

#' @return a list contains 5 elements respectively
#'      (and colud called with the name, see the @example):
#'                  @ a dataframe same as feature dataframe but a extra columns called
#'                    "Clusters" which show clustering names for each cell of interest (CellClusterData)
#'                  @ a heatmap shows z-score (standardized values) of majority of each features in the Clusters (CellHeatmap)
#'                  @ Tsne plot of clusterd cells of interest (CellTsne)
#'                  @ a boxplot of number of different clustered cells of interest across the samples (CellBoxplot)
#'                  @ a barplot of row counts of clustered cells of interest (CellBarplot)

#'  @example
#'  First using "Nclusters_CellLevel" function we can find optimal number of clusters (k) and pass to main function "ClusteringCellLevel"
#'  Dat_Cell <- ClusteringCellLevel(Dat_Feature,
#'                                  NCellTypes=11,
#'                                  k=8,
#'                                  Nseed=1234)
#' Now dataframe of clusters cloud called by Dat_Cell$CellClusterData



# Set some color
Cell_Colors <- colors()[c(27,202,107,62,256,
                           147,79,386,452,491,539,117,142,115,
                           630,636,506,442,55,87)]


Nclusters_CellLevel <- function(Dat_Feature,
                                NCellTypes,
                                Nseed=1234,
                                PlotName){

    Nfeatues <- 6 + 2*NCellTypes + 12 # 6 circel information + celltypes + connectivity of cell types + 12 graph features
    NExistClass <- ncol(Dat_Feature) - Nfeatues
    L1 <- 7
    L2 <- ncol(Dat_Feature)-NExistClass
    Dat_Feature <- na.omit(Dat_Feature)
    Dat_Feature[L1:L2] <- sapply(Dat_Feature[L1:L2], as.numeric)
    Dat_Feature[L1:L2] <- Dat_Feature[L1:L2]/rowSums(Dat_Feature[L1:L2])

    Dat_Feature <- Dat_Feature[,c(L1:L2)]
    set.seed(Nseed)
    Dat_km <- kmeans(Dat_Feature, centers=50, iter.max = 100)
    Dat_Feature$km_cluster <- Dat_km$cluster

    Dat_Feature <- data.table(Dat_Feature)
    Dat_Ave_km <- data.frame(Dat_Feature[, sapply(.SD, function(x) list(mean(x))), by=km_cluster])

    Dat_Ave_km_1 <- Dat_Ave_km[,2:(L2-5)]
    maxK = 15 # maximum number of clusters to try
    mc = ConsensusClusterPlus(as.matrix(t(Dat_Ave_km_1)),
                              maxK=maxK,reps=50,
                              pItem=0.8,
                              pFeature=1,
                              clusterAlg="km", # "hc","km","kmdist","pam"
                              title=PlotName,
                              distance="euclidean",  #"euclidean","pearson","spearman","binary","maximum","canberra","minkowski"
                              innerLinkage="complete",
                              seed=Nseed,
                              plot='pdf')
    }

#Main function
ClusteringCellLevel <- function(Dat_Feature,
                                     NCellTypes,
                                     k,
                                     Nseed=1234){

    Nfeatues <- 6 + 2*NCellTypes + 12 # 6 circel information + celltypes + connectivity of cell types + 12 graph features
    NExistClass <- ncol(Dat_Feature) - Nfeatues

    L1 <- 7
    L2 <- ncol(Dat_Feature)-NExistClass
    Dat_Feature <- na.omit(Dat_Feature)

    Dat_Feature[L1:L2] <- sapply(Dat_Feature[L1:L2], as.numeric)
    Dat_Feature[L1:L2] <- Dat_Feature[L1:L2]/rowSums(Dat_Feature[L1:L2])

    F_dat <- Dat_Feature[,c(L1:L2)]

    set.seed(Nseed)
    res.km <- kmeans(F_dat, centers=k, iter.max = 100)
    F_dat$km_cluster <- res.km$cluster
    F_dat$Clusters <- paste("C", F_dat$km_cluster, sep="")

    Dat_Feature$Clusters <- F_dat$Clusters
    #Clustered data
    Dat_Feature$Clusters <- factor(Dat_Feature$Clusters ,levels=sort(unique(Dat_Feature$Clusters)))

    Dat_Feature1 <- Dat_Feature[L1:Nfeatues]
    Dat_Feature1$Clusters <- Dat_Feature$Clusters
    Dat_heat <- data.table(Dat_Feature1)
    
    Dat_heat <- data.frame(Dat_heat[, sapply(.SD, function(x) list(mean(x))), by=Clusters])
    rownames(Dat_heat) <- Dat_heat$Clusters

    Dat_heat <- Dat_heat[,c(2:ncol(Dat_heat))]
    Dat_heat <- scale(Dat_heat)
    Dat_heat[Dat_heat > 2] =2
    Dat_heat[Dat_heat < -2] =-2

    # Heatmap
    pheat <- pheatmap(as.matrix(Dat_heat),
                      name = "z-score",
                      border_color = "grey30",
                      gaps_col=c(NCellTypes,2*NCellTypes),
                      fontsize = 7,
                      cellwidth = 11,
                      cellheight=11,
                      show_rownames = TRUE,
                      show_colnames = T,
                      fontsize_row = 7,
                      fontsize_col =7,
                      display_numbers = F,
                      cluster_rows=F,
                      cluster_cols=F
                      )
    #Tsne plot    
    tmp <- Dat_Feature1
    myClusters <- Dat_Feature1$Clusters 

    duplicates <- duplicated(tmp) | duplicated(tmp, fromLast = TRUE)
    tsne_out<-Rtsne(tmp[!duplicates,1:ncol(tmp)],max_iter=1000,seed=Nseed, perplexity=30)
    tsne_data <- tmp[!duplicates,]
    tsne_data$TSNE1 <- tsne_out$Y[,1]
    tsne_data$TSNE2 <- tsne_out$Y[,2]
    tsne_data$Clusters <- myClusters[!duplicates]

    tsne_plot <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = Clusters )) +
                 geom_point(size = 1.5,alpha=0.8) +
                 theme_bw()+
                 ylab('t-SNE2')+
                 xlab('t-SNE1')+
                 #scale_color_brewer(palette="Set2")+
                 scale_color_manual(values=Cell_Colors)+
                 theme(axis.text.y = element_text( size = 12 ),
                 axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 12),
                 legend.position = "bottom",
                 legend.text=element_text(size=12))+
                 guides(color=guide_legend(nrow=4,override.aes = list(size=4)))

    # Box plot
    Melted <- melt(table(Dat_Feature$sample_id, Dat_Feature$Clusters))
    colnames(Melted) <- c("sample_id", "Clusters", "Counts")

    p1 <- ggplot(Melted, aes(x=Clusters, y=Counts, fill=Clusters)) +
          geom_boxplot(width=0.5)+
          scale_fill_manual(values=Cell_Colors)+
          theme_classic()+
          labs(title="", y="")+
          theme(legend.position = "none")+
          coord_flip()
                                   
    Melted2 <- melt(table(Dat_Feature$Clusters))
    colnames(Melted2) <- c("Clusters", "NumberGraphs")

    # Bar plot
    p2 <-  ggplot(Melted2, aes(x=Clusters, y=NumberGraphs, fill=Clusters)) + 
           geom_bar(stat="identity",width=0.7)+
           scale_fill_manual(values=Cell_Colors)+
           labs(title="", y="")+
           theme_classic()+  
           theme(legend.position = "none")+
           coord_flip()
                                   
    return(list(CellClusterData=Dat_Feature,
                CellHeatmap=pheat,
                CellTsne=tsne_plot,
                CellBoxplot=p1,
                CellBarplot=p2
                ))
    }

