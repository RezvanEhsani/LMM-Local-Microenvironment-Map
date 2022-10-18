
#' Clustering samples based on proportion of different clusters of cells of interest


#' @param data         dataframe returend by ClusteringCells.py function
#' @param NCellTypes   Number of cell types
#' @param Nseed        Seed to pass on to given clustering method
#' @param k            number of clusters
#' @param ExistClasses Set TRUE if there is some the existing classification to be campered

#' There is a side function to determine optimal number of clusters (k) using 
#'                      "ConsensusClusterPlus" algorithm



#' @return a list contains 5 elements respectively
#'      (and colud called with the name, see the @example):
#'                  @ a dataframe same as input dataframe but a extra columns called 
#'                    "SampleCluster" which show clustering names for each sample (SampleClusterData)
#'                  @ a stack bar plot of proportion of different clusters of cells of interest across samples (StackBar)        
#'                  @ a heatmap shows z-score (standardized values) of majority of each cell of interest clusters in samples (SampleHeatmap)
#'                  @ Tsne plot of clusterd samplest (SampleTsne)
#'                  @ Alluvium plot, if  ExistClasses=TRUE to see agreement between new and existing clusters (SampleAlluvium)



#'  @example
#'  First using "Nclusters_SampleLevel" function we can find optimal number of clusters (k) and pass to main function "ClusteringSampleLevel"
#'  Dat_Sample <- ClusteringSampleLevel(Dat_Cells,
#'                                  NCellTypes=11,
#'                                  k=4,
#'                                  Nseed=1234)
#' Now the dataframe of sample clusters cloud be called by Dat_Sample$SampleClusterData





# Set some color
# Put all the color values (in hex format) from Dark2 into a vector
library(RColorBrewer)
MyDark2 <- brewer.pal(8,"Dark2")
Sample_Colors <- c(colors()[c(133,642)], MyDark2)

#
Nclusters_SampleLevel <- function(Dat_Cells,
                                  maxK = 5,
                                  PlotName,
                                  Nseed=1234){

    Dat_Sample <- dcast(Dat_Cells, sample_id  ~ Clusters)
    L <- ncol(Dat_Sample)
    Dat_Sample[,2:L] <- Dat_Sample[,2:L]/rowSums(Dat_Sample[,2:L])

    mchc = ConsensusClusterPlus(as.matrix(t(Dat_Sample[,2:L])),
                                maxK = maxK,
                                clusterAlg = 'km',
                                distance = "euclidean",
                                reps = 50,
                                pItem = 0.8,
                                pFeature = 1,
                                title = PlotName,
                                innerLinkage = "complete",
                                seed = Nseed,
                                plot = 'pdf')
    }


#
alluvium_plot <- function(Dat, ExistClassName){

    ExistClass <- colnames(Dat)[3]
    ExistClassName <- String(ExistClassName)
    ExistClassName <- ensym(ExistClassName)

    dat <- Dat %>%
      group_by(SampleCluster, !!ExistClassName) %>%
      summarise(freq = n()) %>%
      ungroup()

    k_sample = length(unique(Dat$SampleCluster))
    ggallu <- ggplot(data = dat,aes(axis1 = !!ExistClassName, axis2 = SampleCluster, y = freq)) +
          geom_alluvium(aes(fill = SampleCluster),width =  1/16,alpha=.7) + #, knot.pos = 0
          geom_stratum(width =  1/16 , fill = "gray50") +
          geom_text(stat = "stratum", aes(label = after_stat(stratum)),color = "White",angle = 0) +
          scale_x_discrete(limits = c(ExistClass, "SampleCluster"),expand = c(0.15, 0.05)) +
          scale_fill_manual(values = Sample_Colors[1:k_sample]) +
          theme_void()+
          coord_flip( )+
          theme(legend.position="non")

    return(ggallu)
    }

        
ClusteringSampleLevel <- function(Dat_Cells,
                                  NCellTypes,
                                  k_sample,
                                  Nseed=1234,
                                  ExistClasses=TRUE){
                                      
    Nfeatues <- 6 + 2*NCellTypes + 12
    Dat_Cells <- Dat_Cells[-c(2:Nfeatues)]

    Dat_Sample <- dcast(Dat_Cells, sample_id  ~ Clusters)
    L <- ncol(Dat_Sample)
    Dat_Sample[,2:L] <- Dat_Sample[,2:L]/rowSums(Dat_Sample[,2:L])
    Dat_Sample$sample_id <- factor(Dat_Sample$sample_id,unique(sort(Dat_Sample$sample_id)))

    # Stack bar
    melted <- melt(Dat_Sample, id = c("sample_id"))
    StackBar <- ggplot(data=melted, aes(x=sample_id, y=value,fill=variable)) +
                geom_bar(stat="identity")+
                labs(title = "",y = "", x = "Samples")+
                scale_fill_manual(values=Cell_Colors)+
                theme_bw()
    #
    rownames(Dat_Sample) <- Dat_Sample$sample_id

    set.seed(Nseed)
    res.km <- kmeans(Dat_Sample[,2:L], k_sample, iter.max = 100, nstart = 25)
    
    Dat_Sample$SampleCluster <- paste("SC", res.km$cluster, sep="")
    Dat_Sample$SampleCluster <- factor(Dat_Sample$SampleCluster, levels=sort(unique(Dat_Sample$SampleCluster)))

    Da <- unique(Dat_Cells[1:(ncol(Dat_Cells)-1)])
    Dat_Sample <- merge(Dat_Sample, Da, by=c("sample_id") )
    Dat_Sample <- Dat_Sample[order(Dat_Sample$SampleCluster),]

    Dat_heat <- Dat_Sample[,2:L]
   
    rownames(Dat_heat) <- Dat_Sample$sample_id

    Dat_heat <- scale(Dat_heat)
    Dat_heat[Dat_heat > 2] =2
    Dat_heat[Dat_heat < -2] =-2
    
    annotation_row <- data.frame(SampleCluster=factor(Dat_Sample$SampleCluster))
    rownames(annotation_row) <- Dat_Sample$sample_id

    mycolor <- Sample_Colors[1:k_sample]
    names(mycolor) <- unique(Dat_Sample$SampleCluster)
    annotation_colors <- list(SampleCluster = mycolor)

    # Heatmap
    pheat_samples <- pheatmap(as.matrix(Dat_heat),
                          name = "z-score",
                          fontsize = 7,
                          border_color = "grey60",
                          fontsize_row = 7,
                          fontsize_col =7,
                          cellwidth = 12,
                          display_numbers = F,
                          cluster_rows = F,
                          cluster_cols = F,
                          annotation_row = annotation_row ,
                          annotation_colors = annotation_colors
                          )

    #Tsne plot    
    tmp <- Dat_Sample
    myClusters <- Dat_Sample$SampleCluster

    duplicates <- duplicated(tmp) | duplicated(tmp, fromLast = TRUE)
    tsne_out<-Rtsne(tmp[!duplicates,1:ncol(tmp)],max_iter=1000,seed=Nseed, perplexity=5)
    tsne_data <- tmp[!duplicates,]
    tsne_data$TSNE1 <- tsne_out$Y[,1]
    tsne_data$TSNE2 <- tsne_out$Y[,2]
    tsne_data$SampleCluster <- myClusters[!duplicates]

    tsne_plot <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = SampleCluster )) +
                 geom_point(size = 2.5) +
                 theme_bw()+
                 ylab('t-SNE2')+
                 xlab('t-SNE1')+
                 scale_color_manual(values=Sample_Colors)+
                 theme(axis.text.y = element_text( size = 12 ),
                 axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 12),
                 legend.position = "bottom",
                 legend.text=element_text(size=12))+
                 guides(color=guide_legend(nrow=2,override.aes = list(size=4)))


    # alluvium plot
    if(ExistClasses == TRUE){
        
        alluvium_Plots <- list()
        LC = length(unique(Dat_Cells$Clusters)) 
        EistClassA <- Dat_Sample[-(2:(LC+1))]

        for(i in 3:ncol(EistClassA)){
            EistClassAi <- EistClassA[,c(1,2,i)]
            EistClassAi[,3] <- factor(EistClassAi[,3])
            
            pi <- alluvium_plot(EistClassAi, ExistClassName=colnames(EistClassAi)[3])
            
            alluvium_Plots[[i-2]] <- pi
            }

        Plot_alluvium <- ggarrange(plotlist = alluvium_Plots, nrow = 1, ncol = 1)

        return(list(SampleClusterData = Dat_Sample,
                    StackBar = StackBar,
                    SampleHeatmap = pheat_samples,
                    SampleTsne = tsne_plot,
                    SampleAlluvium = Plot_alluvium
                ))
    }
    else return(list(SampleClusterData = Dat_Sample,
                     StackBar = StackBar,
                     SampleHeatmap = pheat_samples,
                     SampleTsne = tsne_plot
                ))
}
