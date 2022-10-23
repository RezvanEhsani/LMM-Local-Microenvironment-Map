
#' Map clustered cells of interest on the real Images(Samples)

#' @param data                  a dataframe of samples, cells, markers and spatial information
#' @param DatCells              dataframe returend by ClusteringCells.R function
#' @param SampleID_col column   Index of Sample ID or Name
#' @param CellID_col   column   Index of Cell ID
#' @param X_col and Y_col       Column indexs of X and Y position of cells


#' a ggarrange type of ggplot 



library(ggpubr)
library(ggrepel)


# Set some color
Cell_Colors <- colors()[c(27,202,107,62,256,
                           147,79,386,452,491,539,117,142,115,
                           630,636,506,442,55,87)]

MapClusteredCellsOnSamples <- function(Data,
                                      DatCells,
                                      SampleID_col,
                                      CellID_col,
                                      MetaCellType_col,
                                      X_col,
                                      Y_col){
                                          
        colnames(Data)[SampleID_col] <- "sample_id"
        colnames(Data)[CellID_col] <- "cell_id"
        colnames(Data)[MetaCellType_col] <- "MetaCellType"
        colnames(Data)[X_col] <- "Pos_X"
        colnames(Data)[Y_col] <- "Pos_Y"

        colnames(DatCells)[3] <- "cell_id"
        DatCells <- DatCells[c("sample_id", "cell_id","Clusters")]

        Data <- merge(Data, DatCells, by=c("sample_id", "cell_id"), all.x=TRUE)

        Data$Clusters <- as.character(Data$Clusters)
        Data$Clusters[is.na(Data$Clusters)] <- 0
        Data$Clusters[Data$Clusters==0] <- Data$MetaCellType
        Data$Clusters[Data$Clusters==0] <- "NonTumor"
        Data$Clusters[Data$Clusters==1] <- "Tumor"

        Data$Clusters <- factor(Data$Clusters, levels=c('Tumor', 'NonTumor', levels(unique(DatCells$Clusters))))
        cols <- c('#C0C0C0','#F9BFBF',Cell_Colors[1:length(unique(DatCells$Clusters))])
        names(cols) <- c('Tumor', 'NonTumor',levels(unique(DatCells$Clusters)))

        ID <- unique(Data$sample_id)
        plot_all <- list()

        for(i in 1:length(ID)){

            Datai <- Data[(Data$sample_id ==ID[i]),]

            p <- ggplot(Datai, aes(Pos_X, Pos_Y, color=Clusters))+
                 geom_point(size=0.2, alpha=0.6)+
                 scale_colour_manual(values = cols) +
                 guides(color = guide_legend(override.aes = list(size = 3) ) )+
                 theme_bw()+
                 labs(title=ID[i], x ="", y = "") + 
                 theme(legend.position="bottom")+
                 theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                   axis.title.x = element_blank(), axis.title.y = element_blank(),
                   legend.key = element_blank(),
                   plot.title = element_text(size=7),
                   legend.text = element_text(size = 8)
                   )
        plot_all[[i]] <- p
        }

    PlotAll <- ggarrange(plotlist=plot_all, nrow = 3,ncol = 4,common.legend = TRUE, legend="bottom")
    return(PlotAll)

    }
