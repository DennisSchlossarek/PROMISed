# Comparing AHA2 Interactors
library(VennDiagram)
library(eulerr)
library(ggplot2)
setwd("H:/10. Publications/PROMISED/^^ Manuscript PROMISed/DataAnalysis_Gilbert_Schulze/AHA2 Interactors")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") #Switch of log-files while creating venn-diagrams

GetIntersections <- function(x){
  
  overlap <- calculate.overlap(x)
  
  if(length(x) == 2){
    names(overlap) <- c("a1", "a2", "a12")
  } else if(length(x) == 3){
    names(overlap) <- c("a123", "a12", "a13", "a23", "a1", "a2", "a3")
  } else if(length(x) == 4){
    names(overlap) <- c("a1234", "a123", "a124", "a134", "a234", "a12", "a13", "a14", "a23", "a24", "a34", "a1", "a2", "a3", "a4")
  }
  
  maxrow <- c()
  k <- 1
  while(k <= length(overlap)){
    
    maxrow[k] <- length(overlap[[k]])
    
    k <- k + 1
  }
  
  overlap_matrix <- matrix(nrow = max(maxrow), ncol = length(overlap), NA)
  colnames(overlap_matrix) <- names(overlap)
  
  k <- 1
  while(k <= length(overlap)){
    
    if(length(c(overlap[[k]])) > 0){
      overlap_matrix[(1:maxrow[k]), k] <- c(overlap[[k]])  
    }
    k <- k + 1
  }
  
  overlap_matrix
  
}

mycol = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

GilbertSchulze <- read.delim("AHA2_interactors_gilbert.txt", header = TRUE)
GilbertSchulze_pulldown <- read.delim("AHA2_pulldown_gilbert.txt", header = TRUE)

GilbertSchulze_HC <- GilbertSchulze[grep(x = GilbertSchulze$Interaction.type, pattern = "HC"),]


PROMISed_rep1 <- read.delim("Rep1_AHA2_edges.txt", header = TRUE)
PROMISed_rep2 <- read.delim("Rep2_AHA2_edges.txt", header = TRUE)

# Filter for High Confidence
PROMISed_rep1 <- PROMISed_rep1[(PROMISed_rep1$weight >= 0.89),]
PROMISed_rep2 <- PROMISed_rep2[(PROMISed_rep2$weight >= 0.89),]

PROMISed_rep1 <- (unique(data.frame(matrix(unlist(strsplit(PROMISed_rep1[,2], split = "_K_")), ncol = 2, byrow = TRUE))[,1]))
PROMISed_rep2 <- (unique(data.frame(matrix(unlist(strsplit(PROMISed_rep2[,2], split = "_K_")), ncol = 2, byrow = TRUE))[,1]))

AHA2_overlap <- list(Gilbert = na.omit(GilbertSchulze_HC$Protein.ID),
                     PROMISed_1 = na.omit(PROMISed_rep1),
                     PROMISed_2 = na.omit(PROMISed_rep2))

{png("AHA2_interactors_overlap_HC.png", height = 210, width = 210, units = "mm", res = 300, bg = "transparent")
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))

venn_over <- venn.diagram(AHA2_overlap, 
                          filename = NULL,
                          output = TRUE,
                          
                          imagetype = "tiff",
                          height = 600,
                          width = 600,
                          resolution = 256,
                          
                      #    main = "Predicted AHA2 Interactors",
                      #    sub = "",
                          main.cex = 2,
                          
                          category.names = c("Gilbert&Schulze","PROMISed Rep1", "PROMISed Rep2"),
                          cat.cex = 1.5,
                          cat.default.pos = "outer",
                          cat.fontfamily = "sans",
                          cex = 1.5,
                          
                          lwd = 3,
                          col = mycol[1:length(AHA2_overlap)],
                          fill = alpha(mycol[1:length(AHA2_overlap)], 0.6)
)

grid.draw(venn_over)

dev.off()}

png("AHA2_interactors_overlap_HC_euler.png", height = 210, width = 210, units = "mm", res = 300, bg = "transparent")
plot(euler(AHA2_overlap),
     fills = list(fill = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[1:length(AHA2_overlap)], alpha = 0.7),
     quantities = list( cex = 2),
     names = list(cex = 3),
     legend = list(labels = names(AHA2_overlap), cex = 1.5))
dev.off()


my_intersections <- GetIntersections(AHA2_overlap)
colnames(my_intersections) <- c("Center", "GS_PROMISed1", "GS_PROMISed2", "PROMISed1_PROMISed2", "GS", "PROMISed1", "PROMISed2")

write.table(my_intersections, "AHA2_interactions_intersection_HC.txt", sep = "\t")

list_pd <- list(SEC = na.omit(my_intersections[,1]), 
                Pulldown = GilbertSchulze_pulldown$IDENTIFIER)

{png("AHA2_interactors_SEC_pulldown_HC.png", height = 210, width = 210, units = "mm", res = 300, bg = "transparent")
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))

venn_pull <- venn.diagram(list_pd, 
             filename = NULL,
             output = TRUE,
             
             imagetype = "tiff",
             height = 600,
             width = 600,
             resolution = 256,
             
          #   main = "Predicted AHA2 Interactors",
          #   sub = "",
             main.cex = 2,
             
             category.names = c("SEC", "Pulldown"),
             cat.cex = 1.5,
             cat.default.pos = "outer",
             cat.fontfamily = "sans",
             cex = 1.5,
             
             lwd = 3,
             col = mycol[1:2],
             fill = alpha(mycol[1:2], 0.6)
)

grid.draw(venn_pull)

dev.off()}

png("AHA2_interactors_SEC_pulldown_HC_euler.png", height = 210, width = 210, units = "mm", res = 300, bg = "transparent")
plot(euler(list_pd),
     fills = list(fill = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[1:length(list_pd)], alpha = 0.7),
     quantities = list( cex = 2),
     names = list(cex = 3),
     legend = list(labels = names(list_pd), cex = 1.5))
dev.off()

overlap_GS_PROMISed <- list(PROMISed = na.omit(c(my_intersections[,1], my_intersections[,4])), 
                            GilbertSchulze = na.omit(GilbertSchulze_HC$Protein.ID), 
                            Pulldown = GilbertSchulze_pulldown$IDENTIFIER)

my_intersections2 <- GetIntersections(overlap_GS_PROMISed)
colnames(my_intersections2) <- c("Center", "PROMISed_GS", "PROMISed_PD", "GS_PD", "PROMISed", "GS", "PD")

write.table(my_intersections2, "AHA2_PROMISed_GS_intersection_HC.txt", sep = "\t")

{png("AHA2_interactors_PROMISed_pulldown_HC.png", height = 148, width = 148, units = "mm", res = 300, bg = "transparent")
  grid.newpage()
  pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
  
  venn_pull <- venn.diagram(overlap_GS_PROMISed, 
                            filename = NULL,
                            output = TRUE,
                            
                            imagetype = "tiff",
                            height = 600,
                            width = 600,
                            resolution = 256,
                            
                            #   main = "Predicted AHA2 Interactors",
                            #   sub = "",
                            main.cex = 2,
                            
                            category.names = c("PROMISed", "Gilbert&Schulze", "Pull-Down"),
                            cat.cex = 1.5,
                            cat.default.pos = "outer",
                            cat.fontfamily = "sans",
                            cex = 1.5,
                            
                            lwd = 3,
                            col = mycol[1:3],
                            fill = alpha(mycol[1:3], 0.6),
                            rotation.degree = 0
  )
  
  grid.draw(venn_pull)
  
  dev.off()}

png("AHA2_interactors_PROMISed_pulldown_HC_euler.png", height = 130, width = 210, units = "mm", res = 300, bg = "transparent")
plot(venn(overlap_GS_PROMISed),
     fills = list(fill = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[1:length(overlap_GS_PROMISed)], alpha = 0.7),
     quantities = list( cex = 1.5),
     names = list(cex = 3),
     legend = list(labels = names(overlap_GS_PROMISed), cex = 1.5))
dev.off()

