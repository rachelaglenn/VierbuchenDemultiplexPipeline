library(ggplot2)
library(ggrepel)
library(outliers)
library(EnvStats)
#########################doubletdetection#############################
doublet_small <- read.csv("C_P2-5-11-25-30_detect_doublets.txt", sep = "\t", skip = 1)

doubletdata <- read.csv("C_P2-5-11-25-30_detect_doublets_all.txt", sep = "\t", skip = 1)

plot_doublet <- function(db_data,db_small, barcode, show=T){
  print(barcode)
  db_filt <- doubletdata[doubletdata$cell == barcode,]
  plot_data <- data.frame(likelihood = c(db_filt$sampleOneLikelihood[1],db_filt$sampleTwoLikelihood))
  plot_data$sample <- c(db_filt$sampleOne[1],db_filt$sampleTwo)
  plot_data <- plot_data[order(plot_data$likelihood, decreasing = T),]
  plot_data$index <- c(1:nrow(plot_data))
  
  bestsample <- doublet_small[doublet_small$cell == barcode,"bestSample"]
  db_pval <- doublet_small[doublet_small$cell == barcode,"doublet_pval"]
  
  plot <- ggplot(data=plot_data, aes(x = index, y=likelihood, label = sample)) +
    geom_point() +
    geom_line() + 
    geom_text_repel() + 
    ggtitle(paste0(barcode,"---- best Sample: ", bestsample, "  doublet_pval: ", db_pval))
  
  if(show){
    print(plot)
  }
  return(plot)
}

#singlet
plot_doublet(doubletdata, doublet_small, "CATGAGTAGGCCTGCT")

#doublet
plot_doublet(doubletdata,doublet_small, "AAAGAACTCGCGATCG")

#noise
plot_doublet(doubletdata,doublet_small, "AAAGGGCGTTGCTGAT")


pdf("doublet_plots.pdf", width = 8, height = 6)
for(cell in unique(doubletdata$cell)[1:1000]){
  plot_doublet(doubletdata,doublet_small, cell)
  
}
dev.off()

##############own method

assign_cell <- function(db_data,db_small, barcode, show=T){
  #print(barcode)
  db_filt <- doubletdata[doubletdata$cell == barcode,]
  plot_data <- data.frame(likelihood = c(db_filt$sampleOneLikelihood[1],db_filt$sampleTwoLikelihood))
  plot_data$sample <- c(db_filt$sampleOne[1],db_filt$sampleTwo)
  plot_data <- plot_data[order(plot_data$likelihood, decreasing = T),]
  plot_data$index <- c(1:nrow(plot_data))
  
  bestsample <- plot_data$sample[1]
  
  # test <- grubbs.test(plot_data$likelihood)
  # test
  #db_pval <- doublet_small[doublet_small$cell == barcode,"doublet_pval"]
  
  test <- rosnerTest(plot_data$likelihood,
                     k = 2
  )
  test <- test$all.stats
  test$sample <- plot_data$sample[test$Obs.Num]
  
  # if(sum(test$Outlier) == 1 & test$sample[1] == bestsample){
  #   type <- "singlet"
  # }else{
  #   if(sum(test$Outlier) == 2 & all(test$sample[1:2] == plot_data$sample[1:2])){
  #     if(abs(plot_data$likelihood[2] - plot_data$likelihood[1]) < abs(plot_data$likelihood[3] - plot_data$likelihood[2])){
  #       type <- "doublet"
  #     }
  #     
  #   }else{
  #     type <- "unknown"
  #   }
  # }
  
  type <- "unknown"
  if(sum(test$Outlier) == 2 & all(test$sample[1:2] == plot_data$sample[1:2])){
    if(abs(plot_data$likelihood[2] - plot_data$likelihood[1]) < abs(plot_data$likelihood[3] - plot_data$likelihood[2])){
      type <- "doublet"
    }else{
      type <- "singlet"
    }
  }else{
    if(sum(test$Outlier) == 1 & test$sample[1] == bestsample){
      type <- "singlet"
    }
    if(abs(plot_data$likelihood[2] - plot_data$likelihood[1]) > 2 *(abs(plot_data$likelihood[2] - plot_data$likelihood[nrow(plot_data)]) / (nrow(plot_data)-1)) ){
      type <- "singlet"
    }
  }
  
  
  if(show){
    
    plot <- ggplot(data=plot_data, aes(x = index, y=likelihood, label = sample)) +
      geom_point() +
      geom_line() + 
      geom_text_repel() + 
      ggtitle(paste0(barcode,"---- best Sample: ", bestsample, "  type: ", type))
    print(plot)
    
  }
  return(c(bestsample, type))
}

#singlet
assign_cell(doubletdata, doublet_small, "CATGAGTAGGCCTGCT")

#doublet
assign_cell(doubletdata,doublet_small, "AAAGAACTCGCGATCG")

#noise
assign_cell(doubletdata,doublet_small, "AAAGGGCGTTGCTGAT")


pdf("doublet_plots_own_method2.pdf", width = 8, height = 6)
for(cell in unique(doubletdata$cell)[1:500]){
  assign_cell(doubletdata,doublet_small, cell)
  
}
dev.off()

result <- data.frame()
i <- 1
for(cell in unique(doubletdata$cell)){
  print(i)
  res <- assign_cell(doubletdata,doublet_small, cell, show = F)
  result[cell,"sample"] <- res[1]
  result[cell, "type"] <- res[2]
  i <- i + 1
}

write.csv(result, "assigned/result.csv")
table(result$sample)
table(result$type)

data <- as.data.frame(table(result$sample))
colnames(data) <- c("sample", "#cells")

plot <- ggplot(data, aes(x=sample, y=.data[["#cells"]], fill=sample)) +
  geom_bar(stat="identity", width=1) 
plot

result <- read.csv("assigned/result.csv")
