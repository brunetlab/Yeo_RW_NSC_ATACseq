# Figure 3A Marker tSNE
# Matthew Buckley

rm(list=ls())
library(Seurat)
library(ggplot2)
library(tidyverse)


# Change directory for testing purposes. Subsequent paths are relative.
setwd("/Users/robinyeo/Dropbox/RWY_ATAC_Code_Checking/Code_Checking_resubmission/")

# Load Data
load(file = "Original_Data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata
features <- colnames(d)

GO <- read.delim(file="Original_Data/GO_Cell_Adhesion_GO0007155.txt", header=TRUE)
genelist <- tolower(GO$Gene.name)


colnames(d) <- tolower(colnames(d))
adhesion_data <- d[, colnames(d) %in% genelist]
adhesion_data$adhesion_response <- rowSums(adhesion_data)
meta <- d[, colnames(d) %in% c("age", "replicate", "celltype", "tsne_1", "tsne_2")]
d2 <- cbind(meta, adhesion_data)
d2$age <- factor(d2$age, levels = c("y", "o"), ordered=TRUE)


low_col <- "grey"
high_col <- "darkred"
size_range <- c(.8, 1.2)
alpha_range <- c(0.10, .90)

# Saturate color at 99.5% percentile of response
threshold <- quantile(d2$adhesion_response, .99)
adhesion_sig_cap <- d2$adhesion_response
adhesion_sig_cap[adhesion_sig_cap > threshold] <- threshold
d2$adhesion_sig_cap <- adhesion_sig_cap

q <- ggplot(data = filter(d2, age == "y"), aes(x = tsne_1, y = tsne_2,
			size = adhesion_sig_cap, alpha = adhesion_sig_cap, color = adhesion_sig_cap))
q <- q + geom_point() + scale_size(range = size_range) + scale_alpha(range = alpha_range)
q <- q + scale_colour_gradient(low = low_col, high = high_col)
q <- q + theme(legend.position="none", panel.background = element_blank(),axis.line = element_line(colour = "black"), axis.text=element_text(colour="black", size=rel(2)))
q
ggsave("Output_Figs/YOUNG_tSNE_GO_Cell_Adhesion_GO0007155.png", q, height = 7, width = 6)

q <- ggplot(data = filter(d2, age == "o"), aes(x = tsne_1, y = tsne_2,
                                               size = adhesion_sig_cap, alpha = adhesion_sig_cap, color = adhesion_sig_cap))
q <- q + geom_point() + scale_size(range = size_range) + scale_alpha(range = alpha_range)
q <- q + scale_colour_gradient(low = low_col, high = high_col)
q <- q + theme(legend.position="none", panel.background = element_blank(),axis.line = element_line(colour = "black"), axis.text=element_text(colour="black", size=rel(2)))
q
ggsave("Output_Figs/OLD_tSNE_GO_Cell_Adhesion_GO0007155.png", q, height = 7, width = 6)

#######################
#Plotting cell cycle states here. Cleaned up code by Paloma N.

rm(list=ls())
library(Seurat)
library(ggplot2)
library(tidyverse)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Dropbox/Code_Checking_resubmission/")

load(file = "Original_Data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata
features <- colnames(d)
colnames(d) <- tolower(colnames(d))
meta <- d[, colnames(d) %in% c("age", "replicate", "celltype", "tsne_1", "tsne_2", "phase")]
meta <- filter(meta, age == "y")

q <- ggplot(meta, aes(x = tsne_1, y = tsne_2, color = phase))
q <- q + geom_point() 
q <- q + scale_colour_manual(values=c("dodgerblue","orange","red"))
q <- q + theme(legend.position="none", panel.background = element_blank(),axis.line = element_line(colour = "black"), axis.text=element_text(colour="black", size=rel(2)))
q
ggsave("Output_Figs/YOUNG_tSNE_cell_cycle_Stages.png", q, height = 7, width = 6)

#Now for old
load(file = "Original_Data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata
features <- colnames(d)
colnames(d) <- tolower(colnames(d))
meta <- d[, colnames(d) %in% c("age", "replicate", "celltype", "tsne_1", "tsne_2", "phase")]
meta <- filter(meta, age == "o")

q <- ggplot(meta, aes(x = tsne_1, y = tsne_2, color = phase))
q <- q + geom_point() 
q <- q + scale_colour_manual(values=c("dodgerblue","orange","red"))
q <- q + theme(legend.position="none", panel.background = element_blank(),axis.line = element_line(colour = "black"), axis.text=element_text(colour="black", size=rel(2)))
q
ggsave("Output_Figs/OLD_tSNE_cell_cycle_Stages.png", q, height = 7, width = 6)



