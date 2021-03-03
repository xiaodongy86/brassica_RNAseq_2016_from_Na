
library(ggplot2)
library(ggthemes)
setwd("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq")
unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list <- read.csv("unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list_with_anotation_7_11_2020.csv")

####  up regulated gene networks
#  carbohydrate_metabolic
#  development
#  hormone_metabolic_and_response
#  senescence_lite
#  carbohydrate_transport

up_regualted_carbohydrate_metabolic <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/up_regualted_carbohydrate_metabolic.csv")
up_regualted_carbohydrate_metabolic_full <- data.frame(up_regualted_carbohydrate_metabolic,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list[match(up_regualted_carbohydrate_metabolic$ID,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list$sseqid),])
up_regualted_carbohydrate_metabolic_lite <- unique(up_regualted_carbohydrate_metabolic_full[,c(4,1,11,21,2)])
write.csv(up_regualted_carbohydrate_metabolic_lite, file = "up_regualted_carbohydrate_metabolic_lite.csv")
write.csv(up_regualted_carbohydrate_metabolic_full, file = "up_regualted_carbohydrate_metabolic_full.csv")


up_regualted_carbohydrate_transport <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/up_regualted_carbohydrate_transport.csv")
up_regualted_carbohydrate_transport_full <- data.frame(up_regualted_carbohydrate_transport,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list[match(up_regualted_carbohydrate_transport$ID,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list$sseqid),])
up_regualted_carbohydrate_transport_lite <- unique(up_regualted_carbohydrate_transport_full[,c(4,1,11,21,2)])
write.csv(up_regualted_carbohydrate_transport_lite, file = "up_regualted_carbohydrate_transport_lite.csv")
write.csv(up_regualted_carbohydrate_transport_full, file = "up_regualted_carbohydrate_transport_full.csv")


up_regualted_development <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/up_regualted_development.csv")
up_regualted_development_full <- data.frame(up_regualted_development,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list[match(up_regualted_development$ID,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list$sseqid),])
up_regualted_development_lite <- unique(up_regualted_development_full[,c(4,1,11,21,2)])
write.csv(up_regualted_development_lite, file = "up_regualted_development_lite.csv")
write.csv(up_regualted_development_full, file = "up_regualted_development_full.csv")



up_regualted_hormone_metabolic_and_response <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/up_regualted_hormone_metabolic_and_response.csv",header = T)
up_regualted_hormone_metabolic_and_response_full <- data.frame(up_regualted_hormone_metabolic_and_response,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list[match(up_regualted_hormone_metabolic_and_response$ID,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list$sseqid),])
up_regualted_hormone_metabolic_and_response_lite <- unique(up_regualted_hormone_metabolic_and_response_full[,c(4,1,11,21,2)])
write.csv(up_regualted_hormone_metabolic_and_response_lite, file = "./figures/up_regualted_hormone_metabolic_and_response_lite.csv")
write.csv(up_regualted_hormone_metabolic_and_response_full, file = "./figures/up_regualted_hormone_metabolic_and_response_full.csv")



up_regualted_senescence <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/up_regualted_senescence.csv")
up_regualted_senescence_full <- data.frame(up_regualted_senescence,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list[match(up_regualted_senescence$ID,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list$sseqid),])
up_regualted_senescence_lite <- unique(up_regualted_senescence_full[,c(4,1,11,21,2)])
write.csv(up_regualted_senescence_lite, file = "./figures/up_regualted_senescence_lite.csv")
write.csv(up_regualted_senescence_full, file = "./figures/up_regualted_senescence_full.csv")

###########################################################################################################################
####                                               down regulated gene networks
###########################################################################################################################


# glucosinolate
downregulated_glucosinolate <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_glucosinolate.csv")
downregulated_glucosinolate_full <- data.frame(downregulated_glucosinolate,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list[match(downregulated_glucosinolate$ID,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list$sseqid),])
downregulated_glucosinolate_lite <- unique(downregulated_glucosinolate_full[,c(4,1,11,21,2)])
downregulated_glucosinolate_lite$category <- rep(c("glucosinolate_metabolic"),nrow(downregulated_glucosinolate_lite))
downregulated_glucosinolate_lite <- na.omit(downregulated_glucosinolate_lite)
write.csv(downregulated_glucosinolate_lite,file = "/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_glucosinolate_lite.csv")
write.csv(downregulated_glucosinolate_full, file = "downregulated_glucosinolate_full.csv")

### response_to_insect
downregulated_response_to_insect <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_response_to_insect.csv")
downregulated_response_to_insect_full <- data.frame(downregulated_response_to_insect,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list[match(downregulated_response_to_insect$ID,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list$sseqid),])
downregulated_response_to_insect_lite <- unique(downregulated_response_to_insect_full[,c(4,1,11,21,2)])
downregulated_response_to_insect_lite$category <- rep(c("Response_to_insect"),nrow(downregulated_response_to_insect_lite))
#write.csv(downregulated_glucosinolate_lite,file = "/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_glucosinolate_lite.csv")
downregulated_response_to_insect_lite <- na.omit(downregulated_response_to_insect_lite)
write.csv(downregulated_response_to_insect_lite, file = "downregulated_response_to_insect_lite.csv")
write.csv(downregulated_response_to_insect_full, file = "downregulated_response_to_insect_full.csv")


### auxin_biosynthetic
downregulated_auxin_biosynthetic <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_auxin_biosynthetic.csv")
downregulated_auxin_biosynthetic_full <- data.frame(downregulated_auxin_biosynthetic,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list[match(downregulated_auxin_biosynthetic$ID,unique_match_DEGs_DEGs_Bj_CMS_NP_vs_P_complete_list$sseqid),])
downregulated_auxin_biosynthetic_lite <- unique(downregulated_auxin_biosynthetic_full[,c(4,1,11,21,2)])
downregulated_auxin_biosynthetic_lite$category <- rep(c("Auxin_biosynthetic"),nrow(downregulated_auxin_biosynthetic_lite))
#write.csv(downregulated_glucosinolate_lite,file = "/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_glucosinolate_lite.csv")
downregulated_auxin_biosynthetic_lite <- na.omit(downregulated_auxin_biosynthetic_lite)
write.csv(downregulated_auxin_biosynthetic_lite, file = "downregulated_auxin_biosynthetic_lite.csv")
write.csv(downregulated_auxin_biosynthetic_full, file = "downregulated_auxin_biosynthetic_full.csv")


#######  flower_development
downregulated_flower_development <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_flower_development.csv")
downregulated_flower_development$category <- rep(c("Flower_development"),nrow(downregulated_flower_development))



library(ggplot2)
library(viridis)
library(ggthemes)


up_regualted_carbohydrate_metabolic_P <- ggplot(up_regualted_carbohydrate_metabolic_lite, aes(x = category, y = Name, fill = logFC)) +
  #facet_wrap(BJ_HeatMap$category,scales = "free")+
  geom_tile(color = "white", size =0.1) +
  scale_fill_gradient2(name = "", low = "white", high = "#ff4040", midpoint = 0.5,
                       breaks=seq(0,4,1), #breaks in the scale bar
                       limits=c(0,4)) +
  #scale_fill_viridis(name="Log2FC") +
  labs(x=NULL, y=NULL, title= "") +
  geom_text(aes(x=category, y=Name, label=sprintf("%0.1f",logFC,digits=2)),color="black",size=6)+
  theme_tufte(base_family = "Helvetica", ticks = TRUE)+
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=14),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=14))
ggsave(up_regualted_carbohydrate_metabolic_P, width=3, height=12, file="/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/up_regualted_carbohydrate_metabolic_P.tiff", dpi = 600)



upregulated_flower_development <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/upregulated_flower_development.csv")
upregulated_flower_development$category <- rep(c("Flower_development"),nrow(upregulated_flower_development))
upregulated_flower_development <- na.omit(upregulated_flower_development)

upregulated_flower_development_P <- ggplot(upregulated_flower_development, aes(x = category, y = Name, fill = logFC)) +
  #facet_wrap(BJ_HeatMap$category,scales = "free")+
  geom_tile(color = "white", size =0.1) +
  scale_fill_gradient2(name = "", low = "white", high = "#ff4040", midpoint = 0.5,
                       breaks=seq(0,4,1), #breaks in the scale bar
                       limits=c(0,4)) +
  #scale_fill_viridis(name="Log2FC") +
  labs(x=NULL, y=NULL, title= "") +
  geom_text(aes(x=category, y=Name, label=sprintf("%0.1f",logFC,digits=2)),color="black",size=6)+
  theme_tufte(base_family = "Helvetica", ticks = TRUE)+
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=14),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=14))
ggsave(upregulated_flower_development_P, width=2.8, height=10, file="/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/upregulated_flower_development_P.tiff", dpi = 600)


up_regualted_carbohydrate_transport <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/up_regulated_carbohydrate_transport.csv")
up_regualted_carbohydrate_transport$category <- rep(c("Carbohydrate_transport"),nrow(up_regualted_carbohydrate_transport))
#up_regualted_carbohydrate_transport <- na.omit(up_regualted_carbohydrate_transport)

up_regualted_carbohydrate_transport_P <- ggplot(up_regualted_carbohydrate_transport, aes(x = category, y = Name, fill = logFC)) +
  #facet_wrap(BJ_HeatMap$category,scales = "free")+
  geom_tile(color = "white", size =0.1) +
  scale_fill_gradient2(name = "", low = "white", high = "#ff4040", midpoint = 0.5,
                       breaks=seq(0,4,1), #breaks in the scale bar
                       limits=c(0,4)) +
  #scale_fill_viridis(name="Log2FC") +
  labs(x=NULL, y=NULL, title= "") +
  geom_text(aes(x=category, y=Name, label=sprintf("%0.1f",logFC,digits=2)),color="black",size=6)+
  theme_tufte(base_family = "Helvetica", ticks = TRUE)+
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=14),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=14))
ggsave(up_regualted_carbohydrate_transport_P, width=3, height=10, file="/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/up_regualted_carbohydrate_transport_P.tiff", dpi = 600)




#up_regualted_senescence_lite
up_regualted_senescence_lite$category <- rep(c("Senescence"),nrow(up_regualted_senescence_lite))
#up_regualted_carbohydrate_transport <- na.omit(up_regualted_carbohydrate_transport)

up_regualted_senescence_lite_P <- ggplot(up_regualted_senescence_lite, aes(x = category, y = Name, fill = logFC)) +
  #facet_wrap(BJ_HeatMap$category,scales = "free")+
  geom_tile(color = "white", size =0.1) +
  scale_fill_gradient2(name = "", low = "white", high = "#ff4040", midpoint = 0.5,
                       breaks=seq(0,4,1), #breaks in the scale bar
                       limits=c(0,4)) +
  #scale_fill_viridis(name="Log2FC") +
  labs(x=NULL, y=NULL, title= "") +
  geom_text(aes(x=category, y=Name, label=sprintf("%0.1f",logFC,digits=2)),color="black",size=7)+
  theme_tufte(base_family = "Helvetica", ticks = TRUE)+
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=14),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=17))
ggsave(up_regualted_senescence_lite_P, width=3, height=10, file="/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/up_regualted_senescence_lite_P.tiff", dpi = 600)


#up_regualted_hormone_metabolic_and_response_lite

up_regualted_hormone_metabolic_and_response_lite$category <- rep(c("hormone_metabolic_and_response"),nrow(up_regualted_hormone_metabolic_and_response_lite))
up_regualted_hormone_metabolic_and_response_lite <- na.omit(up_regualted_hormone_metabolic_and_response_lite)

up_regualted_hormone_metabolic_and_response_lite_P <- ggplot(up_regualted_hormone_metabolic_and_response_lite, aes(x = category, y = Name, fill = logFC)) +
  #facet_wrap(BJ_HeatMap$category,scales = "free")+
  geom_tile(color = "white", size =0.1) +
  scale_fill_gradient2(name = "", low = "white", high = "#ff4040", midpoint = 0.5,
                       breaks=seq(0,4,1), #breaks in the scale bar
                       limits=c(0,4)) +
  #scale_fill_viridis(name="Log2FC") +
  labs(x=NULL, y=NULL, title= "") +
  geom_text(aes(x=category, y=Name, label=sprintf("%0.1f",logFC,digits=2)),color="black",size=6)+
  theme_tufte(base_family = "Helvetica", ticks = TRUE)+
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=14),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=14))
ggsave(up_regualted_hormone_metabolic_and_response_lite_P, width=3, height=10, file="/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/up_regualted_hormone_metabolic_and_response_lite_P.tiff", dpi = 600)




#   downregulated_glucosinolate_lite
downregulated_glucosinolate_lite$category <- rep(c("glucosinolate metabolic"),nrow(downregulated_glucosinolate_lite))
#write.csv(downregulated_glucosinolate_lite,file = "/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_glucosinolate_lite.csv")
downregulated_glucosinolate_lite <- na.omit(downregulated_glucosinolate_lite)
write.csv(downregulated_glucosinolate_lite,file = "/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_glucosinolate_lite.csv")

downregulated_glucosinolate_P <- ggplot(downregulated_glucosinolate_lite, aes(x = category, y = Name, fill = logFC)) +
  #facet_wrap(BJ_HeatMap$category,scales = "free")+
  geom_tile(color = "white", size =0.1) +
  scale_fill_gradient2(name = "", low = "#00c2c7", high = "white", midpoint = -0.5,
                       limits=c(-2.5,0)) +
  #scale_fill_viridis(name="Log2FC") +
  labs(x=NULL, y=NULL, title= "") +
  geom_text(aes(x=category, y=Name, label=sprintf("%0.1f",logFC,digits=2)),color="black",size=6)+
  theme_tufte(base_family = "Helvetica", ticks = TRUE)+
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=14),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=14))
ggsave(downregulated_glucosinolate_P, width=3.5, height=12, file="/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_glucosinolate_P.tiff", dpi = 600)


#downregulated_response_to_insect_lite

downregulated_response_to_insect_P <- ggplot(downregulated_response_to_insect_lite, aes(x = category, y = Name, fill = logFC)) +
  #facet_wrap(BJ_HeatMap$category,scales = "free")+
  geom_tile(color = "white", size =0.1) +
  scale_fill_gradient2(name = "", low = "#00c2c7", high = "white", midpoint = -0.5,
                       limits=c(-2.5,0))+
  #scale_fill_viridis(name="Log2FC") +
  labs(x=NULL, y=NULL, title= "") +
  geom_text(aes(x=category, y=Name, label=sprintf("%0.1f",logFC,digits=2)),color="black",size=6)+
  theme_tufte(base_family = "Helvetica", ticks = TRUE)+
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=14),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=14))
ggsave(downregulated_response_to_insect_P, width=3.5, height=7, file="/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_response_to_insect_P.tiff", dpi = 600)

# downregulated_auxin_biosynthetic_lite
downregulated_auxin_biosynthetic_P <- ggplot(downregulated_auxin_biosynthetic_lite, aes(x = category, y = Name, fill = logFC)) +
  #facet_wrap(BJ_HeatMap$category,scales = "free")+
  geom_tile(color = "white", size =0.1) +
  scale_fill_gradient2(name = "", low = "#00c2c7", high = "white", midpoint = -0.5,
                       limits=c(-2.5,0))+
  #scale_fill_viridis(name="Log2FC") +
  labs(x=NULL, y=NULL, title= "") +
  geom_text(aes(x=category, y=Name, label=sprintf("%0.1f",logFC,digits=2)),color="black",size=6)+
  theme_tufte(base_family = "Helvetica", ticks = TRUE)+
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=14),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=14))
ggsave(downregulated_auxin_biosynthetic_P, width=3.5, height=4, file="/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_auxin_biosynthetic_P.tiff", dpi = 600)


#downregulated_flower_development

downregulated_flower_development$Name <- c("LSU4","FD","YAB1","CYP78A5","CYP78A5.","GA20OX1","FIE")
downregulated_flower_development_P <- ggplot(downregulated_flower_development, aes(x = category, y = Name, fill = logFC)) +
  #facet_wrap(BJ_HeatMap$category,scales = "free")+
  geom_tile(color = "white", size =0.1) +
  scale_fill_gradient2(name = "", low = "#00c2c7", high = "white", midpoint = -0.5,
                       limits=c(-2.5,0))+
  #scale_fill_viridis(name="Log2FC") +
  labs(x=NULL, y=NULL, title= "") +
  geom_text(aes(x=category, y=Name, label=sprintf("%0.1f",logFC,digits=2)),color="black",size=6)+
  theme_tufte(base_family = "Helvetica", ticks = TRUE)+
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=14),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=14))
ggsave(downregulated_flower_development_P, width=3.5, height=6, file="/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/downregulated_flower_development_P.tiff", dpi = 600)



down_regulated_sucrose_biosynthetic <- read.csv("/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/down_regulated_sucrose_biosynthetic.csv")
down_regulated_sucrose_biosynthetic$Name <- c("SUS6", "SUS6.", "SUS5")
down_regulated_sucrose_biosynthetic_P <- ggplot(down_regulated_sucrose_biosynthetic, aes(x = category, y = Name, fill = logFC)) +
  #facet_wrap(BJ_HeatMap$category,scales = "free")+
  geom_tile(color = "white", size =0.1) +
  scale_fill_gradient2(name = "", low = "#00c2c7", high = "white", midpoint = -0.5,
                       limits=c(-2.5,0))+
  #scale_fill_viridis(name="Log2FC") +
  labs(x=NULL, y=NULL, title= "") +
  geom_text(aes(x=category, y=Name, label=sprintf("%0.1f",logFC,digits=2)),color="black",size=6)+
  theme_tufte(base_family = "Helvetica", ticks = TRUE)+
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=14),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=14))
ggsave(down_regulated_sucrose_biosynthetic_P, width=3, height=3, file="/Users/mackenzie/RProject/Brassica_juncea/2016RNA-seq/figures/down_regulated_sucrose_biosynthetic_P.tiff", dpi = 600)

