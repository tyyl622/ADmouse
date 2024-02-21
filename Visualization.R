library(patchwork)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(tidyr)
library(stringr)
library(stringi)
load('allmarkerfiltRegion.rda')
load('brain_slice1_final.rda')
load('brain_slice1meta.rda')
regionorder <- c('Layer1/2','Layer3/4','Layer5','Layer6','CA1','GrDG','Dendritic','CA2/3',
                  'C0','C1','C2','C3','C4','C5')
regionorder1 <- c('C5','C4','C3','C2','C1','C0','CA2/3','Dendritic',
                'GrDG','CA1','Layer6','Layer5','Layer3/4','Layer1/2')
Idents(brain_slice1_final) <- brain_slice1_final$region
levels(Idents(brain_slice1_final)) <- regionorder1
markplot3 <- lapply(unique(allmarkerfilt$cluster),function(x){
  tmp <- allmarkerfilt[allmarkerfilt$cluster==x,] 
  tmp <- tmp[order(tmp$avg_log2FC,decreasing=T),]
  tmp <- tmp$gene %>% .[1:5]
  return(tmp)
}) %>% unlist(.) %>% unique(.)
mycol <- colorRampPalette(brewer.pal(11, "RdYlBu"))
p2 <- DotPlot(object = brain_slice1_final, features =markplot3,assay='SCT')+
  scale_colour_gradient2(low="#4575B4",mid="white",high="#D73027")+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=15))
ggsave(p2,filename='allmarkerDotTop10.pdf',width=25,height=7)
ggsave(p2,filename='allmarkerDotTop10.tiff',width=25,height=7,dpi=300)

markplot5 <- data.frame()
for (x in regionorder){
  tmp <- allmarkerfilt[allmarkerfilt$cluster==x,] 
  tmp <- tmp[order(tmp$avg_log2FC,decreasing=T),]
  tmp <- tmp[1:10,]
  markplot5 <- rbind(markplot5,tmp)
}
write.table(markplot5,file='markplot5.txt',quote=F,row.names=F,sep='\t')

#complexheatmap
library(ComplexHeatmap)
group.use <- c()
for (re in regionorder){
  tmp <- brain_slice1meta[brain_slice1meta$region==re,'Spot']
  group.use <- c(group.use,tmp)
}

annogene <-  lapply(unique(allmarkerfilt$cluster),function(x){
  tmp <- allmarkerfilt[allmarkerfilt$cluster==x,] 
  tmp <- tmp[order(tmp$avg_log2FC,decreasing=T),]
  tmp <- tmp$gene %>% .[1:2]
  return(tmp)
}) %>% unlist(.) %>% unique(.)

annogene <- unique(c(annogene, 'Calb1','Cux2','Bcl11b','Th','Gad2','Map2k1','Gabra1','Gabra5','Foxo1','Foxp2','Satb2','Tbr1'))
#手动把需要标注的'Calb1','Cux2','Bcl11b'按照其主要归属脑区的所在顺序，加入到markplot3(Newmarplot3.txt)，再读进来
newmarkplot3 <- read.table('Newmarplot3.txt',sep='\t',header=F,check.names=F)
newmarkplot3 <- unique(newmarkplot3$V1)
heatbox <- brain_slice1_final@assays$SCT@scale.data[newmarkplot3,group.use] %>% data.frame(.,check.names=F)
ann <- data.frame(brain_slice1meta[colnames(heatbox),'region'])
colnames(ann) <- 'Region'
colours <- list('Region' = c('Layer1/2'='#98D277','Layer3/4'='#C3AAD2','Layer5'='LemonChiffon',
                              'Layer6'='#A6CEE3','CA1'='Cyan','GrDG'='SlateBlue1','Dendritic'='GreenYellow',
                              'CA2/3'='Thistle1','C0' = 'Red', 'C1' = '#438EC0','C2'='#3BA432','C3'='#63A8A0',
                               'C4'='Gold','C5'='#ED8F47'))
colAnn  <- ComplexHeatmap::HeatmapAnnotation(df = ann,
  which = 'column',
  col = colours)
theatbox <-  t(heatbox) %>% pheatmap:::scale_rows(.)
theatbox <- t(theatbox)
colors2 <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
colors <- colorRampPalette(colors =colors2 )(100)
n_col = length(colors)
lim = max(abs(theatbox), na.rm = TRUE)
library(colorRamp2)
colors2 <- colorRamp2(seq(-lim, lim, length = n_col), colors)
colsplit <- factor(ann$Region,levels=regionorder,labels=c(1:14))
p <- ComplexHeatmap::Heatmap(theatbox,show_column_names=FALSE,cluster_columns = FALSE,cluster_rows = FALSE,
                        show_row_dend = FALSE,show_column_dend = FALSE,
                        show_row_names = FALSE,
                        col = colors2, 
                        border=TRUE,border_gp = gpar(col = "white", lty = 1,lwd=0.5),
                        top_annotation=colAnn,
                        column_split=colsplit,
                        column_title = NULL,
                        )
pos <-  which(rownames(heatbox) %in% annogene)
p <- p + rowAnnotation(link = anno_mark(at = pos, 
    labels = rownames(heatbox)[pos], labels_gp = gpar(fontsize = 25)))
pdf('allmarkerComplexHeatTop10.pdf',width=12,height=25)
# png('allmarkerComplexHeatDotTop10.png',width=1000,height=2000)
p
dev.off()

#marker umap or tsne#################################################
DefaultAssay(brain_slice1_final) <- 'SCT'
allmarkerfilt$pct.diff <- allmarkerfilt$pct.1-allmarkerfilt$pct.2
features1 <-  lapply(c('Layer1/2','Layer3/4','Layer5','Layer6','CA1','GrDG','Dendritic','CA2/3'),function(x){
  tmp <- allmarkerfilt[allmarkerfilt$cluster==x,] 
  tmp <- tmp[order(tmp$pct.diff,decreasing=T),]
  tmp <- tmp$gene %>% .[1]
  return(tmp)
}) %>% unlist(.) %>% unique(.)
features2 <-  lapply(c('C0','C1','C2','C3','C4','C5'),function(x){
  tmp <- allmarkerfilt[allmarkerfilt$cluster==x,] 
  tmp <- tmp[order(tmp$pct.diff,decreasing=T),]
  tmp <- tmp$gene %>% .[1]
  return(tmp)
}) %>% unlist(.) %>% unique(.)
features <-  unique(c(features1,'Satb2','Bcl11b',features2,'Foxp2'))
fix.sc <- scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p1 <- FeaturePlot(object = brain_slice1_final, reduction = "tsne",
                  features = features,combine = FALSE)
p2 <- lapply(p1,function(x)x+fix.sc) 
p3 <- aplot::plot_list(gglist=p2,ncol=4)
ggsave(p3,filename='markertsne.pdf',width=13,height=11)

###feaplot######################################################
DefaultAssay(brain_slice1_final) <- 'SCT'
mySpatial <- function(mat,assay='SCT',fea,threshold=0.95,savetype='tiff',filename,dpi=300){
  DefaultAssay(mat) <- assay
    for ( i in fea){
    test <- mat[i,]
    pp <- quantile(test@assays$SCT@data[1,],threshold)
    pp <- signif(as.numeric(pp),2)
    pp1 <- min(test@assays$SCT@data)
    mid <- (pp+pp1)/2
    test@assays$SCT@data[1,] <- sapply(test@assays$SCT@data[1,], function(x) {ifelse(x>pp,pp,x)})
    rownames(test@assays$SCT@data) <- paste0(i,'   ')
    p1 <- SpatialPlot(test, features = rownames(test@assays$SCT@data),combine=FALSE,
                      stroke = NA,slot = 'data',
                      crop=FALSE,pt.size.factor=1.1,
                      alpha = c(0.8,1))
    p2 <- lapply(p1,function(y){
      return(y+scale_fill_gradientn(colours=rev(brewer.pal(n = 11, name = "RdYlBu")),limits = c(signif(pp1,2),pp),breaks=c(signif(pp1,2),signif(mid,2),pp))+
                theme(legend.position='none',
                plot.margin=margin(0,0,0,-5),
                plot.title=element_blank()))})
    legend <- get_legend(p2[[1]] + guides(color = guide_legend(nrow = 1)) +
                        theme(legend.position = "top",
                        legend.justification='left',
                        legend.direction='horizontal',
                        legend.margin=margin(30,0,0,5),
                        legend.title=element_text(size=20),
                        legend.text=element_text(size=15)))
    plots <- align_plots(legend,plotlist=p2)
    bottom_row <- plot_grid(plotlist=plots[2:length(plots)],nrow=1)
    p4 <-  plot_grid(legend,NULL,bottom_row,nrow=3,rel_widths=c(0.4,0,1),rel_heights = c(0.5,0.05,1))
    width <- 6
    height <- 2.5
    if (savetype=='tiff'){
        ggsave(p4,filename = paste0(i,'_',filename,'.',savetype),width=width,height=height,dpi=dpi, bg = 'white')
    } else if (savetype=='pdf') {
       ggsave(p4,filename = paste0(i,'_',filename,'.',savetype),width=width,height=height, bg = 'white')
    } else {
      ggsave(p4,filename = paste0(i,'_',filename,'.tiff'),width=width,height=height,dpi=dpi, bg = 'white')
      ggsave(p4,filename = paste0(i,'_',filename,'.pdf'),width=width,height=height, bg = 'white')
    }
  }
}
mySpatial(brain_slice1_final,fea=c('Slc2a2','Slc2a4','Slc2a5','Slc2a6','Slc2a8','Slc2a9','Slc2a10','Slc2a12'),threshold=1,filename='plotSlice1',savetype='pdf')

#violin plot############################################################
fea <- c('Map2k1','Eloc','Rbx1','Pfkm','Elob','Mapk1','Pgk1','Aldoa','Rheb','Hras','Prkacb','Atp2a2','Thra','Atp1b1','Actg1') 
region <- 'Dendritic'
x.pos <- 0
feainfo <- data.frame()
dbox <- data.frame()
for (i in fea){
  x.pos <- x.pos+1
  gene <- i
  dbox1 <- brain_slice1exp[gene,brain_slice1meta[brain_slice1meta$region2==region,'Spot']] %>% data.frame(.,check.names=F)
  colnames(dbox1)[1] <- 'expression'
  dbox1$gene <- gene
  dbox1$Type <- brain_slice1meta[rownames(dbox1),'Type']
  print(min(dbox1$expression))
  #prepare stat.test
  stat.test <- DEGpadj[DEGpadj$gene==i & DEGpadj$region==region,]
  if(nrow(stat.test)==4){
      stat.test <- stat.test[-grep('SLvsWT_24M',stat.test$versus),]
      stat.test <- cbind(data.frame(`.y.`=rep('len',3)),stat.test,data.frame(method='Wilcox'))
      stat.test <- stat.test[order(stat.test$versus,decreasing=T),]
      stat.test$group1 <- c('24WT','24SL','24AD')
      stat.test$group2 <- c('4WT','24AD','24WT')
      stat.test <- stat.test[,-c(2,3,4,6,7)]
      stat.test$p.format <- as.character(format(stat.test$pvalue,scientific = TRUE))
      stat.test$p.signif <- ifelse(stat.test$padj<0.05,
                                    ifelse(stat.test$padj<0.01,
                                      ifelse(stat.test$padj<0.001,'***','**'),'*'),'ns')
      colnames(stat.test)[c(2,3)] <- c('p','p.adj')
      stat.test <- stat.test[,c(1,5,6,2,3,7,8)]
      max <- max(dbox1$expression)
      stat.test <- stat.test %>%
                    mutate(y.position = c(max+0.5,max+1.5,max+1))
      stat.test$gene <- i
      stat.test$xpos <- c(x.pos-0.2,x.pos+0.2,x.pos)

  } else{
    stat.test <- data.frame()
  }
    feainfo <- rbind(feainfo,stat.test)
    dbox <- rbind(dbox,dbox1)
}
dbox$Type <- factor(dbox$Type,levels=c('4WT','24WT','24AD','24SL'))
p <- ggviolin(dbox, x = "gene", y = "expression", fill = "Type",
        palette = c("gray", "black","#de1e17",'blue'),
        add = "boxplot",width=1,add.params = list(color='white'),
        title="",font.label=list(size=35),nrow=1,trim=TRUE)+theme_bw()+
theme(#axis.text.x=element_text(size=35,angle=45,hjust=1,vjust=1),
      axis.text.x=element_text(size=35),
      axis.text.y=element_text(size=30),
      axis.title=element_text(size=35),
      legend.position="none",panel.grid = element_blank())+
      coord_cartesian(ylim=c(min(dbox$expression)-0.5, max(dbox$expression)+2))+
ylab('Normalized expression')+xlab("") +
stat_pvalue_manual(feainfo, label = "{p.signif}",size=10,x='xpos') 
ggsave(p,filename = 'Fig6violinDendritic.pdf',width=30,height=8)

#single violin#
group  <- list(c('24SL','24AD'),
              c('24WT','4WT'),
              c('24AD','24WT'))
#Dendritic
#add pvalue manually
library(grid)
fea <- c('Gpx4','Slc25a4','Mif','Coq7','Surf1','Txn1','Sod2')
for (i in fea){
  gene <- i
  region <- 'Dendritic'
  dbox1 <- brain_slice1exp[gene,brain_slice1meta[brain_slice1meta$region2==region,'Spot']] %>% data.frame(.,check.names=F)
  colnames(dbox1)[1] <- 'expression'
  dbox1$gene <- gene
  dbox1$Type <- brain_slice1meta[rownames(dbox1),'Type']
  print(min(dbox1$expression))
  #prepare stat.test
  stat.test <- DEGpadj[DEGpadj$gene==i & DEGpadj$region==region,]
  stat.test <- stat.test[-grep('24SLvs24WT',stat.test$versus),]
  if(nrow(stat.test)==3){
      stat.test <- cbind(data.frame(`.y.`=rep('len',3)),stat.test,data.frame(method='Wilcox'))
      stat.test <- stat.test[order(stat.test$versus,decreasing=T),]
      stat.test$group1 <- c('24WT','24SL','24AD')
      stat.test$group2 <- c('4WT','24AD','24WT')
      stat.test <- stat.test[,-c(2,3,4,6,7)]
      stat.test$p.format <- as.character(format(stat.test$pvalue,scientific = TRUE))
      stat.test$p.signif <- ifelse(stat.test$padj<0.05,
                                    ifelse(stat.test$padj<0.01,
                                      ifelse(stat.test$padj<0.001,'***','**'),'*'),'ns')
      colnames(stat.test)[c(2,3)] <- c('p','p.adj')
      stat.test <- stat.test[,c(1,5,6,2,3,7,8)]
      max <- max(dbox1$expression)
      stat.test <- stat.test %>%
                    mutate(y.position = c(max+0.5,max+1.5,max+1))
      p <- ggviolin(dbox1, x = "Type", y = "expression", fill = "Type",
              palette = c("gray", "black","#de1e17",'blue'),
              add = "boxplot", add.params = list(fill = "white"),
              title=paste0(region," ",gene),font.label=list(size=38),nrow=1,trim = TRUE)+theme_bw()+
        theme(axis.text.x=element_text(size=38,angle=45,hjust=1,vjust=1),
              axis.text.y=element_text(size=38),
              axis.title=element_text(size=38,face='italic'),
              title=element_text(size=38),
              legend.position="none",panel.grid = element_blank(),
              plot.margin = margin(0.5, 2, 0, 2, "cm"))+
              coord_cartesian(ylim=c(0, max(dbox1$expression)+1.5))+
        ylab('Normalized expression')+xlab("") +
        stat_pvalue_manual(stat.test, label = "{p.signif}",size=9)+labs(x=NULL)
  }
}

#enrichment visualize##############################################
#KEGG
#HPd
library(scales)
allenrich1 <- read.table('Dendritic_enrichKEGG_keygene.txt',header=T,check.names=F,sep='\t')
allenrich1$Description <- gsub(' \\- Mus musculus \\(house mouse\\)','',allenrich1$Description)
# allenrich2 <- read.table('HPd_enrichhall_keygene.txt',header=T,check.names=F,sep='\t')
# allenrich <- rbind(allenrich1,allenrich2)
allenrich  <- allenrich1
allenrich <- allenrich[order(allenrich$p.adjust),]
allenrich$`-log10(p.adjust)` <- (-log10(allenrich$p.adjust))
allenrich <- allenrich[,-8]#
allenrich$geneRation <- sapply(allenrich$GeneRatio,function(x){
  tmp <- as.numeric(unlist(str_split(x,'/')))
  return(tmp[1]/tmp[2])
})
boxdata <- allenrich
boxdata <- boxdata[order(boxdata$p.adjust,decreasing = T),]
#define row color
class <- read.delim2('hpdclass.txt',check.names=F,header=T,row.names=1)
rowcol <- class[boxdata$Description,'color']
boxdata$Description <- factor(boxdata$Description,levels = unique(boxdata$Description))
p1 <- ggplot(boxdata,aes(x=Count,y=Description)) +
  geom_bar(stat='identity',aes(fill=`-log10(p.adjust)`,colour=`-log10(p.adjust)`))+
  # scale_fill_gradient(low="Blue1",high="PaleGreen")+
  # scale_colour_gradient(low="Blue1",high="PaleGreen")+
  scale_fill_gradientn(colours=c("#3300CC","#3399FF","white","#FF3333","#CC0000"))+
  scale_colour_gradientn(colours=c("#3300CC","#3399FF","white","#FF3333","#CC0000"))+
  theme_bw()+xlab('Count')+ylab('')+
  theme(axis.title = element_text(size=30),
        axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=30,color=rowcol),
        legend.text = element_text(size=30),
        legend.title = element_text(size=30),
        title=element_text(size=30))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=40))+
  ggtitle('')
ggsave(p1,filename = "enrichKEGG_HPdkeyGenes.pdf",width = 12,height = 23)

#GO#####################
#Dendritic
boxdata <- read.delim2('Dendritic_enrichGO_keygene.txt',header = T,sep='\t',check.names = F)
boxdata$geneRatio <- as.vector(sapply(boxdata$GeneRatio,function(x){
  tmp1 <- as.numeric(str_split_fixed(x,'/',2)[[1]])
  tmp2 <- as.numeric(str_split_fixed(x,'/',2)[[2]])
  ratio <- round(tmp1/tmp2,2)
  return(ratio)
}))
boxdata <- boxdata[order(boxdata$geneRatio,decreasing=T),]
boxdata <- rbind(subset(boxdata,Category=='BP')[1:20,],
                 subset(boxdata,Category=='MF')[1:20,],
                 subset(boxdata,Category=='CC')[1:20,]) %>% na.omit(x)
boxdata <- boxdata[order(boxdata$geneRatio,decreasing = F),]
boxdata$Description <- factor(boxdata$Description,levels = boxdata$Description)
boxdata$LogAdjustPvalue <- as.vector(sapply(boxdata$p.adjust,function(x){-(log(as.numeric(x)))}))
pdf('Dendritic_enrichGO_keygene.pdf',width = 17,height =20)
# png('Dendritic_enrichGO_keygene.tiff',width = 1150,height =1000)
ggplot(boxdata,aes(x=geneRatio,y=Description,fill = -LogAdjustPvalue,color=-LogAdjustPvalue)) +
  geom_vline(xintercept = c(0.05,0.1,0.15,0.2),lty=2, lwd=1,col='Azure3')+
  geom_point(shape = 21,aes(size=geneRatio))+
  scale_fill_gradientn(colours=c("#3300CC","#3399FF","white","#FF3333","#CC0000"))+
   scale_color_gradientn(colours=c("#3300CC","#3399FF","white","#FF3333","#CC0000"))+
  scale_size_continuous(range = c(5,15))+
  theme_bw()+xlab('GeneRatio')+ylab('')+xlim(c(0,0.22))+ggtitle('')+
  theme(axis.title = element_text(size=25),
        axis.text.x = element_text(size=25,angle=45,vjust=1,hjust=1),
        axis.text.y = element_text(size=25),
        legend.text = element_text(size=25),
        legend.title = element_text(size=25),
        title=element_text(size=25))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))+
  facet_grid(Category~., scale='free',space='free')+
  theme(strip.text.y = element_text(size = 25))
dev.off()

#KEGGpathway-gene complex heatmap########################################
setwd('GSVA')
load('brain_slice1exp.rda')
load('brain_slice1meta.rda')
load('wilcox/resultHPF.rda')
#gsva pathway heatmap
library(reshape2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
gsvaint <- read.delim2('GSVAOut.txt',header=T,check.names=F)
rownames(gsvaint) <- gsvaint$id
load('gsvafilt_HPd.rda')
load('gsvafilt_HPs.rda')
hpdpath <- read.table('Dendritic_enrichKEGG_keygene.txt',sep='\t',header=T,check.names=F)
hpdpath$Description <- gsub(' \\- Mus musculus \\(house mouse\\)','',hpdpath$Description)
hpdgsva <- gsvafilt_hpf[gsvafilt_hpf$pathway %in% unique(hpdpath$Description),]
#Dendritic
region <- 'Dendritic'
mypathway <- hpdgsva
gsvaOut3 <- data.frame()
for (i in seq(1,nrow(mypathway),3)){
  for (type in c('24SL','24AD','24WT','4WT')){
    mypa <- mypathway[i,'pathway']
    reg <- mypathway[i,'region']
    gsva_tmp1 <- gsvaint[mypa,c('id',brain_slice1meta[brain_slice1meta$region2==reg & brain_slice1meta$Type==type,'Spot'])]
    gsva_ave <- data.frame(score = mean(as.numeric(gsva_tmp1[1,-1]))) %>% data.frame(.,check.names = F)
    gsva_ave$type <-  paste0(type,'_',reg)
    gsva_ave$pathway <- mypa
    gsvaOut3 <- rbind(gsvaOut3,gsva_ave)
    
  }
}
plotdata2 <- pivot_wider(gsvaOut3,names_from=pathway,values_from = score) %>% data.frame(.,check.names=F)
rownames(plotdata2) <- c('24 SL','24 AD','24 WT','4 WT')
plotdata2 <- plotdata2[,-1]
color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")

# #pathway color
class <- read.delim2('hpdclass.txt',check.names=F,header=T,row.names=1)#this file contain your pathway and corresponding color code
plotdata2 <- plotdata2[,rownames(class)[rownames(class) %in% colnames(plotdata2)]]
plotdata2 <- plotdata2[c(4,3,2,1),]
heatbox <- t(plotdata2) 

#GSVA score heatmap
colAnn <-  HeatmapAnnotation(
        text = anno_text(colnames(heatbox), rot = 90, 
        location = unit(1, "npc"), just = "right",
        gp = gpar(col = "black", fontsize = 30,fontface='bold')),
        annotation_height = max_text_width(colnames(heatbox))
    )
colsplit <- c(1,2,3,4)
rowsplit <- 1:nrow(heatbox)
# scaleheatbox <-  heatbox %>% pheatmap:::scale_mat(.,'column')
scaleheatbox <-  heatbox %>% pheatmap:::scale_rows(.)
scaleheatbox <- data.frame(scaleheatbox,check.names=F)
scaleheatbox <- scaleheatbox[order(scaleheatbox$`24 WT`,decreasing=T),]
rowcol <- class[rownames(scaleheatbox),'color']
colors2 <- c("RoyalBlue1","white","Brown1")
pdf('Somatic_gsvaheat.pdf',width=20,height=15)
p <- ComplexHeatmap::Heatmap(as.matrix(scaleheatbox),show_column_names=FALSE,cluster_columns = FALSE,cluster_rows = FALSE,
                        show_row_dend = FALSE,show_column_dend = FALSE,
                        show_row_names = TRUE,
                        border=TRUE,border_gp = gpar(col = "black",lwd=1.5,lty=1),#border
                        col = colorRampPalette(colors =colors2 )(200), bottom_annotation=colAnn,
                        width = unit(6, "cm"), height = unit(25, "cm"),
                        column_split=colsplit,
                        row_split=rowsplit,
                        row_title = NULL,
                        column_title=region,
                        column_title_gp = gpar(fontsize = 35, fontface = "bold"),
                        row_names_gp = gpar(fontsize = 30,col=rowcol),
                        heatmap_legend_param=list(title='Scaled\nGSVA score\n',
                        labels_gp=gpar(fontsize=25),
                        title_gp=gpar(fontsize=25),
                        legend_height=unit(6,'cm')),
                        row_names_side = 'left'
                        )
draw(p,heatmap_legend_side = "right",padding = unit(c(2, 20, 2, 2), "cm"))#margin
dev.off()
                   
#define a function to plot pathway-gene heatmap############
getheatdata <- function(pathway,region,mypathway){
  gsvaout <- mypathway[mypathway$Description==pathway,]
  overgenefinal <- str_split(mypathway[mypathway$Description==pathway,'geneID'],pattern='/')[[1]]
  submeta <- brain_slice1meta[brain_slice1meta$region2==region,]
  submeta <- submeta[submeta$Type=='4WT'| 
                        submeta$Type=='24WT'| 
                        submeta$Type=='24AD'|
                        submeta$Type=='24SL',]
  heatbox <- brain_slice1exp[overgenefinal,submeta$Spot]
  heatbox <- heatbox[rowSums(heatbox)>0,]
  annotation_col = data.frame(Type = submeta[,'Type'])   ##分组，legend名称：Type
  rownames(annotation_col) = submeta$Spot
  heats <- list(heatbox,annotation_col)
  return(heats)
}

#HPd通路中的基因热图，基因只需要按照SLvsAD排序
hpdpath <- read.table('Dendritic_enrichKEGG_keygene.txt',sep='\t',header=T,check.names=F)
hpdpath$Description <- sapply(hpdpath$Description,function(x){gsub(' \\- Mus musculus \\(house mouse\\)','',x)})                                                            
pathway <- 'Insulin secretion'
pathway <- 'Thyroid hormone signaling pathway'
pathway <- 'HIF-1 signaling pathway'
pathway <- 'Oxidative phosphorylation'
pathway <- 'Long-term depression'
pathway <- 'Long-term potentiation'
pathway <- 'Synaptic vesicle cycle'
pathway <- 'Alzheimer disease'
pathway <- 'Parkinson disease'
pathway <- 'Pathways of neurodegeneration - multiple diseases'
pathway <- 'Huntington disease'
pathway <- 'Prion disease'
pathway <- 'Amyotrophic lateral sclerosis'
pathway <- 'Circadian entrainment'
pathway <- 'Thermogenesis'
pathway <- 'Retrograde endocannabinoid signaling'
pathway <- 'Diabetic cardiomyopathy'
pathway <- 'Cardiac muscle contraction'
pathway <- 'Oocyte meiosis'
pathway <- 'Chemical carcinogenesis - reactive oxygen species'
pathway <- 'Vasopressin-regulated water reabsorption'
pathway <- 'Ferroptosis'
pathway <- 'Collecting duct acid secretion'
pathway <- 'Non-alcoholic fatty liver disease'
pathway <- 'Cardiac muscle contraction'
region <- 'Dendritic'
heatbox <- getheatdata(pathway,region,hpdpath)[[1]]#hpdpath就是HPd AAR 富集分析得到的通路

#order gene
boxorder <- resultHPF[resultHPF$gene %in% rownames(heatbox) & resultHPF$region=='Dendritic' & resultHPF$versus=='SLvsAD',] %>% .[order(.$logFC,decreasing=T),]
heatbox <- heatbox[boxorder$gene,]
#heatmap of genes in pathway
#Thyroid hormone signaling pathway
ann <- data.frame(brain_slice1meta[colnames(heatbox),'Type'])
colnames(ann) <- 'Type'
colours <- list('Type' = c('4WT'="gray",'24WT'='black',
                          '24AD'='#de1e17','24SL'='blue'))
colAnn  <- ComplexHeatmap::HeatmapAnnotation(clsuter=anno_block(
  gp=gpar(fill=c("gray",'black','#de1e17','blue')),
  labels=c('4WT','24WT','24AD','24SL'),
  labels_gp = gpar(col = "white", fontsize = 30,fontface='bold') 
))
colsplit <- c(rep(1,as.numeric(table(ann$Type)['4WT'])),
          rep(2,as.numeric(table(ann$Type)['24WT'])),
          rep(3,as.numeric(table(ann$Type)['24AD'])),
          rep(4,as.numeric(table(ann$Type)['24SL'])))
scaleheatbox <-  heatbox %>% pheatmap:::scale_rows(.)
colors2 <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
colors <- colorRampPalette(colors =colors2 )(100)
n_col = length(colors)
lim = max(abs(scaleheatbox), na.rm = TRUE)
library(colorRamp2)
colors2 <- colorRamp2(seq(-lim, lim, length = n_col), colors)
pdf(paste0(region,'_',pathway,'.pdf'),width=15,height=45)
ComplexHeatmap::Heatmap(as.matrix(scaleheatbox),show_column_names=FALSE,cluster_columns = FALSE,cluster_rows = FALSE,
                      show_row_dend = FALSE,show_column_dend = FALSE,
                      show_row_names = TRUE,
                      border=TRUE,border_gp = gpar(col = "black",lwd=1.5,lty=1),#border
                      col = colors2, top_annotation=colAnn,
                      width = unit(18, "cm"), height = unit(50, "cm"),
                      column_split=colsplit,
                      # row_split=rowsplit,
                      column_title = NULL,
                      # row_title = NULL,
                      # row_names_gp = gpar(fontsize = 25,col=c('red',rep('black',7))),
                      row_names_gp = gpar(fontsize = 25),
                      heatmap_legend_param=list(title='Scaled Expression\n',
                      labels_gp=gpar(fontsize=25),
                      title_gp=gpar(fontsize=25),
                      legend_height=unit(6,'cm'))
                      )
dev.off()
