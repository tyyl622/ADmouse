library(patchwork)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(tidyr)
library(stringr)
library(stringi)

#step1.Getsvg###########################
setwd('svgbit')
st.features <- c()
for (sam in c('AD335','SL298','WT731','WT281')){
    genes <- read.csv(paste0(sam,"_svg.csv"),header=T,check.names=F)
    colnames(genes)[1] <- 'gene'
    delete <- c(grep('^mt-',genes))#,grep('^Hb',genes),grep('^Rp[sl]',genes)
    if (length(delete)>0){
    genes <- genes[-delete] 
    } else{
    genes <- genes
    }
    genes <- genes[order(genes$AI,decreasing=T),] %>% .$gene %>% .[1:500]
    st.features <- c(st.features,genes) %>% unique(.)
}
save(st.features,file='step1.st.features.rda')

#integrate
setwd('preSample')
print('load AD335')
load('AD335_ob.rda')
AD335_ob <- object
AD335_ob$Type <- '24AD'

print('load SL298')
load('SL298_ob.rda')
SL298_ob <- object
SL298_ob$Type <- '24SL'

print('load WT281')
load('WT281_ob.rda')
WT281_ob <- object
WT281_ob$Type <- '24WT'

print('load WT731')
load('WT731_ob.rda')
WT731_ob <- object
WT731_ob$Type <- '4WT'

#filt
filtob <- function(object,nfeature,permt){
    obnew <- PercentageFeatureSet(object ,"^mt-", col.name = "percent.mt")
    obnew <- PercentageFeatureSet(obnew ,"^Hb", col.name = "percent.Hb")
    obnew <- PercentageFeatureSet(obnew ,"^Rp[sl]", col.name = "percent.Rpsl")
    obnew <- obnew[,obnew$nFeature_Spatial > nfeature & 
                            obnew$percent.mt < permt]
    counts <- GetAssayData(obnew, assay = "Spatial")
    de <- grep('^mt-',rownames(counts))
    counts <- counts[-de,]
    obnew <- subset(obnew, features = rownames(counts))
    obnew <- SCTransform(obnew, assay = "Spatial", 
                        verbose = FALSE,vars.to.regress='percent.mt')
    sampName <- obnew$Sample 
    newname <- paste0(sampName,'_',colnames(obnew))
    obnew <- RenameCells(obnew,new.names=newname)   
    return(obnew)
}
print('filt')
WT731_ob <- filtob(WT731_ob,200,30)
WT281_ob<- filtob(WT281_ob,200,30)
AD335_ob <- filtob(AD335_ob,200,30)
SL298_ob <- filtob(SL298_ob,200,30)

#Analysis
print('creat list')
st.list = list(WT731_ob,WT281_ob,AD335_ob,SL298_ob)
options(future.globals.maxSize = 8000 * 1024^2)
save(st.list,file='step2.st.list.rda')

#SCTmodel gene
load('step1.st.features.rda')
sctgene <- lapply(1:4,function(x){
  rownames(st.list[[x]]@assays$SCT@SCTModel.list$model1@feature.attributes)
})
sctuse <- sctgene[[1]]
for (i in 2:4){
  sctuse <- intersect(sctuse,sctgene[[i]])
}
st.features <- intersect(sctuse,st.features)
save(st.features,file='st.features.final.rda')

st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features, 
                              verbose = FALSE)
#integration
print('FindAnchors')
int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT", 
                                      verbose = FALSE, anchor.features = st.features)

print('Integrate')
brain_slice1 <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", 
                                  verbose = FALSE)
save(brain_slice1,file='step2.brain_slice1.rda')

#get clusters
pc <- 18
reso <- 0.5
brain_slice1_new1 <- FindNeighbors(brain_slice1, dims = 1:pc)
brain_slice1_new1 <- FindClusters(brain_slice1_new1, verbose = FALSE,resolution = reso)
brain_slice1_new1 <- RunUMAP(brain_slice1_new1, reduction = "pca",dims = 1:pc)
brain_slice1_new1 <- RunTSNE(brain_slice1_new1, reduction = "pca", dims = 1:pc)
colourCount <- length(unique(brain_slice1_new1$seurat_clusters))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
mycol <- getPalette(colourCount)
names(mycol) <- unique(brain_slice1_new1$seurat_clusters)[order(unique(brain_slice1_new1$seurat_clusters))]
pdf(paste0('Dim_pc',pc,'_reso',reso,'.pdf'),
    width = 15,height = 7)
p1 <- DimPlot(brain_slice1_new1, reduction = "tsne", group.by =  "ident",cols=mycol)
p2 <- DimPlot(brain_slice1_new1, reduction = "umap", group.by =  "ident",cols=mycol)
print(p1+p2)
dev.off()

pdf(paste0('DimGroup_pc',pc,'_reso',reso,'.pdf'),
    width = 15,height = 7)
p1 <- DimPlot(brain_slice1_new1, reduction = "tsne", group.by =  "Type")
p2 <- DimPlot(brain_slice1_new1, reduction = "umap", group.by =  "Type")
print(p1+p2)
dev.off()

metadata <- brain_slice1_new1@meta.data
metadata$region <- as.vector(metadata$seurat_clusters)
# metadata[rownames(remeta),'region'] <- remeta$reregion
metadata$region <- gsub('^0$','C0',metadata$region)
metadata$region <- gsub('^1$','C1',metadata$region)
metadata$region <- gsub('^2$','C1',metadata$region)
metadata$region <- gsub('^3$','Layer1/2',metadata$region)
metadata$region <- gsub('^4$','Dendritic',metadata$region)
metadata$region <- gsub('^5$','C2',metadata$region)
metadata$region <- gsub('^6$','C3',metadata$region)
metadata$region <- gsub('^7$','Layer6',metadata$region)
metadata$region <- gsub('^8$','Layer3/4',metadata$region)
metadata$region <- gsub('^9$','Layer5',metadata$region)
metadata$region <- gsub('^10$','Layer1/2',metadata$region)
metadata$region <- gsub('^11$','C4',metadata$region)
metadata$region <- gsub('^12$','C5',metadata$region)
metadata$region <- gsub('^13$','C2',metadata$region)
metadata$region <- gsub('^14$','GrDG',metadata$region)
metadata$region <- gsub('^15$','CA2/3',metadata$region)
metadata$region <- gsub('^16$','CA1',metadata$region)
metadata$region <- gsub('^17$','C0',metadata$region)

metadata$region2 <- metadata$region
metadata$region2 <- gsub('CA2/3','Somatic',metadata$region2)
metadata$region2 <- gsub('CA1','Somatic',metadata$region2)
metadata$region2 <- gsub('GrDG','Somatic',metadata$region2)
metadata$Spot <- rownames(metadata)
brain_slice1_final <- brain_slice1_new1
brain_slice1_final@meta.data <- metadata
Idents(brain_slice1_final) <- brain_slice1_final$region

pc <- 18
reso <- 0.5
load('/public/home/tangy/AD_ST/result/Seurat/filt/Finalmerge/mycol.rda')
names(mycol)[3] <- 'GrDG'
names(mycol)[6] <- 'Layer1/2'
names(mycol)[7] <- 'Layer3/4'
names(mycol)[8] <- 'Layer5'
names(mycol)[9] <- 'Layer6'
mycol['C4'] <- 'Gold'

Idents(brain_slice1_final) <- brain_slice1_final$region
# pdf(paste0('Dim_pc',pc,'_reso',reso,'_final.pdf'),width = 15,height = 7)
p1 <- DimPlot(brain_slice1_final, reduction = "tsne", group.by =  "ident",cols=mycol,pt.size=1.5)+
                theme(legend.text=element_text(size=27),
                      axis.title=element_text(size=23),
                      axis.text=element_text(size=23))+ggtitle('')
p2 <- DimPlot(brain_slice1_final, reduction = "umap", group.by =  "ident",cols=mycol,pt.size=1.5)+
                theme(legend.text=element_text(size=27),
                      axis.title=element_text(size=23),
                      axis.text=element_text(size=23))+ggtitle('')
ggsave(p1,filename=paste0('tsne_pc',pc,'_reso',reso,'_final.tiff'),width=8.5,height=7,dpi=600)
ggsave(p2,filename=paste0('umap_pc',pc,'_reso',reso,'_final.tiff'),width=8.5,height=7,dpi=600)

p2 <- SpatialDimPlot(brain_slice1_final,stroke=NA,label = FALSE,
                    cols=mycol,pt.size.factor = 1.2,crop=FALSE)+ 
  patchwork::plot_layout(ncol = 4, nrow = 1,guides='collect')
ggsave(p2,filename = paste0('cluster.pc',pc,'_reso',reso,'_final_unlabel.png'),
        width = 15,height = 4)
ggsave(p2,filename = paste0('cluster.pc',pc,'_reso',reso,'_final_unlabel.pdf'),
        width = 15,height = 4)
save(brain_slice1_final,file='brain_slice1_final.rda')

#findallmarker
setwd('allmarker')
Idents(brain_slice1_final) <- brain_slice1_final$region
allmarkers <- FindAllMarkers(brain_slice1_final, only.pos = FALSE, min.pct = 0.2, 
                        logfc.threshold = 0.25,assay='integrated')
save(allmarkers,file='allmarkersRegion.rda')
allmarkerfilt <- allmarkers[allmarkers$p_val_adj<0.05,]
save(allmarkerfilt,file='allmarkerfiltRegion.rda')
write.csv(allmarkerfilt,file='allmarkerfiltRegion.csv',row.names=F,quote=T)

#AAR
brain_slice1exp <- brain_slice1_final@assays$SCT@data
brain_slice1meta <- brain_slice1_final@meta.data
save(brain_slice1exp,file='brain_slice1exp.rda')
save(brain_slice1meta,file='brain_slice1meta.rda')

#####4.1 差异分析#########################
print('4.1')
group <- list(c('24SL','24AD'),
                c('24SL','24WT'),
                c('24AD','24WT'),
                c('24WT','4WT'))
myDEG <- function(expdata,mygenes,mymeta,group){
    DEG_sta <- data.frame()
    for(i in unique(mymeta$region2)){
        for(versus in group ){
            print(paste0(i,'_',versus))
            for(gene in mygenes){
                    DEG_tmp1 <- expdata[gene,mymeta[mymeta$region2==i & mymeta$Type==versus[1],'Spot']]  %>%
                                data.frame(.,check.names = F)
                    colnames(DEG_tmp1)[1] <- 'exp'
                    DEG_tmp1$region <- i
                    DEG_tmp1$Type <- versus[1]
                    DEG_tmp2 <- expdata[gene,mymeta[mymeta$region2==i & mymeta$Type==versus[2],'Spot']]  %>%
                                data.frame(.,check.names = F)
                    colnames(DEG_tmp2)[1] <- 'exp'
                    DEG_tmp2$region <- i
                    DEG_tmp2$Type <- versus[2]
                    DEG_tmp3 <- rbind(DEG_tmp1,DEG_tmp2)
                    sta_value <- wilcox.test(exp ~ Type, DEG_tmp3)
                    sta_tmp <- data.frame(gene=gene,region = i,versus = paste0(versus[1],'vs',versus[2]),
                                        pvalue =sta_value$p.value)
                    mean1 <- mean(DEG_tmp1$exp)
                    mean2 <- mean(DEG_tmp2$exp)
                    sta_tmp$logFC <- log2(mean1/mean2)
                    sta_tmp$trend <- ifelse(mean1 > mean2,"up",'down')
                    DEG_sta <- rbind(DEG_sta,sta_tmp)
            }
        }
    }
    return(DEG_sta)
}
#Hippocampus####################
HPF <- c('Somatic','Dendritic')
HPFMeta <- brain_slice1meta[brain_slice1meta$region2 %in% HPF,]
#filt genes
filtexp <- function(exp,genes,cut){
    tmpexp <- rowSums(exp[genes,])
    tmpexp <- tmpexp[tmpexp<=cut]
    gfilt <- unique(names(tmpexp))
    return(gfilt)
} 
deletebyexp_HPF <- filtexp(brain_slice1exp,rownames(brain_slice1exp),25)
deletes <- c('Gapdh$','Actb','^mt-','^Rp[sl]','^Hb','Xist','Fis1')
allgene <- rownames(brain_slice1exp)
mydelete1 <- sapply(deletes,function(x){
    tmp <- allgene[grep(x,allgene)]
    return(tmp)
}) %>% unlist(.) %>% unique(.)
mydelete3 <- unique(c(mydelete1,deletebyexp_HPF))
if (length(which(allgene %in% mydelete3))>0){
        staingene <- allgene[-which(allgene %in% mydelete3)]}
DEGstep1_HPF <- myDEG(brain_slice1exp,staingene,HPFMeta,group)
save(DEGstep1_HPF,file='wilcox/DEGstep1_HPF.rda')

#Add pdaj
DEGprocess <- function(DEGs1,padjcut,filenames){
    #加padj
    DEGpadj <- data.frame()
    for (region in unique(DEGs1$region)){
        for (versus in group){
            tmp <- DEGs1[DEGs1$region==region & DEGs1$versus==paste0(versus[1],'vs',versus[2]),]
            tmp <- tmp[!is.infinite(tmp$logFC),]
            tmp$padj <- p.adjust(tmp$pvalue,method = 'fdr')
            DEGpadj <- rbind(DEGpadj,tmp)
        }
    }
    save(DEGpadj,file=filenames)
    
    #cut padj
    DEGpadjfilt <- DEGpadj[DEGpadj$padj<padjcut,]
    DEGpadjfilt <- na.omit(DEGpadjfilt)
    return(DEGpadjfilt)
}
padj <- 0.001
DEGstep2_HPF <- DEGprocess(DEGstep1_HPF,padj,'wilcox/DEG_padjAll_HPF.rda')
save(DEGstep2_HPF,file=paste0('wilcox/DEG_padj',padj,'filtHPF.rda'))

#other regions###############
other <- unique(brain_slice1meta$region)[-c(10,11)]
otherMeta <- brain_slice1meta[brain_slice1meta$region2 %in% other,]
deletebyexp_other <- filtexp(brain_slice1exp,rownames(brain_slice1exp),25)
deletes <- c('Gapdh$','Actb','^mt-','^Rp[sl]','^Hb','Xist','Fis1')
allgene <- rownames(brain_slice1exp)
mydelete1 <- sapply(deletes,function(x){
    tmp <- allgene[grep(x,allgene)]
    return(tmp)
}) %>% unlist(.) %>% unique(.)
mydelete3 <- unique(c(mydelete1,deletebyexp_other))
if (length(which(allgene %in% mydelete3))>0){
        staingene <- allgene[-which(allgene %in% mydelete3)]}
DEGstep1_other <- myDEG(brain_slice1exp,staingene,otherMeta,group)
save(DEGstep1_other,file='other/DEGstep1_other.rda')
padj <- 0.001
DEGstep2_other <- DEGprocess(DEGstep1_other,padj,'other/DEG_padjAll_other.rda')
save(DEGstep2_other,file=paste0('other/DEG_padj',padj,'filtother.rda'))
write.csv(DEGstep2_other,'other/差异基因信息other.csv',quote = T,row.names = F)#不用阈值筛选，所有基因差异信息

####find AAR##############################
findres <- function(DEGs2){
    result <- data.frame()
    for (region in unique(DEGs2$region)){
        tmp1 <- DEGs2[DEGs2$region==region,]
            for(gene in unique(tmp1$gene)){
                tmp <- tmp1[tmp1$gene == gene,]
                tmp$trend <- paste0(tmp$versus,'_',tmp$trend)
                if ('WT_24MvsWT_4M_up'  %in% tmp$trend & 'SLvsAD_up' %in% tmp$trend & 'ADvsWT_24M_down' %in% tmp$trend ){
                retain <- rbind(tmp[tmp$trend=='WT_24MvsWT_4M_up',],
                                tmp[tmp$trend=='SLvsAD_up',],
                                tmp[tmp$trend=='ADvsWT_24M_down',])
                
                } else if ('WT_24MvsWT_4M_down'  %in% tmp$trend & 'SLvsAD_down' %in% tmp$trend & 'ADvsWT_24M_up' %in% tmp$trend ){
                retain <- rbind(tmp[tmp$trend=='WT_24MvsWT_4M_down',],
                                tmp[tmp$trend=='SLvsAD_down',],
                                tmp[tmp$trend=='ADvsWT_24M_up',])
                } else {
                retain <- data.frame()
                }
                # print(gene)
                # print(region)
                result <- rbind(result,retain)
        }
    }
    return(result)
}
resultHPF <- findres(DEGstep2_HPF)
save(resultHPF,file='resultHPF.rda')

#enrichment##################################################
library(stringr)
library(clusterProfiler)
library("org.Mm.eg.db")
library("enrichplot")
library("ggplot2")
library(dplyr)
library(msigdbr)
for (i in unique(resultHPF$region)){
  if (length(gene)>0){
    geneID <- mapIds(org.Mm.eg.db,keys = gene,column = 'ENTREZID',
                     keytype = 'SYMBOL',multiVals='filter')
    kkbp <- enrichGO(gene = gene,
                     OrgDb = org.Mm.eg.db, 
                     pvalueCutoff =0.05, 
                     qvalueCutoff = 0.2,
                     ont="BP",
                     readable =F,
                     keyType = 'SYMBOL'
    ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
    save(kkbp,file=paste0(i,"_","enrichGOBP_keygene.rda"))
    
    kkcc <- enrichGO(gene = gene,
                     OrgDb = org.Mm.eg.db, 
                     pvalueCutoff =0.05, 
                     qvalueCutoff = 0.2,
                     ont="CC",
                     readable =F,
                     keyType = 'SYMBOL'
    ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
    save(kkcc,file=paste0(i,"_","enrichGOCC_keygene.rda"))
    
    kkmf <- enrichGO(gene = gene,
                     OrgDb = org.Mm.eg.db, 
                     pvalueCutoff =0.05, 
                     qvalueCutoff = 0.2,
                     ont="MF",
                     readable =F,
                     keyType = 'SYMBOL'
    ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min) 
    save(kkmf,file=paste0(i,"_","enrichGOMF_keygene.rda"))
    
    kkbp <- data.frame(kkbp,check.names = F)
    kkbp$Category <- 'BP'
    kkcc <-  data.frame(kkcc,check.names = F)
    kkcc$Category <- 'CC'
    kkmf <-  data.frame(kkmf,check.names = F)
    kkmf$Category <- 'MF'
    kk4 <- rbind(kkbp,kkcc,kkmf)
    write.table(kk4,file=paste0(i,"_","enrichGO_keygene.txt"),
                quote=F,row.names = F,sep = '\t') 
    k <- enrichKEGG(gene = geneID, organism = "mmu",
                    pvalueCutoff =0.05, qvalueCutoff =0.2,
                    keyType = 'kegg',minGSSize=5) %>%
      setReadable(., OrgDb = org.Mm.eg.db, keyType="ENTREZID") %>%
      data.frame(.,check.names = F)
    write.table(k,file=paste0(i,"_","enrichKEGG_keygene.txt"),
                quote=F,row.names = F,sep = '\t')
  }
}

#######################GSVA#######################################
library(plyr)
library(GSVA)
library(GSEABase)
library(ggplot2)
library(ExperimentHub)
library(org.Mm.eg.db)
library(limma)
load('brain_slice1exp.rda')#SCT matrix of 4 brains
load('brain_slice1meta.rda')#metadata of 4 brains
setwd('GSVA/')
gmt1 <- read.delim2('mousekeggGMT20231020.gmt',sep='\t',header=F,check.names=F)#download from KEGG

#GSVA分析
gmt <- getGmt('gmtForGSVA.gmt')
gsvaOut=gsva(data.matrix(brain_slice1exp), gmt, min.sz=10,###通路中最小基因数目 
             #max.sz=Inf,
             verbose=TRUE,parallel.sz=2)
gsvaOut=rbind(id=colnames(gsvaOut),gsvaOut)
write.table(gsvaOut,file="GSVAOut.txt",sep="\t",quote=F,col.names=F)

#limma统计差异通路
gsvameta <- brain_slice1meta
group <- list(c('24SL','24AD'),
                c('24SL','24WT'),
                c('24AD','24WT'),
                c('24WT','24WT'))
diffgsva <- function(gsvares,gsvameta,myregions){
    gsva_sta <- data.frame()
    for  (i in myregions){
      for(versus in group ){
        gsva_tmp1 <- gsvares[,c('id',gsvameta[gsvameta$region2==i & gsvameta$Type==versus[1],'Spot'])]
        gsva_tmp2 <- gsvares[,c('id',gsvameta[gsvameta$region2==i & gsvameta$Type==versus[2],'Spot'])]
        identical(gsva_tmp1$id,gsva_tmp2$id)
        gsva <- cbind(gsva_tmp1,gsva_tmp2[,-1])
        rownames(gsva) <- gsva$id
        gsva <- gsva[,-1]      
        group_list <- data.frame(cell = colnames(gsva), group = c(rep('V1', (ncol(gsva_tmp1)-1)), rep('V2', (ncol(gsva_tmp2)-1))))
        group_list$group <- factor(group_list$group,levels=c('V1','V2'))
        design <- model.matrix(~ 0 + group_list$group)
        colnames(design) <- levels(group_list$group)
        rownames(design) <- colnames(gsva)
        # design
        # 构建差异比较矩阵
        contrast.matrix <- makeContrasts(V1-V2, levels = design)
        # 差异分析，b vs. a
        fit <- lmFit(gsva, design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
        x$versus <- paste0(versus[1],'vs',versus[2])
        x$trend <- ifelse(x$logFC>0,paste0(x$versus,'_up'),paste0(x$versus,'_down'))
        x$region <- i
        x$pathway <- rownames(x)
        gsva_sta <- rbind(gsva_sta,x)
    }

  }
  return(gsva_sta)
}
gsvaout <- read.table('GSVAOut.txt',sep = '\t',header =T,check.names=F)
gsvarlimma_hpf <- diffgsva(gsvaout,gsvameta,'Dendritic')
save(gsvarlimma_hpf,file='gsvarlimma_HPdentritic.rda')

#filt AAR pathway
gsvafilt <- function(gsvalimma,myregions,pcut,tcut){
  gsva_filt <- data.frame()
  for (region in myregions){
    tmp1 <- gsvalimma[gsvalimma$region==region,]
    for(path in unique(tmp1$pathway)){
      tmp <- tmp1[tmp1$pathway==path,]
      if ('24WTvs4WT_up'  %in% tmp$trend & 
          '24SLvs24AD_up' %in% tmp$trend & 
          '24ADvs24WT_down' %in% tmp$trend &
          tmp$adj.P.Val[grep('24WTvs4WT',tmp$versus)]<pcut &
          tmp$adj.P.Val[grep('24SLvs24AD',tmp$versus)]<pcut &
          tmp$adj.P.Val[grep('24ADvs24WT',tmp$versus)]<pcut &
          abs(tmp$t[grep('24WTvs4WT',tmp$versus)])>=tcut &
          abs(tmp$t[grep('24SLvs24AD',tmp$versus)])>=tcut &
          abs(tmp$t[grep('24ADvs24WT',tmp$trend)])>=tcut ){
        retain <- rbind(tmp[tmp$trend=='24WTvs4WT_up',],
                        tmp[tmp$trend=='24SLvs24AD_up',],
                        tmp[tmp$trend=='24ADvs24WT_down',])
        
      } else if ('24WTvs4WT_down'  %in% tmp$trend & 
                '24SLvs24AD_down' %in% tmp$trend & 
                '24ADvs24WT_up' %in% tmp$trend &
                  tmp$adj.P.Val[grep('24WTvs4WT',tmp$versus)]<pcut &
                  tmp$adj.P.Val[grep('24SLvs24AD',tmp$versus)]<pcut &
                  tmp$adj.P.Val[grep('24ADvs24WT',tmp$versus)]<pcut &
                  abs(tmp$t[grep('24WTvs4WT',tmp$versus)])>=tcut &
                  abs(tmp$t[grep('24SLvs24AD',tmp$versus)])>=tcut &
                  abs(tmp$t[grep('24ADvs24WT',tmp$versus)])>=tcut  ){
        retain <- rbind(tmp[tmp$trend=='24WTvs4WT_down',],
                        tmp[tmp$trend=='24SLvs24AD_down',],
                        tmp[tmp$trend=='24ADvs24WT_up',])
      } else {
        retain <- data.frame()
      }
      # print(gene)
      # print(region)
      gsva_filt <- rbind(gsva_filt,retain)
    }
  }
  return(gsva_filt)
}
gsvafilt_hpf <- gsvafilt(gsvarlimma_hpf,'Dendritic',0.05,2)
save(gsvafilt_hpf,file='gsvafilt_HPd.rda')
write.csv(gsvafilt_hpf,file = 'gsvaFilt_HPdentritic.csv',row.names = F,quote = T)#

#Somatic
gsvarlimma_hps<- diffgsva(gsvaout,gsvameta,'Somatic')
save(gsvarlimma_hps,file='gsvarlimma_Somatic.rda')
gsvafilt_hps <- gsvafilt(gsvarlimma_hps,'Somatic',0.05,2)
save(gsvafilt_hps,file='gsvafilt_HPs.rda')
write.csv(gsvafilt_hps,file = 'gsvaFilt_Somatic.csv',row.names = F,quote = T)#T
