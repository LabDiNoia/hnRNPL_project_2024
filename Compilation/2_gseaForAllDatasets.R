# get Hallmark genesets and format them for fgsea
library(msigdbr)

# function to remove duplicates and NAs from a vector
rm_dups_NAs<-function(vec){
	vec<-vec[!duplicated(vec)]
	vec<-vec[!is.na(vec)]
	return(vec)
}

Hs_hallmark_gs_df<-msigdbr(species = "human", category = "H")
Hs_hallmark_gs<-split(x=Hs_hallmark_gs_df$ensembl_gene,f=Hs_hallmark_gs_df$gs_name)
Hs_hallmark_gs<-lapply(Hs_hallmark_gs,rm_dups_NAs)

Mm_hallmark_gs_df<-msigdbr(species = "mouse", category = "H")
Mm_hallmark_gs<-split(x=Mm_hallmark_gs_df$ensembl_gene,f=Mm_hallmark_gs_df$gs_name)
Mm_hallmark_gs<-lapply(Mm_hallmark_gs,rm_dups_NAs)


# get datasets and order them by foldChange or Wald statistic (stat column in DESeq2 output)
library(data.table)
allNormcounts_dt<-readRDS("./R_rds_files/allNormcounts.rds")
allNormcounts_dt[,gene_id_noVer:=sub('\\..*','',GeneID),by=.I]

datasets<-allNormcounts_dt[,unique(dataset)]

do_fgsea<-function(whichDataset){
	cur_dataset=whichDataset
	# remove commonID column to remove duplicates
	cur_dt<-allNormcounts_dt[dataset==cur_dataset]
	cur_dt<-unique(cur_dt[,.SD,.SDcols=!c("commonID")])
	#remove pseudoautosomal region IDs, cos they cause duplicates downstream
	cur_dt<-cur_dt[!grep("PAR_Y",GeneID)]
	# remove padj=NA to remove low expression values and order according to stat
	stat_dt<-unique(cur_dt[!is.na(padj),.(gene_id_noVer,log2FoldChange)][order(log2FoldChange)])
	# make named vector
	ordered_stat<-stat_dt$log2FoldChange
	names(ordered_stat)<-stat_dt$gene_id_noVer
	
	cat("Performing GSEA on",length(ordered_stat),"genes in the dataset:", cur_dataset,"\n")

	require(fgsea)
	scoretype<-"std"
	cur_genome<-cur_dt[,unique(mouse_or_human)]
	cur_pathways<-if(cur_genome=="human"){Hs_hallmark_gs
			}else if(cur_genome=="mouse"){Mm_hallmark_gs}
	res<-fgsea(pathways=cur_pathways, stats=ordered_stat, scoreType=scoretype, nPermSimple=10000)
	res[,dataset:=cur_dataset]
	return(res)
}

allGSEA_dt<-sapply(datasets,do_fgsea,simplify=F)
allGSEA_dt<-rbindlist(allGSEA_dt)
saveRDS(allGSEA_dt,"./R_rds_files/allGSEA_dt.rds")
require(data.table)
allGSEA_dt<-readRDS("./R_rds_files/allGSEA_dt.rds")

p53_leadingEdge<-allGSEA_dt[dataset=="actb_d1" & pathway=="HALLMARK_P53_PATHWAY", unlist(leadingEdge)]
allNormcounts_dt[dataset=="actb_d1"][gene_id_noVer %in% p53_leadingEdge][order(log2FoldChange,decreasing = T)]
write.table(allNormcounts_dt[dataset=="actb_d1"][gene_id_noVer %in% p53_leadingEdge][order(log2FoldChange,decreasing = T)],
            "./R_rds_files/HM_p53_leadingEdge_actB.txt", sep="\t", quote=F, row.names=F)

# make comparative GSEA heatmaps
#1. create matrix to plot?
# dcast to convert to wide format. take only padj<0.1
gsea_mat<-as.matrix(dcast(allGSEA_dt[padj<0.1], pathway~dataset, value.var = "NES"),rownames = T)
rownames(gsea_mat)<-gsub("_"," ",sub("HALLMARK_","",rownames(gsea_mat)),)
gsea_mat<-cbind(gsea_mat, thymocytes=rep(NA,nrow(gsea_mat)))
gsea_mat[is.na(gsea_mat)]<-0 # replace NAs with 0
gsea_mat<-gsea_mat[rowSums(gsea_mat!=0)>1,] # keep columns with 2 sig values

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
colorScaleFunction<-colorRamp2(c(-2,0,2),c("blue","white","red"))
#complicaed to make square cells. will do in illustrator
Heatmap(gsea_mat,col = colorScaleFunction)
# need to order based on the story. also illustrator?
gsea_mat<-gsea_mat[,c("thymocytes","fetalLiv","actb_d1","kcytes","hepg2","k562","lncap","bjFibro","HEK293T")]
#gsea_mat<-gsea_mat[,c("HEK293T","bjFibro","lncap","k562","hepg2","kcytes","actb_d1","fetalLiv","thymocytes")]
gsea_mat<-gsea_mat[c("P53 PATHWAY","E2F TARGETS","MYC TARGETS V1","MYC TARGETS V2","MITOTIC SPINDLE","G2M CHECKPOINT",
                     "MTORC1 SIGNALING","KRAS SIGNALING DN","KRAS SIGNALING UP","IL2 STAT5 SIGNALING","IL6 JAK STAT3 SIGNALING",
                     "INFLAMMATORY RESPONSE","INTERFERON ALPHA RESPONSE","INTERFERON GAMMA RESPONSE","TGF BETA SIGNALING",
                     "CHOLESTEROL HOMEOSTASIS","PEROXISOME","BILE ACID METABOLISM","OXIDATIVE PHOSPHORYLATION",
                     "UNFOLDED PROTEIN RESPONSE","PROTEIN SECRETION","UV RESPONSE DN",
                     "APICAL JUNCTION","APICAL SURFACE","XENOBIOTIC METABOLISM","COMPLEMENT","COAGULATION",
                     "EPITHELIAL MESENCHYMAL TRANSITION","MYOGENESIS","ESTROGEN RESPONSE EARLY","ANDROGEN RESPONSE"),]
library(stringr)
rownames(gsea_mat)<-str_to_sentence(rownames(gsea_mat))
pdf(file="./figures/GSEA_heatmaps.pdf",useDingbats=F,width=8,height=8)
Heatmap(gsea_mat, col=colorScaleFunction, cluster_columns=F, cluster_rows=F, 
        name="GSEA (NES)", height=unit(6.2,"cm"), width=unit(1.8,"cm"),
        row_names_gp=gpar(fontsize=5.5), column_names_gp=gpar(fontsize=5.5))
dev.off()

# code to extract specific info about leadingEdge genes
test<-res[padj<0.1][order(NES)][1,leadingEdge][[1]]
allNormcounts_dt[gene_id_noVer%in%test & dataset=="kcytes"]
allNormcounts_dt[gene_id_noVer%in%test & dataset=="kcytes",sort(log2FoldChange)]

#---------------------------------------------------#
# make Hallmark GSEA NES vs padj plot for actB dataset using old_code

actb_GSEA_dt<-allGSEA_dt[dataset=="actb_d1"]
actb_GSEA_dt[padj<=0.1,cex_pt:=2]
actb_GSEA_dt[padj>0.1,cex_pt:=0.5]

library(RColorBrewer)
my_palette<-colorRampPalette(c("blue","white","red"))(299)
abs_max_val<-max(abs(min(actb_GSEA_dt[,NES])),max(actb_GSEA_dt[,NES]))
actb_GSEA_dt[,col_pt:=my_palette[findInterval(NES,seq(-abs_max_val,abs_max_val,length.out=299))]]
actb_GSEA_dt[padj>0.1,col_pt:="darkgrey"]
actb_GSEA_dt<-actb_GSEA_dt[order(NES)]

library(ggplot2)
library(grid)
p<-ggplot(transform(actb_GSEA_dt,fct=cut(NES,c(-3,-2.4,-2,1.5,2))),aes(padj,NES))+
  geom_vline(aes(xintercept=0.1),cex=0.5,linetype=2,color="grey")+
  geom_point(size=actb_GSEA_dt[,cex_pt],color=actb_GSEA_dt[,col_pt])+
  geom_text(aes(label=ifelse(padj<0.1,plotLabel,"")),cex=3,hjust=-0.04)+
  xlim(c(0,1))+
  facet_wrap(~fct,strip.position="right",ncol=1,scales="free_y",as.table=F)+
  theme_classic()+
  theme(text=element_text(family="sans"),
        plot.title=element_text(size=8),
        axis.line=element_line(color='black',size=0.25),
        axis.ticks=element_line(colour="black",size = 0.25),
        axis.title=element_text(size=7),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=5),
        strip.background=element_rect(size=0.25),
        strip.text=element_text(size=7),
        panel.border=element_rect(fill=NA,size=NA),
        panel.spacing=unit(0.25,"line"),
        legend.position="none")

gt<-ggplot_gtable(ggplot_build(p))
#gtable::gtable_show_layout(gt)
#gt$heights[7] = 0.15*gt$heights[7]
gt$heights[11] = 0.5*gt$heights[11]
gt$heights[15] = 0.5*gt$heights[15]
gt$heights[19] = 0.5*gt$heights[19]
#gt$heights[23] = 1*gt$heights[23]
#pdf("figures/actB_GSEA_nicePlot.pdf",useDingbats=F,width=2.5,height=2)
grid.draw(gt)
#dev.off()

#get gene names
cat(unique(allNormcounts_dt[gene_id_noVer %in% unlist(allGSEA_dt[pathway=="HALLMARK_XENOBIOTIC_METABOLISM" & dataset=="actb_d1", leadingEdge]), gene_name ]), sep="\n")
cat(unique(allNormcounts_dt[gene_id_noVer %in% unlist(allGSEA_dt[pathway=="HALLMARK_MYC_TARGETS_V1" & dataset=="actb_d1", leadingEdge]), gene_name ]), sep="\n")

#------------------------------------------------------------------------------------------------#
# GSEA using fgsea for other genesets
########################################################################################
# 1. TF motif GSEA
Hs_motif_gs_df<-rbind(msigdbr(species = "human", category = "C3", subcategory = "TFT:GTRD"),
                      msigdbr(species = "human", category = "C3", subcategory = "TFT:TFT_Legacy"))
Hs_motif_gs<-split(x=Hs_motif_gs_df$ensembl_gene,f=Hs_motif_gs_df$gs_name)
Hs_motif_gs<-lapply(Hs_motif_gs,rm_dups_NAs)

Mm_motif_gs_df<-rbind(msigdbr(species = "mouse", category = "C3", subcategory = "TFT:GTRD"),
                      msigdbr(species = "mouse", category = "C3", subcategory = "TFT:TFT_Legacy"))
Mm_motif_gs<-split(x=Mm_motif_gs_df$ensembl_gene,f=Mm_motif_gs_df$gs_name)
Mm_motif_gs<-lapply(Mm_motif_gs,rm_dups_NAs)

do_motif_fgsea<-function(whichDataset){
  cur_dataset=whichDataset
  # remove commonID column to remove duplicates
  cur_dt<-allNormcounts_dt[dataset==cur_dataset]
  cur_dt<-unique(cur_dt[,.SD,.SDcols=!c("commonID")])
  #remove pseudoautosomal region IDs, cos they cause duplicates downstream
  cur_dt<-cur_dt[!grep("PAR_Y",GeneID)]
  # remove padj=NA to remove low expression values and order according to stat
  stat_dt<-unique(cur_dt[!is.na(padj),.(gene_id_noVer,log2FoldChange)][order(log2FoldChange)])
  # make named vector
  ordered_stat<-stat_dt$log2FoldChange
  names(ordered_stat)<-stat_dt$gene_id_noVer
  
  cat("Performing GSEA on",length(ordered_stat),"genes in the dataset:", cur_dataset,"\n")
  
  require(fgsea)
  scoretype<-"std"
  cur_genome<-cur_dt[,unique(mouse_or_human)]
  cur_pathways<-if(cur_genome=="human"){Hs_motif_gs
  }else if(cur_genome=="mouse"){Mm_motif_gs}
  res<-fgsea(pathways=cur_pathways, stats=ordered_stat, scoreType=scoretype, nPermSimple=10000)
  res[,dataset:=cur_dataset]
  return(res)
}

motifGSEA_dt<-sapply(datasets,do_motif_fgsea,simplify=F)
motifGSEA_dt<-rbindlist(motifGSEA_dt)
saveRDS(motifGSEA_dt,"./R_rds_files/all_TFmotif_GSEA_dt.rds")
require(data.table)
motifGSEA_dt<-readRDS("./R_rds_files/all_TFmotif_GSEA_dt.rds")

unique(motifGSEA_dt[padj<0.1,.(dataset,.N),by=c("pathway")][,.(pathway,N)][order(N)])
unique(motifGSEA_dt[padj<0.1,.(dataset,.N),by=c("pathway")][,.(pathway,N)][order(N,decreasing = T)])[N>1][1:20]

########################################################################################
# GSEA using fgsea for other genesets
# 2. GO BP GSEA

Hs_gobp_gs_df<-msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
Hs_gobp_gs<-split(x=Hs_gobp_gs_df$ensembl_gene,f=Hs_gobp_gs_df$gs_name)
Hs_gobp_gs<-lapply(Hs_gobp_gs,rm_dups_NAs)

Mm_gobp_gs_df<-msigdbr(species = "mouse", category = "C5", subcategory = "GO:BP")
Mm_gobp_gs<-split(x=Mm_gobp_gs_df$ensembl_gene,f=Mm_gobp_gs_df$gs_name)
Mm_gobp_gs<-lapply(Mm_gobp_gs,rm_dups_NAs)

do_gobp_fgsea<-function(whichDataset){
  cur_dataset=whichDataset
  # remove commonID column to remove duplicates
  cur_dt<-allNormcounts_dt[dataset==cur_dataset]
  cur_dt<-unique(cur_dt[,.SD,.SDcols=!c("commonID")])
  #remove pseudoautosomal region IDs, cos they cause duplicates downstream
  cur_dt<-cur_dt[!grep("PAR_Y",GeneID)]
  # remove padj=NA to remove low expression values and order according to stat
  stat_dt<-unique(cur_dt[!is.na(padj),.(gene_id_noVer,log2FoldChange)][order(log2FoldChange)])
  # make named vector
  ordered_stat<-stat_dt$log2FoldChange
  names(ordered_stat)<-stat_dt$gene_id_noVer
  
  cat("Performing GSEA on",length(ordered_stat),"genes in the dataset:", cur_dataset,"\n")
  
  require(fgsea)
  scoretype<-"std"
  cur_genome<-cur_dt[,unique(mouse_or_human)]
  cur_pathways<-if(cur_genome=="human"){Hs_gobp_gs
  }else if(cur_genome=="mouse"){Mm_gobp_gs}
  res<-fgsea(pathways=cur_pathways, stats=ordered_stat, scoreType=scoretype, nPermSimple=10000)
  res[,dataset:=cur_dataset]
  return(res)
}


gobpGSEA_dt<-sapply(datasets,do_gobp_fgsea,simplify=F)
gobpGSEA_dt<-rbindlist(gobpGSEA_dt)
saveRDS(gobpGSEA_dt,"./R_rds_files/all_gobp_GSEA_dt.rds")
require(data.table)
gobpGSEA_dt<-readRDS("./R_rds_files/all_gobp_GSEA_dt.rds")

unique(gobpGSEA_dt[padj<0.1,.(dataset,.N),by=c("pathway")][N>=4,.(pathway,N)][order(N,decreasing = T)])
########################################################################################
# GSEA using fgsea for other Mitochondria related genesets
# 3. Mitochondria related genesets
Hs_all_df<-msigdbr(species = "human")
Hs_all_dt<-as.data.table(Hs_all_df)
rm(Hs_all_df)
# Hs_interestingTerms<-rbind(Hs_all_dt[grep("Mitochondri",gs_name,ignore.case = T)],
#                            Hs_all_dt[grep("Autophag",gs_name,ignore.case = T)],
#                            Hs_all_dt[grep("Mitophag",gs_name,ignore.case = T)],
#                            Hs_all_dt[grep("Redox",gs_name,ignore.case = T)],
#                            Hs_all_dt[grep("Oxid",gs_name,ignore.case = T)])

Hs_interestingTerms<-Hs_all_dt[grep("Oxid",gs_name,ignore.case = T)]

Hs_interesting_gs<-split(x=Hs_interestingTerms$ensembl_gene,f=Hs_interestingTerms$gs_name)
Hs_interesting_gs<-lapply(Hs_interesting_gs,rm_dups_NAs)
Hs_interesting_gs<-Hs_interesting_gs[sapply(Hs_interesting_gs,length) >= 20]
Hs_interesting_gs<-Hs_interesting_gs[sapply(Hs_interesting_gs,length) <= 500]
length(Hs_interesting_gs)

Mm_all_df<-msigdbr(species = "mouse")
Mm_all_dt<-as.data.table(Mm_all_df)
rm(Mm_all_df)
# Mm_interestingTerms<-rbind(Mm_all_dt[grep("Mitochondri",gs_name,ignore.case = T)],
#                            Mm_all_dt[grep("Autophag",gs_name,ignore.case = T)],
#                            Mm_all_dt[grep("Mitophag",gs_name,ignore.case = T)],
#                            Mm_all_dt[grep("Redox",gs_name,ignore.case = T)],
#                            Mm_all_dt[grep("Oxid",gs_name,ignore.case = T)])

Mm_interestingTerms<-Mm_all_dt[grep("Oxid",gs_name,ignore.case = T)]

Mm_interesting_gs<-split(x=Mm_interestingTerms$ensembl_gene,f=Mm_interestingTerms$gs_name)
Mm_interesting_gs<-lapply(Mm_interesting_gs,rm_dups_NAs)
Mm_interesting_gs<-Mm_interesting_gs[sapply(Mm_interesting_gs,length) >= 20]
Mm_interesting_gs<-Mm_interesting_gs[sapply(Mm_interesting_gs,length) <= 500]
length(Mm_interesting_gs)
intersect(names(Hs_interesting_gs),names(Mm_interesting_gs))

do_fgsea<-function(whichDataset,humanGeneSets,mouseGeneSets){
  cur_dataset=whichDataset
  # remove commonID column to remove duplicates
  cur_dt<-allNormcounts_dt[dataset==cur_dataset]
  cur_dt<-unique(cur_dt[,.SD,.SDcols=!c("commonID")])
  #remove pseudoautosomal region IDs, cos they cause duplicates downstream
  cur_dt<-cur_dt[!grep("PAR_Y",GeneID)]
  # remove padj=NA to remove low expression values and order according to stat
  stat_dt<-unique(cur_dt[!is.na(padj),.(gene_id_noVer,log2FoldChange)][order(log2FoldChange)])
  # make named vector
  ordered_stat<-stat_dt$log2FoldChange
  names(ordered_stat)<-stat_dt$gene_id_noVer
  
  cat("Performing GSEA on",length(ordered_stat),"genes in the dataset:", cur_dataset,"\n")
  
  require(fgsea)
  scoretype<-"std"
  cur_genome<-cur_dt[,unique(mouse_or_human)]
  cur_pathways<-if(cur_genome=="human"){humanGeneSets
  }else if(cur_genome=="mouse"){mouseGeneSets}
  res<-fgsea(pathways=cur_pathways, stats=ordered_stat, scoreType=scoretype, nPermSimple=1000)
  res[,dataset:=cur_dataset]
  return(res)
}


interesting_GSEA_dt<-sapply(datasets, do_fgsea, simplify=F,
                            humanGeneSets=Hs_interesting_gs, mouseGeneSets=Mm_interesting_gs)
interesting_GSEA_dt<-rbindlist(interesting_GSEA_dt)

saveRDS(interesting_GSEA_dt,"./R_rds_files/interesting_genesets_GSEA_dt.rds")
require(data.table)
interesting_GSEA_dt<-readRDS("./R_rds_files/interesting_genesets_GSEA_dt.rds")

unique(interesting_GSEA_dt[padj<0.1,.(dataset,.N),by=c("pathway")][N>=4,.(pathway,N)][order(N,decreasing = T)])
unique(interesting_GSEA_dt[pval<0.2,.(dataset,.N),by=c("pathway")][N>=4,.(pathway,N)][order(N,decreasing = T)])

interesting_GSEA_dt[dataset=="actb_d1" & pval<0.2][order(pval)]

#------------------------------------------------------------------------------------------------#

