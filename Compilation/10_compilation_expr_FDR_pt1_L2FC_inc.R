#library(rtracklayer)
#library(data.table)
#
## step0: make common IDs for mouse and human genes
#mouse_gtf<-import("../../common/mm10/gencode.vM23.primary_assembly.annotation.gtf")
#human_gtf<-import("../../common/hg38/gencode.v43.primary_assembly.annotation.gtf")
#
#mouse_dt<-data.table(mouse_id=unique(mouse_gtf$gene_id))
#mouse_dt[,mouse_id_noVer:=sub('\\.[0-9]*','',mouse_id,perl=T),by=.I]
#
#human_dt<-data.table(human_id=unique(human_gtf$gene_id))
#human_dt[,human_id_noVer:=sub('\\.[0-9]*','',human_id,perl=T),by=.I]
#
#mouse_human_homolog_dt<-fread("./refrFiles/human_mouse_orthologs.txt")
#colnames(mouse_human_homolog_dt)<-c("mouse_id_noVer","human_id_noVer")
#mouse_human_homolog_dt<-mouse_human_homolog_dt[human_id_noVer!=""]
## only 25719 genes have homologs in mouse and humans.
## verified by getting human->mouse homologs too. i get the same number of genes
#
#commonID_forHomologs_dt<-merge(mouse_dt,mouse_human_homolog_dt,all.x=T)
#commonID_forHomologs_dt<-merge(human_dt,commonID_forHomologs_dt,by="human_id_noVer",all=T)
#commonID_forHomologs_dt[,commonID:=.I]
#commonID_forHomologs_dt
#saveRDS(commonID_forHomologs_dt,"./R_rds_files/commonID_forHomologs_dt.rds")
#library(data.table)
#commonID_forHomologs_dt<-readRDS("./R_rds_files/commonID_forHomologs_dt.rds")

# step2: get gene_ids of differently spliced from the different datasets; translate it to commonID; do upset plot
#folders<-c("../HEK293T/","../actb_d1/","../bjFibro/","../c2c12/","../fetalLiv/","../hepg2/","../k562/","../kcytes/","../lncap/","../thymocytes/")
# removed c2c12 because only one rep: sketchy stats
folders<-c("../HEK293T/","../bjFibro/","../hepg2/","../k562/","../kcytes/","../lncap/","../thymocytes/","../fetalLiv/","../actb_d1/")
datasets<-basename(folders)

# ####################################################################################3
# # goal: upset plot with differentially expressed genes
# 
# # 1. get all diff expr events with padj<0.2 in one table
# get_expr_dt<-function(dataset){
#         # path to file
#         dataset_folder<-paste0("../",dataset)
#         file_path<-paste0(dataset_folder,"/rna_seq/hnrnpl/results/fcounts_deseq/normCounts.txt")
#         # get file
# 	require(data.table)
#         normCounts_dt<-fread(file_path)
# 	normCounts_dt[,dataset:=dataset]
# 
# 	# get info on what the control and hnRNPL_dep genotypes are called
# 	var_file_path<-paste0(dataset_folder,"/rna_seq/hnrnpl/data/pipeline_config_vars.sh")
# 	var_file<-readLines(var_file_path)
# 	ref_group_line<-var_file[grep("REFR_GROUP_NAME",var_file)]
# 	trt_group_line<-var_file[grep("TEST_GROUP_NAME",var_file)]
# 	ref_group_name<-sub("export REFR_GROUP_NAME=","",ref_group_line)
# 	trt_group_name<-sub("export TEST_GROUP_NAME=","",trt_group_line)
# 
# 	# get CTL and hnRNPL_dep columns
# 	normdt_colnames<-colnames(normCounts_dt)
# 	ref_colnames<-normdt_colnames[grep(ref_group_name,normdt_colnames)]
# 	trt_colnames<-normdt_colnames[grep(trt_group_name,normdt_colnames)]
# 
# 	# calculate means for each gene according to genotype
# 	normCounts_dt[,CTL_mean:=rowMeans(.SD),.SDcols=ref_colnames]
# 	normCounts_dt[,hnRNPL_dep_mean:=rowMeans(.SD),.SDcols=trt_colnames]
# 	out_dt<-normCounts_dt[,.SD,.SDcols=!c(ref_colnames,trt_colnames,"leg")]
# 
#         return(out_dt)
# }
# 
# allNormcounts<-sapply(datasets,get_expr_dt,simplify=F)
# allNormcounts<-rbindlist(allNormcounts)
# 
# # get commonID for each gene
# # melt commonID_forHomologs_dt to enable combining
# melted_commonIDs_dt<-melt(commonID_forHomologs_dt, id.vars=c("commonID"), measure.vars=c("mouse_id","human_id"), variable.name="mouse_or_human_id", value.name="GeneID", na.rm=T)
# melted_commonIDs_dt[,mouse_or_human:=sub('_id','',mouse_or_human_id),by=.I]
# melted_commonIDs_dt[,mouse_or_human_id:=NULL]
# 
# allNormcounts<-merge(melted_commonIDs_dt,allNormcounts,all.y=T,all.x=F,by.x="GeneID",by.y="gene_id")
# saveRDS(allNormcounts,"./R_rds_files/allNormcounts.rds")
require(data.table)
allNormcounts<-readRDS("./R_rds_files/allNormcounts.rds")

#filter here
allSig_DEG<-allNormcounts[padj<0.1 & abs(log2FoldChange)>log(1.5,2)]
allSig_DEG[,.(unique(GeneID)),dataset][,.N,dataset] #number of genes differentially expressed per dataset
#allSig_DEG<-allNormcounts[padj<0.2]

#making tables for the paper (common DGEs)
allSig_DEG[log2FoldChange < 0, downreg_in_n_datasets:=.N, by=commonID ]
allSig_DEG[log2FoldChange > 0, upreg_in_n_datasets:=.N, by=commonID ]
allSig_DEG[is.na(downreg_in_n_datasets), downreg_in_n_datasets:=0]
allSig_DEG[is.na(upreg_in_n_datasets), upreg_in_n_datasets:=0]

unique(allSig_DEG[mouse_or_human=="mouse" & log2FoldChange>0, .(status=if(.N==3){"up"}else{""}), by=commonID])[status=="up",.N]
unique(allSig_DEG[mouse_or_human=="mouse" & log2FoldChange<0, .(status=if(.N==3){"up"}else{""}), by=commonID])[status=="up",.N]

allSig_DEG[mouse_or_human=="mouse" & log2FoldChange>0, status:=if(.N==3){"Upreg in all mouse datasets"}else{""}, by=commonID]
allSig_DEG[mouse_or_human=="mouse" & log2FoldChange<0, status:=if(.N==3){"Downreg in all mouse datasets"}else{""}, by=commonID]
allSig_DEG[mouse_or_human=="human" & log2FoldChange>0, status:=if(.N==6){"Upreg in all human datasets"}else{""}, by=commonID]
allSig_DEG[mouse_or_human=="human" & log2FoldChange<0, status:=if(.N==6){"Downreg in all human datasets"}else{""}, by=commonID]

allSig_DEG[status=="Upreg in all mouse datasets" & upreg_in_n_datasets==3, status:="Upreg in all mouse datasets but not in any human_dataset"]
allSig_DEG[status=="Downreg in all mouse datasets" & downreg_in_n_datasets==3, status:="Downreg in all mouse datasets but not in any human dataset"]
allSig_DEG[status=="Upreg in all human datasets" & upreg_in_n_datasets==6, status:="Upreg in all human datasets but not in any mouse dataset"]
allSig_DEG[status=="Downreg in all human datasets" & downreg_in_n_datasets==6, status:="Downreg in all human datasets but not in any mouse dataset"]


allSig_DEG[status=="Upreg in all mouse datasets but not in any human_dataset"]
allSig_DEG[status=="Downreg in all mouse datasets but not in any human dataset"]
allSig_DEG[status=="Upreg in all human datasets but not in any mouse dataset"]
allSig_DEG[status=="Downreg in all human datasets but not in any mouse dataset"]


allSig_DEG<-allSig_DEG[order(commonID,mouse_or_human,dataset)]
allSig_DEG<-allSig_DEG[order(downreg_in_n_datasets, upreg_in_n_datasets, decreasing=T)]
write.table(allSig_DEG, "./R_rds_files/allSig_DEG_FDR_L2FC_filt.txt", sep="\t", quote=T, row.names=F)

#---------------------------------------------------------------------#
# 2. upset plot
# i. plotting all DGEs
# extract relevant info for upset plot
infoForUpset<-unique(allSig_DEG[,.(mouse_or_human,deg_count=length(unique(commonID))), by=dataset])
# order datasets in the order for display
infoForUpset<-infoForUpset[order(mouse_or_human,deg_count,decreasing = T)]

# make objects for Upset and also to intersection comparisons
DT_toUseForUpset<-allSig_DEG
DT_toUseForUpset[,present:=1]

# make matrix for upset plot
wideDT_forUpset<-dcast(DT_toUseForUpset, formula = commonID ~ dataset, value.var = "present", fill=0)
wideDT_forUpset<-unique(wideDT_forUpset)
mat_forUpset<-as.matrix(wideDT_forUpset,rownames="commonID")
#DGE_list<-sapply(datasets,function(x){allSig_DEG[dataset==x,unique(commonID)]})

# defining non-zero intersections to plot
dataset_combns<-unlist(lapply(9:6,function(x){combn(datasets,x,simplify = F)}),recursive = F)

in_upset_intersect<-function(dataset_combn){
  # check if the combination has a non-zero number of events
  n_in_cur_combn<-rowSums(mat_forUpset[,dataset_combn])
  n_in_all<-rowSums(mat_forUpset)
  # events present in all datasets of current combination and only in current combination
  rowWise_status_upset_intersect<- n_in_all==length(dataset_combn) & n_in_cur_combn==n_in_all
  return(rowWise_status_upset_intersect)
}

# filter combinations to include or not
combn_to_include<-sapply(dataset_combns, function(dataset_combn){
  # does the combination have at least 2 mouse and 2 human datasets?
  in2datasets_perSpecies<-min(infoForUpset[mouse_or_human=="mouse" & dataset %in% dataset_combn, .N],
                              infoForUpset[mouse_or_human=="human" & dataset %in% dataset_combn, .N]) >=1
  if(in2datasets_perSpecies){
    nonZero_combn<-sum(in_upset_intersect(dataset_combn))!=0
    return(nonZero_combn)
  } else{
    return(FALSE)
  }
})

combns_with_conserved_events<-dataset_combns[combn_to_include]

#specify the required interactions
combns_for_upset<-c(combns_with_conserved_events,
                    list(infoForUpset[mouse_or_human=="human",dataset]),
                    list(infoForUpset[mouse_or_human=="mouse",dataset]))

pdf(file="./figures/allDGE_upset_FDR_pt1_L2FC_inc.pdf",useDingbats=F,width=3,height=3)
library(UpSetR)
upset(as.data.frame(mat_forUpset), intersections = rev(combns_for_upset),
      sets=infoForUpset[,dataset], keep.order=T,
      text.scale=1, mb.ratio=c(0.55,0.45),
      point.size=1.5, line.size = 0.3)
dev.off()
#----------------------------------------------------------#
# ii. plotting downreg DGEs. keep same order as before
# make objects for Upset and also to intersection comparisons
DT_toUseForUpset<-allSig_DEG[log2FoldChange<0]
DT_toUseForUpset[,present:=1]
# make matrix for upset plot
wideDT_forUpset<-dcast(DT_toUseForUpset, formula = commonID ~ dataset, value.var = "present", fill=0)
wideDT_forUpset<-unique(wideDT_forUpset)
mat_forUpset<-as.matrix(wideDT_forUpset,rownames="commonID")

# defining non-zero intersections to plot
dataset_combns<-unlist(lapply(9:6,function(x){combn(datasets,x,simplify = F)}),recursive = F)

# filter combinations to include or not
combn_to_include<-sapply(dataset_combns, function(dataset_combn){
  # does the combination have at least 2 mouse and 2 human datasets?
  in2datasets_perSpecies<-min(infoForUpset[mouse_or_human=="mouse" & dataset %in% dataset_combn, .N],
                              infoForUpset[mouse_or_human=="human" & dataset %in% dataset_combn, .N]) >=1
  if(in2datasets_perSpecies){
    nonZero_combn<-sum(in_upset_intersect(dataset_combn))!=0
    return(nonZero_combn)
  } else{
    return(FALSE)
  }
})

combns_with_conserved_events<-dataset_combns[combn_to_include]

#specify the required interactions
combns_for_upset<-c(combns_with_conserved_events,
                    list(infoForUpset[mouse_or_human=="human",dataset]),
                    list(infoForUpset[mouse_or_human=="mouse",dataset]))

pdf(file="./figures/downreg_DGE_upset_FDR_pt1_L2FC_inc.pdf",useDingbats=F,width=3,height=3)
require(UpSetR)
upset(as.data.frame(mat_forUpset), intersections = rev(combns_for_upset),
      sets=infoForUpset[,dataset], keep.order=T,
      text.scale=1, mb.ratio=c(0.5,0.5),
      point.size=1.5, line.size = 0.3)
dev.off()

#----------------------------------------------------------#
# iii. plotting upreg DGEs. again keep orders

# make objects for Upset and also to intersection comparisons
DT_toUseForUpset<-allSig_DEG[log2FoldChange>0]
DT_toUseForUpset[,present:=1]
# make matrix for upset plot
wideDT_forUpset<-dcast(DT_toUseForUpset, formula = commonID ~ dataset, value.var = "present", fill=0)
wideDT_forUpset<-unique(wideDT_forUpset)
mat_forUpset<-as.matrix(wideDT_forUpset,rownames="commonID")

# defining non-zero intersections to plot
dataset_combns<-unlist(lapply(9:6,function(x){combn(datasets,x,simplify = F)}),recursive = F)

# filter combinations to include or not
combn_to_include<-sapply(dataset_combns, function(dataset_combn){
  # does the combination have at least 2 mouse and 2 human datasets?
  in2datasets_perSpecies<-min(infoForUpset[mouse_or_human=="mouse" & dataset %in% dataset_combn, .N],
                              infoForUpset[mouse_or_human=="human" & dataset %in% dataset_combn, .N]) >=1
  if(in2datasets_perSpecies){
    nonZero_combn<-sum(in_upset_intersect(dataset_combn))!=0
    return(nonZero_combn)
  } else{
    return(FALSE)
  }
})

combns_with_conserved_events<-dataset_combns[combn_to_include]

#specify the required interactions
combns_for_upset<-c(combns_with_conserved_events,
                    list(infoForUpset[mouse_or_human=="human",dataset]),
                    list(infoForUpset[mouse_or_human=="mouse",dataset]))

pdf(file="./figures/upreg_DGE_upset_FDR_pt1_L2FC_inc.pdf",useDingbats=F,width=3,height=3)
require(UpSetR)
upset(as.data.frame(mat_forUpset), intersections = rev(combns_for_upset),
      sets=infoForUpset[,dataset], keep.order=T,
      text.scale=1, mb.ratio=c(0.55,0.45),
      point.size=1.5, line.size = 0.3)
dev.off()

#---------------------------------------------------------------------#
# 3. make "inset" graphs to display continuous values for different datasets
# define commonName for plot labels
allSig_DEG[,commonName:=paste(unique(gene_name),collapse="|"),by=commonID]

# hits in at least 6 / 9 datasets
allSig_DEG[,commonID,by=dataset][,.N,by=commonID][N>=6]
common_hits<-allSig_DEG[,commonID,by=dataset][,.N,by=commonID][N>=6,commonID]
allSig_DEG[commonID%in%common_hits][order(commonID),.(dataset,commonID,gene_name,log2FoldChange,padj,CTL_mean,hnRNPL_dep_mean)]

# color values
colorsForDatasets<-c("#d7403c","#775fa9","#92ab3c","#29b34a","#2db99c","#0bb9e2","#6f94cc","#b37bb5","#d29329")
names(colorsForDatasets)<-c("actb_d1","thymocytes","fetalLiv","HEK293T","hepg2","k562","kcytes","lncap","bjFibro")

# exploring data to find consistently changed genes
p1<-ggplot(data=allSig_DEG[commonID%in%common_hits],aes(x=as.character(commonID),y=log2FoldChange, shape=mouse_or_human))+
  geom_jitter(width=0.1,aes(color=dataset),size=5)+
  geom_hline(yintercept=0)+
  scale_shape_manual(values=c(20,18))+
  scale_color_manual(values=colorsForDatasets)+
  scale_x_discrete(labels=common_hits_names)+
  theme_classic()+
  theme(text=element_text(family="sans"),
        plot.title=element_text(size=8,hjust=0.5),
        axis.line=element_line(color='black',linewidth=0.25),
        axis.ticks=element_line(colour="black",linewidth = 0.25),
        axis.title=element_text(size=7),
        axis.text.x=element_text(angle=90,size=8,vjust=0.5),
        axis.text.y=element_text(size=8,hjust=0.5))
p1

# downreg hits in at least 6 / 9 datasets - assign to the relevant commonIDs to common_hits
common_hits<-allSig_DEG[log2FoldChange<0,commonID,by=dataset][,.N,by=commonID][N>=6,commonID]
allSig_DEG[commonID%in%common_hits][order(commonID),.(dataset,commonID,gene_name,log2FoldChange,padj,CTL_mean,hnRNPL_dep_mean)]

pdf(file="./figures/downregHitsIn6of9datasets_FDR_L2FC_Filt.pdf.pdf",useDingbats=F,width=3,height=3)
require(ggplot2)
require(egg)
p1<-ggplot(data=allSig_DEG[commonID%in%common_hits],aes(x=as.character(commonID),y=log2FoldChange, shape=mouse_or_human))+
  geom_jitter(width=0.2,aes(color=dataset),size=2.4)+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(-4,1))+
  scale_shape_manual(values=c(20,18))+
  scale_color_manual(values=colorsForDatasets)+
  scale_x_discrete(labels=common_hits_names)+
  theme_classic()+
  theme(text=element_text(family="sans"),
        plot.title=element_text(size=8,hjust=0.5),
        axis.line=element_line(color='black',linewidth=0.25),
        axis.ticks=element_line(colour="black",linewidth = 0.25),
        axis.title=element_text(size=7),
        axis.text.x=element_text(angle=90,size=5,vjust=0.5),
        axis.text.y=element_text(size=5,hjust=0.5))
grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(0.8,"cm"),height=unit(1.75,"cm")))
dev.off()

#------------------------------------#
# upreg hits in at least 6 / 9 datasets - assign to the relevant commonIDs to common_hits
common_hits<-allSig_DEG[log2FoldChange>0,commonID,by=dataset][,.N,by=commonID][N>=6,commonID]
allSig_DEG[commonID%in%common_hits][order(commonID),.(dataset,commonID,gene_name,log2FoldChange,padj,CTL_mean,hnRNPL_dep_mean)]
# only ZBTB40 () is also present in one 

# to change x-axis labels from commonID to commonName
common_hits_names<-sapply(common_hits,function(x){allSig_DEG[commonID==x,unique(commonName)]})
names(common_hits_names)<-common_hits

pdf(file="./figures/upregHitsIn6of9datasets_FDR_L2FC_Filt.pdf",useDingbats=F,width=3,height=3)
require(ggplot2)
require(egg)
p1<-ggplot(data=allSig_DEG[commonID%in%common_hits],aes(x=as.character(commonID),y=log2FoldChange, shape=mouse_or_human))+
  geom_jitter(width=0.2,aes(color=dataset),size=2.4)+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(-1,4))+
  scale_shape_manual(values=c(20,18))+
  scale_color_manual(values=colorsForDatasets)+
  scale_x_discrete(labels=common_hits_names)+
  theme_classic()+
  theme(text=element_text(family="sans"),
        plot.title=element_text(size=8,hjust=0.5),
        axis.line=element_line(color='black',linewidth=0.25),
        axis.ticks=element_line(colour="black",linewidth = 0.25),
        axis.title=element_text(size=7),
        axis.text.x=element_text(angle=90,size=5,vjust=0.5),
        axis.text.y=element_text(size=5,hjust=0.5))
grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(2.8,"cm"),height=unit(1.75,"cm")))
dev.off()
#---------------------------------------------------------------------#
