library(data.table)
commonID_forHomologs_dt<-readRDS("./R_rds_files/commonID_forHomologs_dt.rds")

# step2: get gene_ids of differently spliced from the different datasets; translate it to commonID; do upset plot
#folders<-c("../HEK293T/","../actb_d1/","../bjFibro/","../c2c12/","../fetalLiv/","../hepg2/","../k562/","../kcytes/","../lncap/","../thymocytes/")
# removed c2c12 because only one rep: sketchy stats
folders<-c("../HEK293T/","../bjFibro/","../hepg2/","../k562/","../kcytes/","../lncap/","../thymocytes/","../fetalLiv/","../actb_d1/")
datasets<-basename(folders)
splicing_type<-c("SE")

####################################################################################
# goal: upset plot with splicing events not genes (with arbitrary filters)
# 1. load "merged_allSig_SE" data.table and calculate mean counts
merged_allSig_SE<-readRDS("./R_rds_files/merged_allSig_SE.rds")

merged_allSig_SE[, IncLev_max_diff:=mapply(function(IncL1,IncL2){
  IncLs<-c(as.numeric(unlist(strsplit(c(IncL1),","))),as.numeric(unlist(strsplit(c(IncL2),","))))
  max(IncLs)-min(IncLs)
}, IncL1=IncLevel1, IncL2=IncLevel2, USE.NAMES = F ),
by=.I]

merged_allSig_SE[, IJC_counts:=paste(IJC_SAMPLE_1,IJC_SAMPLE_2, sep=","), .I]
merged_allSig_SE[, SJC_counts:=paste(SJC_SAMPLE_1,SJC_SAMPLE_2, sep=","), .I]
merged_allSig_SE[, all_counts:=paste(IJC_SAMPLE_1,IJC_SAMPLE_2,SJC_SAMPLE_1,SJC_SAMPLE_2, sep=","), .I]

merged_allSig_SE[, IJC_mean:=sapply( IJC_counts, function(x){mean(as.numeric(unlist(strsplit(c(x),","))))} ), .I]
merged_allSig_SE[, SJC_mean:=sapply( SJC_counts, function(x){mean(as.numeric(unlist(strsplit(c(x),","))))} ), .I]
merged_allSig_SE[, mean_count:=sapply( all_counts, function(x){2 * mean(as.numeric(unlist(strsplit(c(x),","))))} ), .I]


merged_allSig_SE

merged_allSig_SE[grep("Sirt1", geneSymbol, ignore.case = T)]
merged_allSig_SE[commonSpliceID=="37776_1"]

####################################################################################
# 2. get junction specific counts from the "*/SE.MATS.JC.txt" files
get_junc_count_dt<-function(cur_dataset,alt_splice_tpye){
  # path to file
  dataset_folder<-paste0("../",cur_dataset)
  jc_file_path<-paste0(dataset_folder,"/rna_seq/hnrnpl/results/rMATS/",alt_splice_tpye,".MATS.JC.txt")
  
  jc_rmats_dt<-fread(jc_file_path)
  # remove one of the duplicated ID columns
  jc_rmats_dt[,ID:=NULL]
  # filter by ID to get only relevant rows
  jc_rmats_dt<-jc_rmats_dt[ID %in% merged_allSig_SE[dataset==cur_dataset,ID]]
  
  jc_rmats_dt[, dataset:=cur_dataset]
  # calculate means of counts
  jc_rmats_dt[, junc_IJC_counts:=paste(IJC_SAMPLE_1,IJC_SAMPLE_2, sep=","), .I]
  jc_rmats_dt[, junc_SJC_counts:=paste(SJC_SAMPLE_1,SJC_SAMPLE_2, sep=","), .I]
  jc_rmats_dt[, junc_all_counts:=paste(IJC_SAMPLE_1,IJC_SAMPLE_2,SJC_SAMPLE_1,SJC_SAMPLE_2, sep=","), .I]
  
  jc_rmats_dt[, junc_IJC_mean:=sapply( junc_IJC_counts, function(x){mean(as.numeric(unlist(strsplit(c(x),","))))} ), .I]
  jc_rmats_dt[, junc_SJC_mean:=sapply( junc_SJC_counts, function(x){mean(as.numeric(unlist(strsplit(c(x),","))))} ), .I]
  jc_rmats_dt[, junc_mean_count:=sapply( junc_all_counts, function(x){2 * mean(as.numeric(unlist(strsplit(c(x),","))))} ), .I]
  
  jc_count_dt<-jc_rmats_dt[,.(ID, dataset, junc_IJC_counts, junc_SJC_counts, junc_all_counts,
                              junc_IJC_mean, junc_SJC_mean, junc_mean_count)]
  
  return(jc_count_dt)
}

jc_allSig_SE<-sapply(datasets,get_junc_count_dt,"SE",simplify=F)
jc_allSig_SE<-rbindlist(jc_allSig_SE)

merged_allSig_SE<-merge(merged_allSig_SE, jc_allSig_SE, by=c("ID","dataset"), all.x=T)
merged_allSig_SE[is.na(junc_IJC_mean), junc_IJC_mean:=0]
merged_allSig_SE[is.na(junc_SJC_mean), junc_SJC_mean:=0]
merged_allSig_SE[is.na(junc_mean_count), junc_mean_count:=0]

saveRDS(merged_allSig_SE, "./R_rds_files/merged_allSig_SE_withJunctionReads.rds")

merged_allSig_SE[grep("Sirt1", geneSymbol, ignore.case = T)]
merged_allSig_SE[commonSpliceID=="37776_1"]

# Filter!
filt_merged_allSig_SE<-merged_allSig_SE[junc_mean_count>=10 & FDR<0.1 & abs(IncLevelDifference)>=0.1 & abs(IncLevelDifference)<=0.95 & IncLev_max_diff>0.05]

filt_merged_allSig_SE[grep("NA", IncLevel1)]
filt_merged_allSig_SE[grep("NA", IncLevel2)]

filt_merged_allSig_SE[,.N,commonSpliceID][N>=6][order(N,commonSpliceID)]
cat(sort(filt_merged_allSig_SE[commonID %in% filt_merged_allSig_SE[,.N,by=c("commonID","commonSpliceID")][N>=6][,commonID]][, unique(toupper(geneSymbol))]) ,sep="\n")

filt_merged_allSig_SE[,.N,commonSpliceID][N>=7][order(N,commonSpliceID)]
cat(sort(filt_merged_allSig_SE[commonID %in% filt_merged_allSig_SE[,.N,by=c("commonID","commonSpliceID")][N>=7][,commonID]][, unique(toupper(geneSymbol))]) ,sep="\n")

filt_merged_allSig_SE[,.N,commonSpliceID][N>=8][order(N,commonSpliceID)]
cat(sort(filt_merged_allSig_SE[commonID %in% filt_merged_allSig_SE[,.N,by=c("commonID","commonSpliceID")][N>=8][,commonID]][, unique(toupper(geneSymbol))]) ,sep="\n")

#make tables for the paper
filt_merged_allSig_SE[,alt_spliced_in_n_datasets:=.N,by=commonSpliceID]
filt_merged_allSig_SE[mouse_or_human=="mouse", status:=if(.N==3){"Alt spliced in all mouse datasets"}else{""}, by=commonSpliceID]
filt_merged_allSig_SE[mouse_or_human=="human", status:=if(.N==6){"Alt spliced in all human datasets"}else{""}, by=commonSpliceID]

filt_merged_allSig_SE[status=="Alt spliced in all mouse datasets" & alt_spliced_in_n_datasets==3, status:="Alt spliced in all mouse datasets but not in any human dataset"]
filt_merged_allSig_SE[status=="Alt spliced in all human datasets" & alt_spliced_in_n_datasets==6, status:="Alt spliced in all human datasets but not in any mouse dataset"]


filt_merged_allSig_SE[status=="Alt spliced in all mouse datasets but not in any human dataset"]
filt_merged_allSig_SE[status=="Alt spliced in all human datasets but not in any mouse dataset"]

filt_merged_allSig_SE<-filt_merged_allSig_SE[order(commonSpliceID,mouse_or_human,dataset)]
filt_merged_allSig_SE<-filt_merged_allSig_SE[order(alt_spliced_in_n_datasets, decreasing=T)]
#write.table(filt_merged_allSig_SE, "./R_rds_files/allSig_SE_DSE_datatable_filt.txt", sep=",", quote=T, row.names=F)

# see individual events
filt_merged_allSig_SE[commonSpliceID=="41568_1"] #Baz2b
filt_merged_allSig_SE[commonSpliceID=="39483_2"] #NSD2
filt_merged_allSig_SE[commonSpliceID=="39483_3"] #NSD2
filt_merged_allSig_SE[commonSpliceID=="37776_1"] #Sirt1

####################################################################################3
# make objects for Upset and also to intersection comparisons
DT_toUseForUpset<-filt_merged_allSig_SE
DT_toUseForUpset[,present:=1]

# make matrix for upset plot
wideDT_forUpset<-dcast(DT_toUseForUpset, formula = commonSpliceID ~ dataset, value.var = "present", fill=0)
wideDT_forUpset<-unique(wideDT_forUpset)
mat_forUpset<-as.matrix(wideDT_forUpset,rownames="commonSpliceID")

infoForUpset<-unique(DT_toUseForUpset[,.(mouse_or_human,spl_count=length(unique(commonSpliceID))), by=dataset])
# order datasets in the order for display
infoForUpset<-infoForUpset[order(mouse_or_human,spl_count,decreasing = T)]

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

# get conserved events (commonSpliceIDs)
conserved_events<-lapply(combns_with_conserved_events, function(dataset_combn){
  rownames(mat_forUpset)[in_upset_intersect(dataset_combn)]
})
names(conserved_events)<-sapply(combns_with_conserved_events, paste0, collapse="|")
conserved_events

lapply(conserved_events,function(commonSpliceIDs){
  DT_toUseForUpset[commonSpliceID %in% commonSpliceIDs,unique(toupper(geneSymbol))]})

conserved_event_gene_commonIDs<-sort(unique(sapply(unlist(conserved_events,use.names = F),
       function(x) {strsplit(x, "_")[[1]][1]}, USE.NAMES = F)))
paste(length(unlist(conserved_events)), "consevered events from", length(conserved_event_gene_commonIDs), "genes")

DT_toUseForUpset[commonID %in% conserved_event_gene_commonIDs, sort(unique(toupper(geneSymbol)))]
cat( DT_toUseForUpset[commonID %in% conserved_event_gene_commonIDs, sort(unique(toupper(geneSymbol)))], sep="\n")

combns_for_upset<-c(combns_with_conserved_events,
                  list(infoForUpset[mouse_or_human=="human",dataset]),
                  list(infoForUpset[mouse_or_human=="mouse",dataset]))

lapply(combns_for_upset, function(dataset_combn){
  commonSpliceIDs<-rownames(mat_forUpset)[in_upset_intersect(dataset_combn)]
  DT_toUseForUpset[commonSpliceID %in% commonSpliceIDs,unique(toupper(geneSymbol))]
})

# calculate upset counts
counts_in_combn<-lapply(combns_for_upset, function(dataset_combn){
  n_in_cur_combn<-rowSums(mat_forUpset[,dataset_combn])
  n_in_all<-rowSums(mat_forUpset)
  # events present in all datasets of current combination and only in current combination
  status_cur_combn<- n_in_all==length(dataset_combn) & n_in_cur_combn==n_in_all
  # number of events present in only the current combination
  sum(status_cur_combn)
})
names(counts_in_combn)<-sapply(combns_for_upset, paste0, collapse = "|")
unlist(counts_in_combn)

df_forUpset<-as.data.frame(mat_forUpset)[infoForUpset[,dataset]]

# do Upset plot
library(UpSetR)
pdf(file="./figures/splice_SE_upset_FDRpt1_IncPt1_cov10_filt.pdf",useDingbats=F,width=3.5,height=3)
upset(df_forUpset, intersections = rev(combns_for_upset),
      sets=infoForUpset[,dataset], keep.order=T,
      text.scale=1, mb.ratio=c(0.55,0.45),
      point.size=1.5, line.size = 0.3)
dev.off()
#------------------------------#
# get FASTA sequences for events in the upset plot (conserved events + species-specific conserved events)
# will take sequences from the beginning of the upstream exon to the end of the downstream exon
# using the above code for doing UpSet plot to get event IDs

# 1. get inter-species conserved events (see above for definition)
events_for_FASTA<-filt_merged_allSig_SE[commonSpliceID %in% unlist(conserved_events)]
# 2. get intra-species conserved events 
intra_sp_combns<-c(human=list(infoForUpset[mouse_or_human=="human",dataset]),
                  mouse=list(infoForUpset[mouse_or_human=="mouse",dataset]))

intra_sp_events<-lapply(intra_sp_combns, function(dataset_combn){
  rownames(mat_forUpset)[in_upset_intersect(dataset_combn)]
})

#coerce to data.table
intra_sp_event_dt<-data.table(data.frame(commonSpliceID=unlist(intra_sp_events)), keep.rownames = "mouse_or_human")
intra_sp_event_dt[,mouse_or_human:=gsub("[0-9]","",mouse_or_human)]

#add intra-species events to the events_for_FASTA data.table
events_for_FASTA<-rbind(events_for_FASTA, 
                        merge(intra_sp_event_dt, filt_merged_allSig_SE, all.x=T, all.y=F))

unique(events_for_FASTA[,.(commonSpliceID, mouse_or_human, chr, upstreamES, downstreamEE, geneSymbol, strand)])[order(chr,upstreamES)]
# several events have the same upstreamES and downstreamEE and these duplicates can be removed
# These events are in NSD2 (2 events), Nsd2 (2 events), Vps13d (3), Fgfr1op2 (3) and Sec31a (2)
coords_for_cons_events<-unique(events_for_FASTA[,.(mouse_or_human, chr, upstreamES, downstreamEE, geneSymbol, strand)])[order(chr,upstreamES)]
coords_for_cons_events[,fasta_name:=paste0(geneSymbol,"_ev",1:.N), by=geneSymbol]

# get FASTA files
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)

humanGR<-makeGRangesFromDataFrame(coords_for_cons_events[mouse_or_human=="human"], keep.extra.columns = T,
                                  start.field="upstreamES", end.field="downstreamEE")
mouseGR<-makeGRangesFromDataFrame(coords_for_cons_events[mouse_or_human=="mouse"], keep.extra.columns = T,
                                  start.field="upstreamES", end.field="downstreamEE")

human_cons_seqs<-getSeq(Hsapiens, humanGR)
mouse_cons_seqs<-getSeq(Hsapiens, mouseGR)

names(human_cons_seqs)<-humanGR$fasta_name
names(mouse_cons_seqs)<-mouseGR$fasta_name

dir.create("fasta_files_cons_seq")
writeXStringSet(human_cons_seqs, "./fasta_files_cons_seq/cons_seqs.fa")
writeXStringSet(mouse_cons_seqs, "./fasta_files_cons_seq/cons_seqs.fa", append = T)

#------------------------------#
# motifs
library(strinr)
str_count(as.character(human_cons_seqs),"(?=(ACACACA|ACACAAA))")
str_count(as.character(mouse_cons_seqs),"(?=(ACACACA|ACACAAA))")

sum(str_count(as.character(human_cons_seqs),"(?=(ACACACA|ACACAAA))")==0) #9
sum(str_count(as.character(mouse_cons_seqs),"(?=(ACACACA|ACACAAA))")==0) #16

sum(str_count(as.character(human_cons_seqs),"(?=(ACAC|CACA).{0,4}(ACAC|CACA))")==0) #7
sum(str_count(as.character(mouse_cons_seqs),"(?=(ACAC|CACA).{0,4}(ACAC|CACA))")==0) #11

sum(str_count(as.character(human_cons_seqs),"(?=(ACAC|CACA|CAAA).{0,4}(ACAC|CACA|CAAA))")==0) #3
sum(str_count(as.character(mouse_cons_seqs),"(?=(ACAC|CACA|CAAA).{0,4}(ACAC|CACA|CAAA))")==0) #3

sum(str_count(as.character(human_cons_seqs), "(?=CA.{0,2}CA.{0,2}CA.{0,2}CA)")==0) #3
sum(str_count(as.character(mouse_cons_seqs), "(?=CA.{0,2}CA.{0,2}CA.{0,2}CA)")==0) #6
sum(str_count(as.character(human_cons_seqs), "(?=AC.{0,2}AC.{0,2}AC.{0,2}AC)")==0) #10
sum(str_count(as.character(mouse_cons_seqs), "(?=AC.{0,2}AC.{0,2}AC.{0,2}AC)")==0) #11

sum(str_count(as.character(human_cons_seqs), "(?=CA.{0,2}CA.{0,2}CA)")==0) #0
sum(str_count(as.character(mouse_cons_seqs), "(?=CA.{0,2}CA.{0,2}CA)")==0) #0
sum(str_count(as.character(mouse_cons_seqs), "(?=AC.{0,2}AC.{0,2}AC)")==0) #0
sum(str_count(as.character(human_cons_seqs), "(?=AC.{0,2}AC.{0,2}AC)")==0) #0

#------------------------------#
# get sequence for skipped exon only
coords_for_cons_SEs<-unique(events_for_FASTA[,.(mouse_or_human, chr, exonStart_0base, exonEnd, geneSymbol, strand)])[order(chr,exonStart_0base)]
coords_for_cons_SEs[,fasta_name:=paste0(geneSymbol,"_ev",1:.N), by=geneSymbol]

humanGR<-makeGRangesFromDataFrame(coords_for_cons_SEs[mouse_or_human=="human"], keep.extra.columns = T,
                                  start.field="exonStart_0base", end.field="exonEnd")
mouseGR<-makeGRangesFromDataFrame(coords_for_cons_SEs[mouse_or_human=="mouse"], keep.extra.columns = T,
                                  start.field="exonStart_0base", end.field="exonEnd")
human_cons_seqs<-getSeq(Hsapiens, humanGR)
mouse_cons_seqs<-getSeq(Hsapiens, mouseGR)
names(human_cons_seqs)<-humanGR$fasta_name
names(mouse_cons_seqs)<-mouseGR$fasta_name
writeXStringSet(human_cons_seqs, "./fasta_files_cons_seq/cons_SEonly_seqs.fa")
writeXStringSet(mouse_cons_seqs, "./fasta_files_cons_seq/cons_SEonly_seqs.fa", append = T)


#------------------------------#
# getting and plotting exon inclusion levels

merged_allSig_SE<-readRDS("./R_rds_files/merged_allSig_SE_withJunctionReads.rds")

get_incLevel<-function(cur_commonSpliceID){
  cur_dt<-merged_allSig_SE[commonSpliceID==cur_commonSpliceID]
  cur_dt[,common_name:=paste(unique(geneSymbol),collapse="|")]
  cur_dt[,CTL:=lapply(strsplit(IncLevel1,","), function(x){mean(as.numeric(x))}),by=dataset]
  cur_dt[,hnRNPL_dep:=lapply(strsplit(IncLevel2,","), function(x){mean(as.numeric(x))}),by=dataset]
  cur_dt<-melt(cur_dt, id.vars=c("dataset","mouse_or_human","common_name"), measure.vars=c("CTL","hnRNPL_dep"), variable.name="genotype", value.name="IncLevel")
  return(cur_dt)
}

# color values
colorsForDatasets<-c("#d7403c","#775fa9","#92ab3c","#29b34a","#2db99c","#0bb9e2","#6f94cc","#b37bb5","#d29329")
names(colorsForDatasets)<-c("actb_d1","thymocytes","fetalLiv","HEK293T","hepg2","k562","kcytes","lncap","bjFibro")

# make function for plotting
plot_incLevel<-function(cur_commonSpliceID, saveToPDF=FALSE, geneNameForPDF=c()){
  cur_dt<-get_incLevel(cur_commonSpliceID)
  plot_name<-cur_dt[,unique(common_name)]
  require(ggplot2)
  require(egg)
  p1<-ggplot(data=cur_dt, aes(x=genotype, y=IncLevel*100, group=dataset, color=dataset, shape=mouse_or_human))+
    geom_point()+
    geom_line(linewidth=0.25)+
    ggtitle(plot_name)+
    scale_y_continuous("Exon inclusion (%)")+
    scale_shape_manual(values=c(20,18))+
    scale_color_manual(values=colorsForDatasets)+
    theme_classic()+
    theme(text=element_text(family="sans"),
          plot.title=element_text(size=8),
          axis.line=element_line(color='black',linewidth=0.125),
          axis.ticks=element_line(colour="black",linewidth = 0.125),
          axis.title=element_text(size=7),
          axis.text.y=element_text(size=5),
          axis.text.x=element_text(size=5),
          strip.background=element_rect(linewidth=0.25),
          strip.text=element_text(size=7),
          panel.border=element_rect(fill=NA,linewidth=0.25),
          panel.spacing=unit(0,"line"))
  grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(1.5,"cm"),height=unit(2.5,"cm")))
  
  if(saveToPDF==T){
    output_folder<-"./figures/"
    if(is.null(geneNameForPDF)){
      pdf_name<-paste0(output_folder, "IncLevels_", sub("\\|.*","",plot_name), ".pdf")
    }else{
      pdf_name<-paste0(output_folder, "IncLevels_", geneNameForPDF, ".pdf")
    }
    pdf(pdf_name,useDingbats=F,height=5,width=5)
    grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(1.1,"cm"),height=unit(2.5,"cm")))
    dev.off()
  }
}

# conserved in 6/9 datasets at least 2 mouse datasets
plot_incLevel("36496_3") #GPBP1
plot_incLevel("45369_2") #KDM6A
plot_incLevel("41936_3") #HNRNPR
plot_incLevel("39483_2") #NSD2
plot_incLevel("39483_3") #NSD2 - event 2
plot_incLevel("39758_1") #FGFR1OP2
plot_incLevel("44005_2") #SEC13A
plot_incLevel("44005_3") #SEC13A
plot_incLevel("36441_2") #PPP1R12A
plot_incLevel("39236_1") #PPP3CB
plot_incLevel("41099_1") #EPC1
plot_incLevel("37776_1") #SIRT1
plot_incLevel("43316_1") #MAP3K7
plot_incLevel("35825_6") #LUC7L

# conserved in 6/9 datasets at least 1 mouse dataset
plot_incLevel("41936_2") #HNRNPR - event 2
plot_incLevel("41568_1") #BAZ2B
plot_incLevel("37427_3") #GNAS
plot_incLevel("44864_6") #UBAP2L


