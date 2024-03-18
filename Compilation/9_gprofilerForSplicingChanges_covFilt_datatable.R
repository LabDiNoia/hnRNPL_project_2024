library(data.table)

folders <- c("../HEK293T/","../bjFibro/","../hepg2/","../k562/","../kcytes/","../lncap/","../thymocytes/","../fetalLiv/","../actb_d1/")
datasets <- basename(folders)
mouse_or_human_dt <-  data.table(dataset = datasets, mouse_or_human = c(rep("human", 6), rep("mouse", 3)))

#cur_dataset<-datasets[1]
info_dt <- fread("refrFiles/datasets_readLen_SEPE_info.txt")

#-------------------------------------------------------------------------------#

num_from_commaSepStr<-function(commaSepStr){
  vec<-unlist(strsplit(c(commaSepStr),","))
  vec[vec=="NA"]<-NA
  vec<-as.numeric(vec)
  vec
}

allSplice_dt<-data.table()
for(cur_dataset in datasets){
  cur_rMATS_folder_path<-paste0("../", cur_dataset, "/rna_seq/hnrnpl/results/rMATS/")
  splicing_files<-list.files(cur_rMATS_folder_path, pattern="*.MATS.JCEC.txt", full.names = T)
  names(splicing_files)<-sub(".MATS.JCEC.txt", "", basename(splicing_files))

  # get all splicing events for all splicing types and all datasets
  cur_splice_dt<-lapply(splicing_files,
                       function(file_path){
                         require(data.table)
                         cur_splicing_type<-sub(".MATS.JCEC.txt", "", basename(file_path))
                         print(paste("Getting", cur_splicing_type,"splicing events for", cur_dataset, "dataset..."))
                         jcec_dt<-fread(file_path)
                         jcec_dt[,ID:=NULL] #remove duplicate ID column
                         jcec_dt[,splicing_type:=cur_splicing_type]
                         jcec_dt[, IncLev_max_diff:=mapply(function(IncL1,IncL2){
                           IncLs<-c(num_from_commaSepStr(IncL1), num_from_commaSepStr(IncL2))
                           return(max(IncLs)-min(IncLs))
                         }, IncL1=IncLevel1, IncL2=IncLevel2, USE.NAMES = F ),
                         by=.I]
                         jcec_dt[is.na(IncLev_max_diff),IncLev_max_diff:=0]
                         
                         
                         #jcec_dt<-jcec_dt[FDR<0.1 & abs(IncLevelDifference)>0.1]
                         
                         jc_file_path<-sub(".MATS.JCEC.txt", ".MATS.JC.txt", file_path)
                         jc_dt<-fread(jc_file_path)
                         jc_dt[,ID:=NULL] #remove duplicate ID column
                         jc_dt[, junc_mean_counts:=mapply(function(IJC_s1,IJC_s2,SJC_s1,SJC_s2){
                           counts<-c(num_from_commaSepStr(IJC_s1), num_from_commaSepStr(IJC_s2),
                                     num_from_commaSepStr(SJC_s1), num_from_commaSepStr(SJC_s2))
                           mean(counts)*2
                         }, IJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_1, SJC_SAMPLE_2, USE.NAMES = F),
                         by=.I]
                         jcec_dt<-merge(jcec_dt, jc_dt[,.(ID, junc_mean_counts)], all.x=T)
                         jcec_dt[is.na(junc_mean_counts),junc_mean_counts:=0]
                         jcec_dt[,.(splicing_type, ID, GeneID, geneSymbol,
                                    FDR, IncLevelDifference, IncLev_max_diff, junc_mean_counts)]
                       })
  cur_splice_dt<-rbindlist(cur_splice_dt)
  cur_splice_dt[,dataset:=cur_dataset]
  allSplice_dt<-rbind(allSplice_dt, cur_splice_dt)
}
allSplice_dt
allSplice_dt<-merge(mouse_or_human_dt,allSplice_dt)
saveRDS(allSplice_dt, "./R_rds_files/allSplice_dt.rds")
allSplice_dt<-readRDS("./R_rds_files/allSplice_dt.rds")


do_GOBP_gprofiler<-function(filteredSplice_dt){
  filteredSplice_dt[,gene_id_noVer:=sub("\\.[0-9]*","",GeneID)]
  datasets<-filteredSplice_dt[,unique(dataset)]
  lapply(datasets, function(cur_dataset){
    cur_splice_dt<-filt_allSplice_dt[dataset==cur_dataset]
    cat("Doing GO analysis for", length(cur_splice_dt[,unique(gene_id_noVer)]),
        "alternatively spliced (all) genes from the", cur_dataset, "dataset\n")
    
    #do GO
    require(gprofiler2)
    #same database as in August for reproducibility
    set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e109_eg56_p17")
    if(unique(cur_splice_dt[,mouse_or_human])=="human"){
      org="hsapiens"
    }else if(unique(cur_splice_dt[,mouse_or_human])=="mouse"){
      org="mmusculus"
    }
    
    gp_res<-gost(unique(cur_splice_dt[,gene_id_noVer]), organism = org, user_threshold = 0.05,
                 sources = "GO:BP", evcodes = T)  
    gp_dt<-data.table(gp_res$result)
    gp_dt<-gp_dt[term_size>=20 & term_size<=500]
    cat("Found", nrow(gp_dt) ,"terms enrinched.\n")
    gp_dt[,dataset:=cur_dataset]
    gp_dt
  })
}

# FILTERS
# coverage filter
juncMeanCount_cov_filt<-10
# FDR and incLevel filters
FDR_filt<-0.1
abs_IncLevelDiff_filt<-0.1
IncLevMaxDiff_filt<-0.05


filt_allSplice_dt<-allSplice_dt[junc_mean_counts>=juncMeanCount_cov_filt & 
                                  FDR<FDR_filt & 
                                  abs(IncLevelDifference)>=abs_IncLevelDiff_filt & 
                                  abs(IncLevelDifference)<=0.95 &
                                  IncLev_max_diff>IncLevMaxDiff_filt]

filt_allSplice_dt[dataset=="actb_d1"][,.N,splicing_type]
filt_allSplice_dt[dataset=="actb_d1"][IncLevelDifference<0,.N,splicing_type]
filt_allSplice_dt[dataset=="actb_d1"][IncLevelDifference>0,.N,splicing_type]
filt_allSplice_dt[dataset=="actb_d1",.(genes=unique(GeneID)),splicing_type][,.N,splicing_type]


commonID_forHomologs_dt<-readRDS("./R_rds_files/commonID_forHomologs_dt.rds")
commonID_forHomologs_dt<-melt(commonID_forHomologs_dt, id.vars = "commonID" , value.name = "GeneID_noVer",
     variable.name = "mouse_or_human_id", measure.vars = c("mouse_id_noVer", "human_id_noVer"))

filt_allSplice_dt[,GeneID_noVer:=sub("\\.[0-9]*","",GeneID)]
filt_allSplice_dt<-merge(filt_allSplice_dt, commonID_forHomologs_dt[,.(GeneID_noVer,commonID)], by="GeneID_noVer", all.x=T)

filt_allSplice_dt

common_spliceGenes_dt<-filt_allSplice_dt[,.(commonID=unique(commonID)),by=dataset][,.(in_x_datasets=.N),commonID]
common_spliceGenes_dt<-common_spliceGenes_dt[order(in_x_datasets, decreasing = T)]

filt_allSplice_dt[commonID %in% common_spliceGenes_dt[in_x_datasets>=8,commonID], sort(unique(toupper(geneSymbol)))]




gp_covfilt10<-do_GOBP_gprofiler(filt_allSplice_dt)
gp_covfilt10<-rbindlist(gp_covfilt10)
gp_covfilt10[dataset=="actb_d1",1:12]
gp_covfilt10[dataset=="actb_d1",1:12][1:12]
gp_covfilt10[,.N,by=term_name][N>=7][order(N,decreasing=T)]
cat(gp_covfilt10[,.N,by=term_id][N>=7][order(N,decreasing=T),term_id], sep="\n") #for REVIGO to cluster terms

gp_covfilt10[term_name=="histone modification",c(1:12,17)]


new_dataset_order<-c("actb_d1","fetalLiv","thymocytes","kcytes","HEK293T","k562","hepg2","lncap","bjFibro")
gp_covfilt10$dataset<-factor(gp_covfilt10$dataset, levels = new_dataset_order)

gp_covfilt10[, num_datasets := .N, by = term_id]
terms_to_plot<-unique(gp_covfilt10[num_datasets>=7, .(term_name, term_id, num_datasets)][order(num_datasets, decreasing = T)])

pval_dt<-gp_covfilt10[, .SD[match(terms_to_plot$term_id, term_id), .(term_name, nL10p=-log(p_value,10))], by = dataset]
pval_dt<-dcast(na.omit(pval_dt), term_name~dataset, value.var = "nL10p")[match(terms_to_plot$term_name, term_name)]
pval_dt<-pval_dt[match(term_name, terms_to_plot$term_name)]

write.table(pval_dt, "./R_rds_files/GOBP_forHeatmaps_filtCov10_datatable.txt", sep="\t", row.names = F, quote = F)


# FILTERS
# NO coverage filter - for showing reviewers
juncMeanCount_cov_filt<-0
# FDR and incLevel filters
FDR_filt<-0.1
abs_IncLevelDiff_filt<-0.1
IncLevMaxDiff_filt<-0


filt_allSplice_dt<-allSplice_dt[junc_mean_counts>=juncMeanCount_cov_filt & 
                                  FDR<FDR_filt & 
                                  abs(IncLevelDifference)>=abs_IncLevelDiff_filt & 
                                  abs(IncLevelDifference)<=0.95 &
                                  IncLev_max_diff>IncLevMaxDiff_filt]

filt_allSplice_dt
filt_allSplice_dt[dataset=="actb_d1"][,.N,splicing_type]
filt_allSplice_dt[dataset=="actb_d1"][IncLevelDifference<0,.N,splicing_type]
filt_allSplice_dt[dataset=="actb_d1"][IncLevelDifference>0,.N,splicing_type]
filt_allSplice_dt[dataset=="actb_d1",.(genes=unique(GeneID)),splicing_type][,.N,splicing_type]

gp_covfilt0<-do_GOBP_gprofiler(filt_allSplice_dt)
gp_covfilt0<-rbindlist(gp_covfilt0)
gp_covfilt0[dataset=="actb_d1",1:12]
gp_covfilt0[dataset=="actb_d1",1:12][1:12]
# gp_covfilt0[,.N,by=term_name][N>=7][order(N,decreasing=T)]
# cat(gp_covfilt0[,.N,by=term_id][N>=7][order(N,decreasing=T),term_id], sep="\n") #for REVIGO to cluster terms
pval_dt_cov0<-gp_covfilt0[, .SD[match(terms_to_plot$term_id, term_id), .(term_name, nL10p=-log(p_value,10))], by = dataset]
pval_dt_cov0<-dcast(na.omit(pval_dt_cov0), term_name~dataset, value.var = "nL10p")[match(terms_to_plot$term_name, term_name)]
pval_dt_cov0<-pval_dt_cov0[match(term_name, terms_to_plot$term_name)]
pval_dt_cov0<-pval_dt_cov0[,c("term_name", new_dataset_order),with=F]
write.table(pval_dt_cov0, "./R_rds_files/GOBP_forHeatmaps_NOfiltCov_datatable.txt", sep="\t", row.names = F, quote = F)


#-----------------------------------------------------------------#
# plotting trends - depth, read length, splicing changes

info_dt <- fread("refrFiles/datasets_readLen_SEPE_info.txt")
all_mqc_star <- data.table()
for (cur_dataset in datasets) {
  cur_mqc_star_file_path <- paste0("../",cur_dataset,"/rna_seq/hnrnpl/results/qc_stats/multiQC_report_data/multiqc_star.txt")
  mqc_star <- fread(cur_mqc_star_file_path)
  mqc_star <- mqc_star[grep("2pass", Sample)]
  cur_mqc_star <- mqc_star[, .(Sample, uniquely_mapped)]
  cur_mqc_star[, dataset := cur_dataset]
  all_mqc_star <- rbind(all_mqc_star, cur_mqc_star)
}

all_mqc_star
all_mqc_star[, .(median_counts = median(uniquely_mapped)), dataset]
all_mqc_star[, .(mean_mapped_counts = mean(uniquely_mapped)), dataset]

info_dt <- merge(info_dt, all_mqc_star[, .(mean_mapped_counts = mean(uniquely_mapped)), dataset])
info_dt$Description <- c("HEK293T", "Activated B", "BJ", "Fetal liver", "HepG2", "K562", "Keratinocytes", "LNCaP", "Thymocytes")

# calculate number of splicing changes
allSplice_dt<-readRDS("./R_rds_files/allSplice_dt.rds")

#coverage filter of 0
filt_allSplice_dt<-allSplice_dt[junc_mean_counts>=0 & 
                                  FDR<0.1 & 
                                  abs(IncLevelDifference)>=0.1 & 
                                  abs(IncLevelDifference)<=0.95 &
                                  IncLev_max_diff>0.05]
info_dt<-merge(info_dt, filt_allSplice_dt[ ,.(n_in_covFilt0=.N, ngenes_in_covFilt0=length(unique(GeneID))),dataset])

#coverage filter of 10
filt_allSplice_dt<-allSplice_dt[junc_mean_counts>=10 & 
                                  FDR<0.1 & 
                                  abs(IncLevelDifference)>=0.1 & 
                                  abs(IncLevelDifference)<=0.95 &
                                  IncLev_max_diff>0.05]

info_dt<-merge(info_dt, filt_allSplice_dt[ ,.(n_in_covFilt10=.N, ngenes_in_covFilt10=length(unique(GeneID))),dataset])

#write.table(info_dt, "R_rds_files/splicingEventCounts_vs_factors_datatableFilt.txt", sep = "\t", row.names = F, quote = F)

#----------------------------------#
# numbers for table S5. Listing number of changed splicing events and genes with splicing changes for each dataset
dcast(filt_allSplice_dt[,.N,by=c("splicing_type","dataset")], splicing_type ~ dataset, value.var = "N")
dcast(unique(filt_allSplice_dt[,.(GeneID,dataset,splicing_type)])[,.N,by=c("splicing_type","dataset")],
      splicing_type ~ dataset, value.var = "N")

#----------------------------------#
# plotting
for_fig_dt <-melt(info_dt, variable.name="cov_filter", value.name="event_counts", id.vars=1:7, measure.vars=c(8,10))

#pdf("figures/splicingEventCounts_vs_factors_datatableFilt.pdf",useDingbats=F,height=8,width=8)
library(ggplot2)
library(egg)
require(ggrepel)
p1 <- ggplot(for_fig_dt, aes(mean_mapped_counts, event_counts, shape = SE_or_PE, color = as.factor(read_length))) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 20000), name = "Number of splicing changes") +
  scale_x_continuous(limits = c(0, 100000000), name = "Averarge mapped read counts") +
  scale_color_manual(values = c("light grey", "grey", "orange", "dark orange", "red")) +
  geom_text_repel(aes(label = Description), color = "black", size = 3, box.padding = 0.125) +
  facet_wrap( ~ cov_filter, scales = "free") +
  theme_classic() +
  theme(text = element_text(family = "sans"),
        axis.line = element_line(color = 'black', linewidth = 0.125),
        axis.ticks = element_line(colour = "black", linewidth = 0.125),
        axis.title = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6))
#p1
grid.arrange(grobs = lapply(list(p1), set_panel_size, width = unit(4, "cm"), height = unit(4, "cm")))
#dev.off()

#pdf("figures/splicingEventCounts_vs_factors_datatable_covFilt10_only.pdf",useDingbats=F,height=8,width=8)
library(ggplot2)
library(egg)
require(ggrepel)
p1 <- ggplot(for_fig_dt[cov_filter=="n_in_covFilt10"], aes(mean_mapped_counts, event_counts, shape = SE_or_PE, color = as.factor(read_length))) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 12500), name = "Number of splicing changes" , breaks = rep(0:5)*2500) +
  scale_x_continuous(limits = c(0, 100000000), name = "Averarge mapped read counts") +
  scale_color_manual(values = c("light grey", "grey", "orange", "dark orange", "red")) +
  geom_text_repel(aes(label = Description), color = "black", size = 3, box.padding = 0.125) +
  theme_classic() +
  theme(text = element_text(family = "sans"),
        axis.line = element_line(color = 'black', linewidth = 0.125),
        axis.ticks = element_line(colour = "black", linewidth = 0.125),
        axis.title = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6))
grid.arrange(grobs = lapply(list(p1), set_panel_size, width = unit(4, "cm"), height = unit(4, "cm")))
#dev.off()

#-----------------------------------------------------------------# 
# # make figure 6A
# #pdf("figures/splicing_changes_covFilt10_plot.pdf",useDingbats=F,height=5,width=5)
# library(ggplot2)
# library(egg)
# library(ggrastr)
# p1<-ggplot(sig_mats_dt,aes(splice_type,IncLevelDifference,color=log2FoldChange))+
#   rasterize(geom_jitter(width=0.3, size=0.8, shape=16),dpi=2400)+
#   #geom_jitter(width=0.3, size=0.8, shape=16)+
#   #scale_color_gradient2(low="blue",mid="white",high="red",na.value="lightgrey",midpoint=0)+
#   #scale_color_fermenter(palette = "RdBu", n.breaks=5, breaks=c(-1.5,-0.5,0.5,1.5), na.value="lightgrey")+
#   binned_scale(aesthetics = "color",
#                scale_name = "stepsn",
#                palette = function(x) c("#0072b0", "#92c4dd", "#d4d3d3", "#f2a481", "#c92027"),
#                breaks=c(-1.5,-0.5,0.5,1.5),
#                limits = c(-2.5, 2.5),
#                show.limits = TRUE,
#                na.value="white",
#                guide = "colorsteps")+
#   ylim(c(-1,1))+
#   theme_classic()+
#   theme(text=element_text(family="sans"),
#         plot.title=element_text(size=8),
#         axis.line=element_line(color='black',linewidth=0.125),
#         axis.ticks=element_line(colour="black",linewidth = 0.125),
#         axis.title=element_text(size=7),
#         axis.text.y=element_text(size=5),
#         axis.text.x=element_text(size=7),
#         strip.background=element_rect(size=0.25),
#         strip.text=element_text(size=7),
#         panel.border=element_rect(fill=NA,size=NA),
#         panel.spacing=unit(0.25,"line"))
# #p1
# grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(2.75,"cm"),height=unit(3.25,"cm")))
# #dev.off()
# 
# 
# actb_d1_gobp<-gp_covfilt10[dataset=="actb_d1",1:12][1:15]
# actb_d1_gobp[,nLog10:=(-log(p_value,10))]
# 
# write.table(actb_d1_gobp,"./Rds_files/actB_GOBP_splicingChanges_covFilt10_datatable.txt", sep="\t", quote = F, row.names = F)
# gp_covfilt0[dataset=="actb_d1",1:12][1:15]

