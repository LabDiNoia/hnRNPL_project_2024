# some custom figures that did not fit in other places.
# for the sake of organisation, this file will have all such files.

#------------------------------------------------------------------------------#
# lncRNA analysis
# get lncRNA IDs

library(rtracklayer)
mm10_gtf<-import("../../common/mm10/gencode.vM23.primary_assembly.annotation.gtf")
unique(mm10_gtf$gene_type)
mm10_gtf[mm10_gtf$gene_type=="lncRNA" & mm10_gtf$type=="gene"]
mm10_lncRNA_geneIDs<-mm10_gtf[mm10_gtf$gene_type=="lncRNA" & mm10_gtf$type=="gene"]$gene_id
mm10_lncRNA_geneIDs_noVer<-sub("\\.[0-9]*","",mm10_lncRNA_geneIDs)
length(mm10_lncRNA_geneIDs_noVer)
rm(mm10_gtf)

hg38_gtf<-import("../../common/hg38/gencode.v43.primary_assembly.annotation.gtf")
unique(hg38_gtf$gene_type)
hg38_gtf[hg38_gtf$gene_type=="lncRNA" & hg38_gtf$type=="gene"]
hg38_lncRNA_geneIDs<-hg38_gtf[hg38_gtf$gene_type=="lncRNA" & hg38_gtf$type=="gene"]$gene_id
hg38_lncRNA_geneIDs_noVer<-sub("\\.[0-9]*","",hg38_lncRNA_geneIDs)
length(hg38_lncRNA_geneIDs_noVer)
rm(hg38_gtf)

#-------------------------------------#
# calculate percentages of up or downreg genes that are lncRNAs
library(data.table)
allNormcounts_dt<-readRDS("./R_rds_files/allNormcounts.rds")
allNormcounts_dt[,gene_id_noVer:=sub('\\..*','',GeneID),by=.I]

datasets<-allNormcounts_dt[,unique(dataset)]

get_lncProps<-function(dataset){
  cur_dataset<-dataset
  temp_dt<-allNormcounts_dt[dataset==cur_dataset]
  upreg_ids<-temp_dt[padj<0.1 & log2FoldChange>log(1.5,2), gene_id_noVer]
  downreg_ids<-temp_dt[padj<0.1 & log2FoldChange< (-log(1.5,2)), gene_id_noVer]
  
  mouse_or_human<-temp_dt[,unique(mouse_or_human)]
  if(mouse_or_human=="human"){
    lnc_IDs<-hg38_lncRNA_geneIDs_noVer
  } else if (mouse_or_human=="mouse") {
    lnc_IDs<-mm10_lncRNA_geneIDs_noVer
  }
  
  n_up<-length(upreg_ids)
  n_down<-length(downreg_ids)
  n_up_lnc<-sum(lnc_IDs %in% upreg_ids)
  n_down_lnc<-sum(lnc_IDs %in% downreg_ids)
  
  
  
  prop_dt<-data.table(DGE_direction = "Upreg",
                      n_deg = n_up,
                      n_deg_lnc = n_up_lnc,
                      lncRNA_perc = n_up_lnc/n_up*100)
  prop_dt<-rbind(prop_dt,
                 data.table(DGE_direction="Downreg",
                            n_deg = n_down,
                            n_deg_lnc = n_down_lnc,
                            lncRNA_perc = n_down_lnc/n_down*100))
  prop_dt[,dataset:=cur_dataset]
  prop_dt[,mouse_or_human:=mouse_or_human]
  return(prop_dt)
}

prop_dt<-lapply(datasets, get_lncProps)
prop_dt<-rbindlist(prop_dt)
prop_dt

write.table(prop_dt, "R_rds_files/lncRNA_DEG_summary.txt", row.names=F, sep="\t", quote=F)
require(data.table)
prop_dt<-fread("R_rds_files/lncRNA_DEG_summary.txt")

#-------------------------------------#
# make bar graphs to show the difference
# color values
colorsForDatasets<-c("#d7403c","#775fa9","#92ab3c","#29b34a","#2db99c","#0bb9e2","#6f94cc","#b37bb5","#d29329")
names(colorsForDatasets)<-c("actb_d1","thymocytes","fetalLiv","HEK293T","hepg2","k562","kcytes","lncap","bjFibro")

prop_dt$DGE_direction<-factor(prop_dt$DGE_direction, levels=c("Upreg","Downreg"))
prop_dt$mouse_or_human<-factor(prop_dt$mouse_or_human, levels=c("mouse","human"))
prop_dt$dataset<-factor(prop_dt$dataset,
                        levels=prop_dt[DGE_direction=="Upreg"][order(mouse_or_human,lncRNA_perc), dataset])

library(ggplot2)
library(egg)
p1<-ggplot(prop_dt,aes(x=DGE_direction, y=lncRNA_perc, fill=dataset))+
  geom_bar(position="dodge", stat="identity", width=0.8)+
  scale_fill_manual(values=colorsForDatasets)+
  scale_y_continuous(expand=expansion(mult = c(0, 0.1)))+
  theme_classic()+
  theme(text=element_text(family="sans"),
        plot.title=element_text(size=8),
        axis.line=element_line(color='black',linewidth=0.125),
        axis.ticks=element_line(colour="black",linewidth = 0.125),
        axis.title=element_text(size=7),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=7),
        strip.background=element_rect(size=0.25),
        strip.text=element_text(size=7),
        panel.border=element_rect(fill=NA,size=NA),
        panel.spacing=unit(0.25,"line"))

#pdf("figures/lncRNA_percentages.pdf",useDingbats=F,height=5,width=5)
grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(4,"cm"),height=unit(2.5,"cm")))
#dev.off()

#------------------------------------------------------------------------------#
# calculate and plot hnRNPL:hnRNPLL expression ratios

library(data.table)
allNormcounts_dt<-readRDS("./R_rds_files/allNormcounts.rds")
allNormcounts_dt[,gene_id_noVer:=sub('\\..*','',GeneID),by=.I]

datasets<-allNormcounts_dt[,unique(dataset)]

dataset<-datasets[1]

get_L_LL_ratio<-function(dataset){
  cur_dataset<-dataset
  temp_dt<-allNormcounts_dt[dataset==cur_dataset]
  L_val<-temp_dt[grep("\\<Hnrnpl\\>",gene_name,ignore.case=T), CTL_mean]
  LL_val<-temp_dt[grep("\\<Hnrnpll\\>",gene_name,ignore.case=T), CTL_mean]
  L_LL_ratio<-L_val/LL_val
  return(L_LL_ratio)
}

L_LL_ratios<-sapply(datasets,get_L_LL_ratio)
L_LL_ratios_dt<-data.table(dataset=names(L_LL_ratios), L_LL_ratio=L_LL_ratios)
L_LL_ratios_dt<-merge(L_LL_ratios_dt, unique(allNormcounts_dt[,.(dataset,mouse_or_human)]))

# color values
colorsForDatasets<-c("#d7403c","#775fa9","#92ab3c","#29b34a","#2db99c","#0bb9e2","#6f94cc","#b37bb5","#d29329")
names(colorsForDatasets)<-c("actb_d1","thymocytes","fetalLiv","HEK293T","hepg2","k562","kcytes","lncap","bjFibro")

require(ggplot2)
require(egg)
require(ggbreak)
p1<-ggplot(L_LL_ratios_dt, aes(y=L_LL_ratio, color=dataset, shape=mouse_or_human))+
  geom_point(aes(x="test"))+
  scale_shape_manual(values=c(20,18))+
  scale_color_manual(values=colorsForDatasets)+
  ylim(0,600)+
  scale_y_cut(breaks=10, which=1, scales=0.25)+
  theme_classic()+
  theme(text=element_text(family="sans"),
        plot.title=element_text(size=8),
        axis.line=element_line(color='black',linewidth=0.125),
        axis.ticks=element_line(colour="black",linewidth = 0.125),
        axis.title=element_text(size=7),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=7),
        strip.background=element_rect(size=0.25),
        strip.text=element_text(size=7),
        panel.border=element_rect(fill=NA,size=NA),
        panel.spacing=unit(0.25,"line"))


pdf("figures/hnRNPL_LL_ratios.pdf",useDingbats=F,height=2.5,width=2.5)
p1
dev.off()
#------------------------------------------------------------------------------#
