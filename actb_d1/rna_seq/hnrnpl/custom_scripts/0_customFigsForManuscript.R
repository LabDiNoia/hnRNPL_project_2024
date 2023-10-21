# some custom figures that did not fit in other places.
# for the sake of organisation, this file will have all such files.

#------------------------------------------------------------------------------#
#check numbers of up- and down-regulated genes
#make volcano plot with names written for some of the hits
library(data.table)
actb_norm_dt<-fread("../results/fcounts_deseq/normCounts.txt")

# checking the number of genes with thresholds
actb_norm_dt[padj<0.2 & log2FoldChange > 0] #2797
actb_norm_dt[padj<0.2 & log2FoldChange < 0] #2344

actb_norm_dt[padj<0.1 & log2FoldChange > 0] #2072
actb_norm_dt[padj<0.1 & log2FoldChange < 0] #1623

actb_norm_dt[padj<0.2 & log2FoldChange > log(1.5,2)] #1837
actb_norm_dt[padj<0.2 & log2FoldChange < (-log(1.5,2))] #982

# let's stick to this one here
actb_norm_dt[padj<0.1 & log2FoldChange > log(1.5,2)] #1499
actb_norm_dt[padj<0.1 & log2FoldChange < (-log(1.5,2))] #817

actb_norm_dt[padj<0.05 & log2FoldChange > log(1.5,2)] #1280
actb_norm_dt[padj<0.05 & log2FoldChange < (-log(1.5,2))] #681


actb_norm_dt[,leg:=NULL]
actb_norm_dt[is.na(padj), leg:="0_unexpr"]
actb_norm_dt[padj<0.1 & log2FoldChange > log(1.5,2), leg:="2_upreg"] #1499
actb_norm_dt[padj<0.1 & log2FoldChange < (-log(1.5,2)), leg:="3_downreg"] #817
actb_norm_dt[is.na(leg), leg:="1_unchanged"]

actb_norm_dt[,gene_id_noVer:=sub("\\.[0-9]*","",gene_id)]
actb_norm_dt$leg<-factor(actb_norm_dt$leg,levels=c("1_unchanged","0_unexpr","2_upreg","3_downreg"))

#pdf("figures/actb_volcano_plot.pdf",useDingbats=F,height=3,width=5)
library(ggplot2)
library(egg)
library(ggrastr)
library(ggrepel)
p1<-ggplot(data=actb_norm_dt[leg!="0_unexpr"], aes(x=log2FoldChange, y=-log(padj,10), color=leg))+
  #rasterize(geom_point(size=0.8,shape=16),dpi=2400)+
  geom_point(size=0.8,shape=16)+
  scale_y_continuous("-log10(adj p-value)")+
  scale_x_continuous("log2 fold change")+
  scale_color_manual(values=c("grey","red","blue"))+
  geom_text_repel(data=subset(actb_norm_dt,-log(padj,10)>50),
                  aes(label=gene_name),
                  size=2, box.padding = 0.05)+
  theme_classic()+
  theme(text=element_text(family="sans"),
        axis.line=element_line(color='black',linewidth=0.125),
        axis.ticks=element_line(colour="black",linewidth = 0.125),
        axis.title=element_text(size=7),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6))
#p1
grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(4,"cm"),height=unit(3,"cm")))
#dev.off()

#------------------------------------------------------------------------------#
# plot differentially spliced genes and overlap with DEGs

# get file names
splicing_files<-list.files("../results/rMATS/", pattern="*.MATS.JCEC.txt", full.names = T)
names(splicing_files)<-sub(".MATS.JCEC.txt", "", basename(splicing_files))

# keep only info relevant for now
# filter for FDR<0.1 & IncLevelDifference>abs(0.1), keep only gene ID/name, FDR and IncLevelDifference 
actb_allSplice_dt<-lapply(splicing_files,
                          function(file_path){
                            require(data.table)
                            temp_dt<-fread(file_path)
                            temp_dt[,splicing_type:=sub(".MATS.JCEC.txt", "", basename(file_path))]
                            temp_dt<-temp_dt[FDR<0.1 & abs(IncLevelDifference)>0.1]
                            temp_dt[,.(splicing_type,GeneID,geneSymbol,FDR,IncLevelDifference)]
                          })
actb_allSplice_dt<-rbindlist(actb_allSplice_dt)
actb_allSplice_dt # 5801 rows
actb_allSplice_dt[,.N,by=splicing_type]
#   splicing_type    N
# 1:          A3SS  657
# 2:          A5SS  543
# 3:           MXE  796
# 4:            RI  740
# 5:            SE 3065
actb_allSplice_dt[IncLevelDifference<0,.N,by=splicing_type] #number of events
actb_allSplice_dt[IncLevelDifference>0,.N,by=splicing_type] #number of events
actb_allSplice_dt[,length(unique(GeneID)),by=splicing_type] #number of genes
actb_allSplice_dt[,gene_id_noVer:=sub("\\.[0-9]*","",GeneID)]

actb_allSplice_dt<-merge(actb_allSplice_dt, actb_norm_dt[padj<0.1,.(gene_id_noVer,expr_log2FC=log2FoldChange)], all.x=T)
actb_allSplice_dt<-rbind(actb_allSplice_dt[is.na(expr_log2FC)], actb_allSplice_dt[!is.na(expr_log2FC)])

#pdf("figures/splicing_changes_plot.pdf",useDingbats=F,height=5,width=5)
library(ggplot2)
library(egg)
library(ggrastr)
p1<-ggplot(actb_allSplice_dt,aes(splicing_type,IncLevelDifference,color=expr_log2FC))+
  #rasterize(geom_jitter(width=0.3, size=0.8, shape=16),dpi=2400)+
  geom_jitter(width=0.3, size=0.8, shape=16)+
  #scale_color_gradient2(low="blue",mid="white",high="red",na.value="lightgrey",midpoint=0)+
  scale_color_fermenter(palette = "RdBu", n.breaks=5, breaks=c(-1.5,-0.5,0.5,1.5), na.value="lightgrey")+
  ylim(c(-1,1))+
  theme_classic()+
  theme(text=element_text(family="sans"),
        plot.title=element_text(size=8),
        axis.line=element_line(color='black',linewidth=0.125),
        axis.ticks=element_line(colour="black",linewidth = 0.125),
        axis.title=element_text(size=7),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=7),
        strip.background=element_rect(size=0.25),
        strip.text=element_text(size=7),
        panel.border=element_rect(fill=NA,size=NA),
        panel.spacing=unit(0.25,"line"))
#p1
grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(2.75,"cm"),height=unit(3.25,"cm")))
#dev.off()

#------------------------------------------------------------------------------#
#calculating overlaps between splicing changes and DEGs

length(unique(actb_allSplice_dt[,gene_id_noVer])) #2836
length(actb_norm_dt[leg %in% c("2_upreg","3_downreg"), gene_id_noVer]) #2316
intersect(actb_allSplice_dt[,gene_id_noVer], actb_norm_dt[leg %in% c("2_upreg","3_downreg"), gene_id_noVer]) #367
2836-367
2316-367

dir.create("./genelists")

#print intersecting genes
sink("./genelists/367_intersecting_splicing_DGE_genes.txt")
cat(intersect(actb_allSplice_dt[,gene_id_noVer], actb_norm_dt[leg %in% c("2_upreg","3_downreg"), gene_id_noVer]),sep="\n")
sink()

#print splicing changes only genes
sink("./genelists/2469_genes_with_splicingChanges_only.txt")
cat(setdiff(actb_allSplice_dt[,gene_id_noVer], actb_norm_dt[leg %in% c("2_upreg","3_downreg"), gene_id_noVer]),sep="\n")
sink()

#print upregulated genes
sink("./genelists/1499_upregulated_genes.txt")
cat(unique(actb_norm_dt[leg=="2_upreg", gene_id_noVer]),sep="\n")
sink()

#print downregulated genes
sink("./genelists/817_downregulated_genes.txt")
cat(unique(actb_norm_dt[leg=="3_downreg", gene_id_noVer]),sep="\n")
sink()

#print all genes
sink("./genelists/55421_all_genes.txt")
cat(unique(actb_norm_dt[, gene_id_noVer]),sep="\n")
sink()

#print genes with changes in 'SE' (skipped exon)
sink("./genelists/1747_genes_with_SE_changes.txt")
cat(actb_allSplice_dt[splicing_type=="SE",unique(gene_id_noVer)],sep="\n")
sink()

#------------------------------------------------------------------------------#
# plotting gene expression changes

actb_norm_dt[gene_name=="Hnrnpl"]

temp_dt<-actb_norm_dt[gene_name=="Hnrnpl"]
# set 1
genes_to_plot<-c("Bcl6","Aicda",
                 "Myc","E2f1","E2f2","E2f3","E2f4","Ccnd2","Ccnd3",
                 "Trp53","Cdkn1a","Trp53cor1","Bcl2l11",
                 "Bbc3","Pmaip1","Bax","Bad",
                 "Bcl2","Mcl1","Bcl2l1",
                 "Sirt1","Nsd2","Setd2",
                 "Ighm","Ighg1", "Foxo3")
temp_dt<-actb_norm_dt[gene_name %in% genes_to_plot]

forPlot_dt<-melt(temp_dt, id.vars=c("gene_name","gene_id"),
                 measure.vars=c("WT1","WT2","WT3","KO1","KO2","KO3"),
                 variable.name = "sample", value.name = "expr")
forPlot_dt[,genotype:=sub("[0-9]","",sample)]
forPlot_dt$genotype<-factor(forPlot_dt$genotype, levels=c("WT","KO"))
forPlot_dt$gene_name<-factor(forPlot_dt$gene_name, levels=genes_to_plot)

max_y_val<-1.1*forPlot_dt[,max(expr)]

p1<-ggplot(forPlot_dt,aes(genotype, expr, fill=genotype))+
  stat_summary(geom = "bar", fun = mean, width=1, position = "dodge")+
  scale_fill_manual(values=c("#B9B8B8","#477FB7"))+
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_y_continuous(expand=expansion(mult = c(0, 0.05)))+
  geom_jitter(width=0.2, size=0.8, shape=16)+
  facet_wrap(~gene_name, scales = "free_y", nrow=1)+
  theme_classic()+
  theme(text=element_text(family="sans"),
        plot.title=element_text(size=8),
        axis.line=element_line(color='black',linewidth=0.125),
        axis.ticks=element_line(colour="black",linewidth = 0.125),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=7),
        axis.text.y=element_text(size=5),
        axis.text.x=element_blank(),
        strip.background=element_rect(linewidth=0.25),
        strip.text=element_text(size=7),
        panel.border=element_rect(fill=NA,size=NA),
        panel.spacing=unit(0.25,"line"))

#pdf("figures/rna_seq_indv_genes_expr.pdf",useDingbats=F,height=5,width=15)
grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(0.45,"cm"),height=unit(2.5,"cm")))
#dev.off()


# set 2
genes_to_plot<-c("Mfn1","Mfn2","Opa1","Dnm1l",
                 "Tfam","Ppargc1a","Ppargc1b",
                 "Cat","Sod1","Sod2","Gpx1",
                 "Ndufs1","Ndufs5","Uqcrc1","Cycs",
                 "Cox5b","Cox6c","Atp5e","Atp5b")
temp_dt<-actb_norm_dt[gene_name %in% genes_to_plot]

forPlot_dt<-melt(temp_dt, id.vars=c("gene_name","gene_id"),
                 measure.vars=c("WT1","WT2","WT3","KO1","KO2","KO3"),
                 variable.name = "sample", value.name = "expr")
forPlot_dt[,genotype:=sub("[0-9]","",sample)]
forPlot_dt$genotype<-factor(forPlot_dt$genotype, levels=c("WT","KO"))
forPlot_dt$gene_name<-factor(forPlot_dt$gene_name, levels=genes_to_plot)

max_y_val<-1.1*forPlot_dt[,max(expr)]

require(ggplot2)
require(egg)
p1<-ggplot(forPlot_dt,aes(genotype, expr, fill=genotype))+
  stat_summary(geom = "bar", fun = mean, width=1, position = "dodge")+
  scale_fill_manual(values=c("#B9B8B8","#477FB7"))+
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_y_continuous(expand=expansion(mult = c(0, 0.05)))+
  geom_jitter(width=0.2, size=0.8, shape=16)+
  facet_wrap(~gene_name, scales = "free_y", nrow=1)+
  theme_classic()+
  theme(text=element_text(family="sans"),
        plot.title=element_text(size=8),
        axis.line=element_line(color='black',linewidth=0.125),
        axis.ticks=element_line(colour="black",linewidth = 0.125),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=7),
        axis.text.y=element_text(size=5),
        axis.text.x=element_blank(),
        strip.background=element_rect(linewidth=0.25),
        strip.text=element_text(size=7),
        panel.border=element_rect(fill=NA,size=NA),
        panel.spacing=unit(0.25,"line"))

#pdf("figures/rna_seq_indv_genes_expr_mitoch_set.pdf",useDingbats=F,height=5,width=15)
grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(0.45,"cm"),height=unit(2.5,"cm")))
#dev.off()



