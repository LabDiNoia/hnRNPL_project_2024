# compare 3'UTR lengths in unchanged, up- and downregulated genes
library(rtracklayer)
library(data.table)
library(GenomicFeatures)

mouse_gtf<-import("../../../../../common/mm10/gencode.vM23.primary_assembly.annotation.gtf")
mouse_txdb<-makeTxDbFromGRanges(mouse_gtf)
threeUTRs<-threeUTRsByTranscript(mouse_txdb, use.names=T)
width(ranges(threeUTRs))
threeUTR_lengths<-as.data.table( width(ranges(threeUTRs)) )
threeUTR_lengths<-threeUTR_lengths[,.(UTR3_length=sum(value)),by=group_name]
colnames(threeUTR_lengths)[1]<-"transcript_id"

gene_tx_id_lookUp<-unique(as.data.table(mouse_gtf)[!is.na(transcript_id),.(gene_id,transcript_id,gene_type)])
gene_tx_id_lookUp[,.N,gene_type]
threeUTR_lengths<-merge(threeUTR_lengths, gene_tx_id_lookUp, all.x=T, all.y=F)
threeUTR_lengths[,.N,by=gene_type] # no lncRNAs here. a couple of polymorphic_pseudogene, TR_C_gene and IG_C_gene.

# merging with differential gene expression info
actb_norm_dt<-fread("../results/fcounts_deseq/normCounts.txt")
actb_norm_dt[,leg:=NULL]
actb_norm_dt[is.na(padj), leg:="0_unexpr"]
actb_norm_dt[padj<0.1 & log2FoldChange > log(1.5,2), leg:="2_upreg"] #1499
actb_norm_dt[padj<0.1 & log2FoldChange < (-log(1.5,2)), leg:="3_downreg"] #817
actb_norm_dt[is.na(leg), leg:="1_unchanged"]
actb_norm_dt[,gene_id_noVer:=sub("\\.[0-9]*","",gene_id)]
actb_norm_dt$leg<-factor(actb_norm_dt$leg,levels=c("1_unchanged","0_unexpr","2_upreg","3_downreg"))

threeUTR_lengths<-merge(threeUTR_lengths, actb_norm_dt[,.(gene_id,leg)], by="gene_id", all.x=T, all.y=F)
threeUTR_lengths<-threeUTR_lengths[gene_type=="protein_coding"]

threeUTR_lengths[,.N,by=leg]
#           leg     N
# 1: 1_unchanged 29391
# 2:    0_unexpr 18541
# 3:     2_upreg  2780
# 4:   3_downreg  1607

unique(threeUTR_lengths[,.(gene_id,leg)])[,.N,by=leg] 
#           leg     N
# 1: 1_unchanged 10303
# 2:    0_unexpr  9231
# 3:     2_upreg   987
# 4:   3_downreg   605

threeUTR_lengths[, mean(UTR3_length), by=leg]

library(ggplot2)
library(egg)
#ggplot(threeUTR_lengths[leg!="0_unexpr"],aes(x=leg, y=UTR3_length))+geom_violin()+scale_y_continuous(trans = "log10")+geom_boxplot()

p1<-ggplot(threeUTR_lengths[leg!="0_unexpr"],aes(x=leg, y=UTR3_length))+
  geom_violin(aes(fill=leg,color=leg), width=0.8)+
  scale_fill_manual(values=c("grey","red","blue"))+
  scale_color_manual(values=c("grey","red","blue"))+
  scale_y_continuous(trans = "log10")+
  theme_classic()+
  theme(text=element_text(family="sans"),
        axis.line=element_line(color='black',linewidth=0.125),
        axis.ticks=element_line(colour="black",linewidth = 0.125),
        axis.title=element_text(size=7),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6))

#p1

#pdf("figures/3pUTR_lengths_violinPlot.pdf",useDingbats=F,height=5,width=15)
grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(1.5,"cm"),height=unit(2.5,"cm")))
#dev.off()


p1<-ggplot(threeUTR_lengths[leg!="0_unexpr"],aes(x=UTR3_length, color=leg))+
  geom_density()+
  scale_color_manual(values=c("grey","red","blue"))+
  scale_x_continuous(trans = "log10")+
  theme_classic()+
  theme(text=element_text(family="sans"),
        axis.line=element_line(color='black',linewidth=0.125),
        axis.ticks=element_line(colour="black",linewidth = 0.125),
        axis.title=element_text(size=7),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6))
#p1
#pdf("figures/3pUTR_length_densityPlot.pdf",useDingbats=F,height=5,width=15)
grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(2.5,"cm"),height=unit(2.5,"cm")))
#dev.off()

threeUTR_lengths<-merge(threeUTR_lengths,actb_norm_dt[,.(gene_id,gene_name)], all.x=T, all.y=F)
#dir.create("infoTables")
#write.table(threeUTR_lengths,"infoTables/3pUTR_lengths.txt", sep="\t", quote=F, row.names=F)
