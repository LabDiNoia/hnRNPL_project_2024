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
library(data.table)
commonID_forHomologs_dt<-readRDS("./R_rds_files/commonID_forHomologs_dt.rds")

# step2: get gene_ids of differently spliced from the different datasets; translate it to commonID; do upset plot
#folders<-c("../HEK293T/","../actb_d1/","../bjFibro/","../c2c12/","../fetalLiv/","../hepg2/","../k562/","../kcytes/","../lncap/","../thymocytes/")
# removed c2c12 because only one rep: sketchy stats
folders<-c("../HEK293T/","../bjFibro/","../hepg2/","../k562/","../kcytes/","../lncap/","../thymocytes/","../fetalLiv/","../actb_d1/")
datasets<-basename(folders)
splicing_type<-c("SE")

####################################################################################3
# goal: upset plot with splicing events not genes

# 1. get all splicing events with FDR<0.2 in one table
get_spl_dt<-function(dataset,alt_splice_tpye){
        # path to file
        dataset_folder<-paste0("../",dataset)
        file_path<-paste0(dataset_folder,"/rna_seq/hnrnpl/results/rMATS/",alt_splice_tpye,".MATS.JCEC.txt")
        # get file
        rmats_dt<-fread(file_path)
        # filter
        rmats_dt<-rmats_dt[FDR<0.2]
	rmats_dt[,dataset:=dataset]
	rmats_dt[,ID:=NULL]
        return(rmats_dt)
}
allSig_SE<-sapply(datasets,get_spl_dt,"SE",simplify=F)
allSig_SE<-rbindlist(allSig_SE)

# get commonID for each gene
# melt commonID_forHomologs_dt to enable combining
melted_commonIDs_dt<-melt(commonID_forHomologs_dt, id.vars=c("commonID"), measure.vars=c("mouse_id","human_id"), variable.name="mouse_or_human_id", value.name="GeneID", na.rm=T)
melted_commonIDs_dt[,mouse_or_human:=sub('_id','',mouse_or_human_id),by=.I]
melted_commonIDs_dt[,mouse_or_human_id:=NULL]

allSig_SE<-merge(melted_commonIDs_dt,allSig_SE,all.y=T,all.x=F)
#---------------------------------------------------------------------#
# 2. find homologous splicing events and assign a commonSpliceID to such splicing events
# libraries to get sequences
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)

# keep only columns till "downstreamEE". keep only unique events (regardless of dataset)
uniqSig_SE<-unique(allSig_SE[,1:match("downstreamEE",colnames(allSig_SE))])
# assign spliceID for each splice event
uniqSig_SE[,spliceID:=paste(commonID,1:.N,sep="_"),by=commonID]
# is the homolog also differentially spliced? also remove non-standardChromosomes
uniqSig_SE[,homolog_diffSpliced:=length(unique(mouse_or_human))==2,by=commonID]
# also remove non-standardChromosomes, because getting sequences is a pain
stdChr<-unique(c(standardChromosomes(Mmusculus),standardChromosomes(Hsapiens)))
uniqSig_SE[commonID %in% uniqSig_SE[!chr %in% stdChr, commonID], homolog_diffSpliced:=FALSE]
# number of genes differentially spliced in both human and mouse in any dataset is 3966
length(unique(uniqSig_SE[homolog_diffSpliced==T,commonID])) # 3966

# workhorse function for calling pairwiseAlignment for the 3 exons involved in each event
get_alnInfo<-function(DTforGR,whichExon){
	# define exonStartEndColDefText variable. to define start and end in GR
	# depending on the strand, upstream and downstream exons are different
	if(whichExon=="current"){
		exonStartEndColDefText="start=exonStart_0base,end=exonEnd"
	} else if(whichExon=="upstream"){
		DTforGR[strand=="+",upStart:=upstreamES,by=.I]
		DTforGR[strand=="+",upEnd:=upstreamEE,by=.I]
		DTforGR[strand=="-",upStart:=downstreamES,by=.I]
                DTforGR[strand=="-",upEnd:=downstreamEE,by=.I]
		exonStartEndColDefText="start=upStart,end=upEnd"
	} else if(whichExon=="downstream"){
		DTforGR[strand=="+",downStart:=downstreamES,by=.I]
                DTforGR[strand=="+",downEnd:=downstreamEE,by=.I]
                DTforGR[strand=="-",downStart:=upstreamES,by=.I]
                DTforGR[strand=="-",downEnd:=upstreamEE,by=.I]
                exonStartEndColDefText="start=downStart,end=downEnd"
	}

	human_gr_expr<-paste0('human_gr<-
		makeGRangesFromDataFrame(DTforGR[mouse_or_human=="human", .(chr,strand,',
		exonStartEndColDefText,
		',spliceID)], keep.extra.columns=T)')
	eval(parse(text=human_gr_expr))

	mouse_gr_expr<-paste0('mouse_gr<-
		makeGRangesFromDataFrame(DTforGR[mouse_or_human=="mouse", .(chr,strand,',
		exonStartEndColDefText,
		',spliceID)], keep.extra.columns=T)')
	eval(parse(text=mouse_gr_expr))

	if(whichExon=="current"){
		assign("human_gr", human_gr, envir=.GlobalEnv)
		assign("mouse_gr", mouse_gr, envir=.GlobalEnv)
	}

	# get sequences
	human_seqs<-Views(Hsapiens,human_gr)
	mouse_seqs<-Views(Mmusculus,mouse_gr)

	# getting pairwiseAlingment score matrix. maybe think about changing gapExtension
	exonAlnInfo_mat<-sapply(human_seqs,function(x){
		pairwiseAlignment(pattern=mouse_seqs,subject=x,type="global")})
	return(exonAlnInfo_mat)
}

extractAlnDetails<-function(alnInfo_matrixList,functionToApply){
		sapply(names(alnInfo_matrixList),function(x){
			matList<-matrix(sapply(alnInfo_matrixList[[x]],functionToApply), 
					ncol=length(human_gr$spliceID))
			colnames(matList)<-human_gr$spliceID
			rownames(matList)<-mouse_gr$spliceID
			return(matList)}, simplify=F)
}

commonIDs_homolog_diffSpliced<-unique(uniqSig_SE[homolog_diffSpliced==T,commonID])
nhomologs<-length(commonIDs_homolog_diffSpliced)

get_commonSpliceID_lookup_dt<-function(cur_ID){
	spilce_homolog_status<-FALSE
        start_time<-Sys.time()
        cur_SE_dt<-uniqSig_SE[commonID==cur_ID]
        alnInfo_matList<-sapply(c("current","upstream","downstream"),function(x){
                                get_alnInfo(cur_SE_dt,x)}, simplify=F)

        # get percentage identity of the sequences
        alnPID<-extractAlnDetails(alnInfo_matList,pid)
        # reduce to get truth table for PID>70% for all three exons
        aln_over70p_TT<-Reduce(f=function(x,y){x&y}, x=sapply(alnPID,function(x)x>70,simplify=F))

        #only proceed to commonSpliceID_lookup_dt, if there are any homologous events
        if(sum(aln_over70p_TT)>0){
                # get score and keep scores for homologous splicing events only
                alnScore<-extractAlnDetails(alnInfo_matList,score)
                alnScore_homMat<-Reduce(f=function(x,y){x+y}, x=alnScore) * aln_over70p_TT

                # get and write commonSpliceIDs (homologous splice event) to spliceID.
                # the logic is to:
                # 1. make 0 = -Inf. then find index of max score
                # 2. get row- and colnames of the maxvalue (i.e., mouse and human spliceIDs)
                # 3. assign human spliceID as the commonSpliceID to mouse spliceID
                # 4. remove this row and column from matrix
                # 5. repeat steps 1-4 till only -Inf is left
                # 6. assign spliceID as commonSpliceID for ones without homologs
                # 7. rbind spliceID,commonSpliceID table
                alnScore_mat_forExtr<-alnScore_homMat
                alnScore_mat_forExtr[alnScore_mat_forExtr==0]=-Inf
                while( sum(alnScore_mat_forExtr!=-Inf) > 0){
                        max_ind<-which(alnScore_mat_forExtr==max(alnScore_mat_forExtr),arr.ind=T)
                        human_spliceID<-colnames(alnScore_mat_forExtr)[max_ind[1,2]]
                        mouse_spliceID<-rownames(alnScore_mat_forExtr)[max_ind[1,1]]
                        cur_SE_dt[spliceID==mouse_spliceID,commonSpliceID:=human_spliceID]
                        alnScore_mat_forExtr<-alnScore_mat_forExtr[-max_ind[1,1],-max_ind[1,2],drop=FALSE]
			spilce_homolog_status<-TRUE
                }
        }
        counter<-match(cur_ID,commonIDs_homolog_diffSpliced)
        cat("Done with",counter ,"of",nhomologs ,":", cur_ID,"(",cur_SE_dt[1,geneSymbol],"). ",sep=" ")
        cat("Took",Sys.time()-start_time,"seconds","\n")

	if(spilce_homolog_status){
	cur_SE_dt[is.na(commonSpliceID),commonSpliceID:=spliceID,by=.I]
	return(cur_SE_dt[,.(spliceID,commonSpliceID)])
	}
}

library(parallel)
thread_count<-as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if(is.na(thread_count)){
	thread_count<-min(detectCores(),16)
}
commonSpliceID_lookup_dt<-mclapply(commonIDs_homolog_diffSpliced,get_commonSpliceID_lookup_dt,mc.cores=thread_count)
commonSpliceID_lookup_dt<-rbindlist(commonSpliceID_lookup_dt)

saveRDS(commonSpliceID_lookup_dt,"./R_rds_files/commonSpliceID_lookup_dt_FDR_lt_pt2.rds")
commonSpliceID_lookup_dt<-readRDS("./R_rds_files/commonSpliceID_lookup_dt_FDR_lt_pt2.rds")

# merge commonSpliceID_lookup_dt and uniqSig_SE
temp<-copy(uniqSig_SE)
uniqSig_SE<-merge(uniqSig_SE,commonSpliceID_lookup_dt,by="spliceID",all.x=T)
uniqSig_SE[is.na(commonSpliceID),commonSpliceID:=spliceID,by=.I]
nrow(uniqSig_SE) # 60365
length(unique(uniqSig_SE[,commonSpliceID])) # 59155

# merge uniqSig_SE with commonSpliceIDs and allSig_SE
merged_allSig_SE<-merge(allSig_SE, uniqSig_SE, by=colnames(allSig_SE)[ 1:match("downstreamEE",colnames(allSig_SE)) ])
saveRDS(merged_allSig_SE,"./R_rds_files/merged_allSig_SE.rds")
require(data.table)
merged_allSig_SE<-readRDS("./R_rds_files/merged_allSig_SE.rds")

#make tables for the paper
merged_allSig_SE[,alt_spliced_in_n_datasets:=.N,by=commonSpliceID]
merged_allSig_SE[mouse_or_human=="mouse", status:=if(.N==3){"Alt spliced in all mouse datasets"}else{""}, by=commonSpliceID]
merged_allSig_SE[mouse_or_human=="human", status:=if(.N==6){"Alt spliced in all human datasets"}else{""}, by=commonSpliceID]

merged_allSig_SE[status=="Alt spliced in all mouse datasets" & alt_spliced_in_n_datasets==3, status:="Alt spliced in all mouse datasets but not in any human dataset"]
merged_allSig_SE[status=="Alt spliced in all human datasets" & alt_spliced_in_n_datasets==6, status:="Alt spliced in all human datasets but not in any mouse dataset"]


merged_allSig_SE[status=="Alt spliced in all mouse datasets but not in any human dataset"]
merged_allSig_SE[status=="Alt spliced in all human datasets but not in any mouse dataset"]

merged_allSig_SE<-merged_allSig_SE[order(commonSpliceID,mouse_or_human,dataset)]
merged_allSig_SE<-merged_allSig_SE[order(alt_spliced_in_n_datasets, decreasing=T)]
write.table(merged_allSig_SE, "./R_rds_files/allSig_SE_DSE.txt", sep="\t", quote=T, row.names=F)

merged_allSig_SE[alt_spliced_in_n_datasets>=8 | status!=""]


# hits in at least 8 / 9 datasets
merged_allSig_SE[,commonSpliceID,by=dataset][,.N,by=commonSpliceID][N>=8][order(N,decreasing=T)]

# make list for upset plot
SE_list<-sapply(datasets,function(x){merged_allSig_SE[dataset==x,unique(commonSpliceID)]})


####################################################################################3
library(UpSetR)
# extract relevant info for upset plot
infoForUpset<-unique(merged_allSig_SE[,.(mouse_or_human,spl_count=length(unique(commonSpliceID))), by=dataset])
# order datasets in the order for display
infoForUpset<-infoForUpset[order(mouse_or_human,spl_count,decreasing = T)]

#change SE_list order to 
intersect_list<-c(unlist(lapply(9:8,function(x){combn(datasets,x,simplify = F)}),recursive = F),
     list(infoForUpset[mouse_or_human=="human",dataset]),
     list(infoForUpset[mouse_or_human=="mouse",dataset]))

pdf(file="./figures/splice_SE_upset_FDR_lt_pt2.pdf",useDingbats=F,width=3,height=3)
upset(fromList(SE_list), intersections = rev(intersect_list),
      sets=infoForUpset[,dataset], keep.order=T,
      text.scale=1, mb.ratio=c(0.55,0.45),
      point.size=1.5, line.size = 0.3)
dev.off()

#------------------------------#
# getting and plotting exon inclusion levels
get_incLevel<-function(cur_commonSpliceID){
  cur_dt<-merged_allSig_SE[commonSpliceID==cur_commonSpliceID]
  cur_dt[,common_name:=paste(unique(geneSymbol),collapse="|")]
  cur_dt[,CTL:=lapply(strsplit(IncLevel1,","), function(x){mean(as.numeric(x))}),by=dataset]
  cur_dt[,hnRNPL_dep:=lapply(strsplit(IncLevel2,","), function(x){mean(as.numeric(x))}),by=dataset]
  cur_dt<-melt(cur_dt, id.vars=c("dataset","mouse_or_human","common_name"), measure.vars=c("CTL","hnRNPL_dep"), variable.name="genotype", value.name="IncLevel")
  return(cur_dt)
}

# for hits in at least 9 datasets
merged_allSig_SE[,commonSpliceID,by=dataset][,.N,by=commonSpliceID][N==9]

get_incLevel("36496_3") #GPBP1
get_incLevel("37776_1") #SIRT1
get_incLevel("42830_1") #ZRANB2

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

plot_incLevel("36496_3", saveToPDF=T) #GPBP1
plot_incLevel("37776_1", saveToPDF=T) #SIRT1
plot_incLevel("42830_1", saveToPDF=T) #ZRANB2

merged_allSig_SE[,commonSpliceID,by=dataset][,.N,by=commonSpliceID][N==8]
plot_incLevel("38235_3") #HM13
plot_incLevel("39483_2", saveToPDF=T) #NSD2 - strong
plot_incLevel("39483_3", saveToPDF=T, geneNameForPDF="NSD2_2") #NSD2 - strong
plot_incLevel("41099_1", saveToPDF=T) #EPC1
plot_incLevel("41936_3") #hnRNPR
plot_incLevel("42612_4") #RBM39 - weak
plot_incLevel("42830_2") #ZRANB2
plot_incLevel("43591_3") #SMPD4
plot_incLevel("44005_2") #SEC31A
plot_incLevel("44005_3") #SEC31A
plot_incLevel("45369_2", saveToPDF=T) #KDM6A
plot_incLevel("49981_3", saveToPDF=T) #BPTF - really strong in humans
plot_incLevel("53387_1", saveToPDF=T) #MORF4L1 - very weak

in8of9<-merged_allSig_SE[,commonSpliceID,by=dataset][,.N,by=commonSpliceID][N==8,commonSpliceID]
in8of9_names<-sapply(in8of9, function(x){merged_allSig_SE[commonSpliceID%in%x & mouse_or_human=="human",unique(geneSymbol)]})
data.frame(in8of9_names)
# 38235_3         HM13      missing in k562
# 39483_2         NSD2      missing in bjFibro
# 39483_3         NSD2      missing in bjFibro
# 41099_1         EPC1      missing in bjFibro
# 41936_3       HNRNPR      missing in lncap
# 42612_4        RBM39      missing in thymocytes
# 42830_2       ZRANB2      missing in thymocytes
# 43591_3        SMPD4      missing in thymocytes
# 44005_2       SEC31A      missing in lncap
# 44005_3       SEC31A      missing in lncap
# 45369_2        KDM6A      missing in bjFibro
# 49981_3         BPTF      missing in thymocytes
# 53387_1      MORF4L1      missing in thymocytes

# to find which is missing in which dataset and update the table above
sapply(in8of9, function(x){
  temp<-merged_allSig_SE[commonSpliceID%in%x,c(sort(dataset))]
  setdiff(datasets,temp)})

#------------------------------------------------------------#
# getting percentile expression values

# get datasets
library(data.table)
allNormcounts_dt<-readRDS("./R_rds_files/allNormcounts.rds")

get_expr_perc<-function(dataset, commonID){
  cur_dataset=dataset
  cur_commonID=commonID
  # keeping only "expressed" counts
  cur_dt<-allNormcounts_dt[dataset==cur_dataset][!is.na(padj)]
  # ecdf outputs a function. again same as above, removing 0s
  ecdf_func_CTL<-ecdf(cur_dt[CTL_mean!=0,CTL_mean])
  ecdf_func_Ldep<-ecdf(cur_dt[hnRNPL_dep_mean!=0,hnRNPL_dep_mean])
  # use the outputted function to get the percentil value
  CTL_perc<-ecdf_func_CTL(cur_dt[commonID==cur_commonID,CTL_mean])
  Ldep_perc<-ecdf_func_Ldep(cur_dt[commonID==cur_commonID,hnRNPL_dep_mean])
  # make data.table formatted long for easy plotting in ggplot
  cur_mouse_or_human<-cur_dt[,unique(mouse_or_human)]
  out_dt<-rbind(data.table(genotype="CTL", percentile=CTL_perc),
                    data.table(genotype="hnRNPL_dep", percentile=Ldep_perc))
  out_dt[,dataset:=cur_dataset]
  out_dt[,mouse_or_human:=cur_mouse_or_human]
  return(out_dt)
}

get_expr_perc("actb_d1","37776")
Sirt1_percExpr_vals<-rbindlist(sapply(datasets,get_expr_perc,"37776",simplify=F))
hnRNPL_percExpr_vals<-rbindlist(sapply(datasets,get_expr_perc,"38727",simplify=F))

# plot CTL vs hnRNPL_dep perc expr plots
#pdf("figures/Sirt1_perc_expr.pdf",useDingbats=F,height=5,width=5)
require(ggplot2)
require(egg)
p1<-ggplot(data=hnRNPL_percExpr_vals, aes(x=genotype, y=percentile*100, group=dataset, color=dataset, shape=mouse_or_human))+
  geom_point()+
  geom_line()+
  scale_y_continuous("Expression levels (percentile)")+
  scale_shape_manual(values=c(20,18))+
  scale_color_manual(values=colorsForDatasets)+
  theme_classic()+
  theme(text=element_text(family="sans"),
        axis.line=element_line(color='black',size=0.25),
        axis.ticks=element_line(colour="black",size = 0.25),
        axis.title=element_text(size=7),
        axis.text.y=element_text(angle=90,size=5,hjust=0.5),
        axis.text.x=element_text(size=5),
        strip.background=element_rect(size=0.25),
        strip.text=element_text(size=7),
        panel.border=element_rect(fill=NA,size=0.25),
        panel.spacing=unit(0,"line"))
grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(1.5,"cm"),height=unit(2.5,"cm")))
#dev.off()

#------------------------------------------------------------#
# try qsmooth instead of percentile?
library(qsmooth)
# make a matrix for qsmooth input
melted_allNorm<-melt(allNormcounts_dt, id.vars = c("commonID","dataset"),
                     measure.vars=c("CTL_mean","hnRNPL_dep_mean"),
                     variable.name="genotype", value.name="normCounts")
melted_allNorm[genotype=="CTL_mean",forColnames:=paste(dataset,"CTL",sep="_"),by=.I]
melted_allNorm[genotype=="hnRNPL_dep_mean",forColnames:=paste(dataset,"Ldep",sep="_"),by=.I]
melted_allNorm<-merge(melted_allNorm,unique(allNormcounts_dt[,.(dataset,mouse_or_human)]))
melted_allNorm<-melted_allNorm[order(forColnames)]

for_qs_mat<-as.matrix(dcast(melted_allNorm, formula= commonID~forColnames, value.var="normCounts"),rownames=T)
# for_qs_mat_group<-unique(melted_allNorm[,.(dataset,forColnames)])[,dataset]
for_qs_mat_group<-unique(melted_allNorm[,.(mouse_or_human,forColnames)])[,mouse_or_human]
# make NAs zero
for_qs_mat[is.na(for_qs_mat)]<-0
qs_obj<-qsmooth(for_qs_mat,for_qs_mat_group)
qs_data<-qsmoothData(qs_obj)
qsmoothPlotWeights(qs_obj)

saveRDS(qs_data,"./R_rds_files/allQsmoothNormcounts_mat.rds")
qs_data<-readRDS("./R_rds_files/allQsmoothNormcounts_mat.rds")

get_qsmoothDT_forPlot<-function(cur_commonID){
  cur_vec<-qs_data[cur_commonID,]
  cur_dt<-data.table(qsNormCounts=cur_vec,forColnames=names(cur_vec))
  cur_dt[grep("CTL",forColnames),genotype:="CTL"]
  cur_dt[grep("Ldep",forColnames),genotype:="hnRNPL_dep"]
  cur_dt[genotype=="CTL",dataset:=sub("_CTL","",forColnames)]
  cur_dt[genotype=="hnRNPL_dep",dataset:=sub("_Ldep","",forColnames)]
  cur_dt<-merge(cur_dt,unique(allNormcounts_dt[,.(dataset,mouse_or_human)]))
  common_name<-allNormcounts_dt[commonID==cur_commonID,paste(unique(gene_name),collapse="|")]
  cur_dt[,common_name:=common_name]
  return(cur_dt)
}

forPlot_dt<-get_qsmoothDT_forPlot("38727") #hnRNPL
forPlot_dt<-get_qsmoothDT_forPlot("36496") #GPBP1
forPlot_dt<-get_qsmoothDT_forPlot("42830") #ZRANB2
forPlot_dt<-get_qsmoothDT_forPlot("37776") #SIRT1


plot_qsNorm_expr<-function(cur_commonID, saveToPDF=FALSE){
  cur_dt<-get_qsmoothDT_forPlot(cur_commonID)
  plot_name<-cur_dt[,unique(common_name)]
  require(ggplot2)
  require(egg)
  p1<-ggplot(data=cur_dt, aes(x=genotype, y=log2(qsNormCounts), group=dataset, color=dataset, shape=mouse_or_human))+
    geom_point()+
    geom_line()+
    ggtitle(plot_name)+
    scale_y_continuous("Quantile normalized\nexpression levels (log2)")+
    scale_shape_manual(values=c(20,18))+
    scale_color_manual(values=colorsForDatasets)+
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
          panel.border=element_rect(fill=NA,size=0.25),
          panel.spacing=unit(0,"line"))
  grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(1.5,"cm"),height=unit(2.5,"cm")))

  if(saveToPDF==T){
    output_folder<-"./figures/"
    pdf_name<-paste0(output_folder, "qsNormExprLevels_", sub("\\|.*","",plot_name), ".pdf")
    pdf(pdf_name,useDingbats=F,height=5,width=5)
    grid.arrange(grobs=lapply(list(p1),set_panel_size,width=unit(1.1,"cm"),height=unit(2.5,"cm")))
    dev.off()
  }
}

plot_qsNorm_expr("38727", saveToPDF=T) #hnRNPL
plot_qsNorm_expr("36496", saveToPDF=T) #GPBP1
plot_qsNorm_expr("42830", saveToPDF=T) #ZRANB2
plot_qsNorm_expr("37776", saveToPDF=T) #SIRT1
plot_qsNorm_expr("44924", saveToPDF=T) #hnRNPLL
plot_qsNorm_expr("41748", saveToPDF=T) #CDKN1A
plot_qsNorm_expr("43667", saveToPDF=T) #Myc
plot_qsNorm_expr("38269", saveToPDF=T) #E2F1
plot_qsNorm_expr("35837", saveToPDF=T) #E2F2
plot_qsNorm_expr("39834", saveToPDF=T) #E2F3
plot_qsNorm_expr("58860", saveToPDF=T) #E2F4
plot_qsNorm_expr("43022") #E2F5
plot_qsNorm_expr("49281") #E2F6
plot_qsNorm_expr("48455") #E2F7

plot_qsNorm_expr("39483", saveToPDF=T) #NSD2
plot_qsNorm_expr("41099", saveToPDF=T) #EPC1
plot_qsNorm_expr("45369", saveToPDF=T) #KDM6A
plot_qsNorm_expr("49981", saveToPDF=T) #BPTF
plot_qsNorm_expr("53387", saveToPDF=T) #MORF4L1

plot_qsNorm_expr("38235") #HM13
plot_qsNorm_expr("41936") #hnRNPR
plot_qsNorm_expr("42612") #RBM39
plot_qsNorm_expr("43591") #SMPD4
plot_qsNorm_expr("44005") #SEC31A

#------------------------------#
#getting intersect values like in upset plot. can eventually make a function if necessary
datasets
SE_list

#getting common diff spliced hits in all mouse datasets
interested_sets<-c("actb_d1","fetalLiv","thymocytes")
other_sets<-setdiff(datasets,interested_sets)

temp<-Reduce(intersect,SE_list[interested_sets])
temp2<-setdiff(temp,unlist(SE_list[other_sets]))
length(temp2)
merged_allSig_SE[commonSpliceID%in%temp2,sort(unique(geneSymbol))]
allMouse_commonIDs<-merged_allSig_SE[commonSpliceID%in%temp2,sort(unique(commonID))]

#getting common diff spliced hits in all human datasets
interested_sets<-setdiff(datasets,c("actb_d1","fetalLiv","thymocytes"))
other_sets<-setdiff(datasets,interested_sets)

temp<-Reduce(intersect,SE_list[interested_sets])
temp2<-setdiff(temp,unlist(SE_list[other_sets]))
length(temp2)
merged_allSig_SE[commonSpliceID%in%temp2,sort(unique(geneSymbol))]
allHuman_commonIDs<-merged_allSig_SE[commonSpliceID%in%temp2,sort(unique(commonID))]

#getting common diff spliced hits in mouse and human datasets, separately, but not exactly the same splice event
intersect(allHuman_commonIDs,allMouse_commonIDs)
plot_incLevel("44864_1") #UBAP2L
plot_incLevel("44864_11") #Ubap2l
# checked the sequence for these events.
# the current and upstream exons are identical between human and mouse.
# But the downstream exon is not the same, hence not called as homologous.
# Regardless, the included exon makes a +17aa insertion near the end of the protein, whose full length is ~1087aa.

# make a list for GO analysis. To test the suspicion that there are a lot of epigenetic modifiers here..
commonSpeciesSpecific_IDs<-commonID_forHomologs_dt[commonID%in%c(allHuman_commonIDs,allMouse_commonIDs), unique(human_id_noVer)]

cat(sub("\\..*","",commonSpeciesSpecific_IDs), sep="\n")
