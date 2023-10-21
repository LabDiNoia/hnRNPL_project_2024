library(rtracklayer)
library(data.table)

# step1: make common IDs for mouse and human genes
mouse_gtf<-import("../../common/mm10/gencode.vM23.primary_assembly.annotation.gtf")
human_gtf<-import("../../common/hg38/gencode.v43.primary_assembly.annotation.gtf")

mouse_dt<-data.table(mouse_id=unique(mouse_gtf$gene_id))
mouse_dt[,mouse_id_noVer:=sub('\\.[0-9]*','',mouse_id,perl=T),by=.I]

human_dt<-data.table(human_id=unique(human_gtf$gene_id))
human_dt[,human_id_noVer:=sub('\\.[0-9]*','',human_id,perl=T),by=.I]

mouse_human_homolog_dt<-fread("human_mouse_orthologs.txt")
colnames(mouse_human_homolog_dt)<-c("mouse_id_noVer","human_id_noVer")
mouse_human_homolog_dt<-mouse_human_homolog_dt[human_id_noVer!=""]
# only 25719 genes have homologs in mouse and humans.
# verified by getting human->mouse homologs too. i get the same number of genes

commonID_forHomologs_dt<-merge(mouse_dt,mouse_human_homolog_dt,all.x=T)
commonID_forHomologs_dt<-merge(human_dt,commonID_forHomologs_dt,by="human_id_noVer",all=T)
commonID_forHomologs_dt[,commonID:=.I]
commonID_forHomologs_dt
saveRDS(commonID_forHomologs_dt,"commonID_forHomologs_dt.rds")
library(data.table)
commonID_forHomologs_dt<-readRDS("commonID_forHomologs_dt.rds")

# step2: get gene_ids of differently spliced from the different datasets; translate it to commonID; do upset plot
#folders<-c("../HEK293T/","../actb_d1/","../bjFibro/","../c2c12/","../fetalLiv/","../hepg2/","../k562/","../kcytes/","../lncap/","../thymocytes/")
# removed c2c12 because only one rep: sketchy stats
folders<-c("../HEK293T/","../bjFibro/","../hepg2/","../k562/","../kcytes/","../lncap/","../thymocytes/","../fetalLiv/","../actb_d1/")
datasets<-basename(folders)
splicing_type<-c("SE")

get_diff_spl_genes<-function(dataset,alt_splice_tpye){
	# path to file
	dataset_folder<-paste0("../",dataset)
	file_path<-paste0(dataset_folder,"/rna_seq/hnrnpl/results/rMATS/",alt_splice_type,".MATS.JCEC.txt")
	# get file
	rmats_dt<-fread(file_path)
	# filter
	#rmats_dt<-rmats_dt[FDR<0.1]
	rmats_dt<-rmats_dt[PValue<0.01]
	# get geneIDs and convert to commonIDs
	commonID_vector<-commonID_forHomologs_dt[human_id %in% rmats_dt[,GeneID] | mouse_id %in% rmats_dt[,GeneID], commonID]
	return(commonID_vector)
}

SE_list<-sapply(datasets,get_diff_spl_genes,"SE")
saveRDS(SE_list,"SE_list.rds")
SE_list<-readRDS("SE_list.rds")

# common ones
Reduce(intersect,SE_list)
# there are 16 hits here!
hitsInAll<-commonID_forHomologs_dt[commonID %in% Reduce(intersect,SE_list)]
cat(hitsInAll[,human_id_noVer],sep="\n")

# for the common ones, print the splicing position information by species?
get_spl_info<-function(dataset,alt_splice_tpye,common_gene_id){
        # path to file
        dataset_folder<-paste0("../",dataset)
        file_path<-paste0(dataset_folder,"/rna_seq/hnrnpl/results/rMATS/",alt_splice_tpye,".MATS.JCEC.txt")
        # get file
        rmats_dt<-fread(file_path)
        # filter
	rmats_dt<-rmats_dt[PValue<0.01]
	# get human and mouse ids
	relevant_ids<-commonID_forHomologs_dt[commonID==common_gene_id,c(human_id,mouse_id)]
	# get relevant lines
        relevant_rmats_dt<-rmats_dt[GeneID %in% relevant_ids]
	return(relevant_rmats_dt)
}
get_spl_info("HEK293T","SE","37776")

SE_hitsInAll_info<-sapply(datasets,get_spl_info,"SE",hitsInAll[1,commonID],simplify=F)
# tried for all 16 hits
# the conserved ones are: Sirt1, GPBP1 and UBAP2L
####################################################################################3
# goal: upset plot with splicing events not genes

# 1. get all splicing events with PValue<0.01 in one table
get_spl_dt<-function(dataset,alt_splice_tpye){
        # path to file
        dataset_folder<-paste0("../",dataset)
        file_path<-paste0(dataset_folder,"/rna_seq/hnrnpl/results/rMATS/",alt_splice_tpye,".MATS.JCEC.txt")
        # get file
        rmats_dt<-fread(file_path)
        # filter
        rmats_dt<-rmats_dt[PValue<0.01]
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
# 2. assign IDs to splicing events ??
# keep only columns till "downstreamEE"
uniqSig_SE<-unique(allSig_SE[,1:match("downstreamEE",colnames(allSig_SE))])
# assign spliceEvent ID
uniqSig_SE[,spliceID:=paste(commonID,1:.N,sep="_"),by=commonID]
# is the homolog also differentially spliced?
uniqSig_SE[,homolog_diffSpliced:=length(unique(mouse_or_human))==2,by=commonID]

# difficult bit now...
# if mouse and human sequences are homologous, assign same spliceID...
# number of genes differentially spliced in both human and mouse in any dataset is 3340
length(unique(uniqSig_SE[homolog_diffSpliced==T,commonID])) # 3340

#testing...............
# libraries to get sequences
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)

# work by commonID. get one eg commonID
test<-uniqSig_SE[commonID=="53305"]

get_alnInfo<-function(DTforGR,whichExon){
	# define exonStartEndColDefText variable. will be used to define start-end
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
		makeGRangesFromDataFrame(test[mouse_or_human=="human", .(chr,strand,',
		exonStartEndColDefText,
		',spliceID)], keep.extra.columns=T)')
	eval(parse(text=human_gr_expr))

	mouse_gr_expr<-paste0('mouse_gr<-
		makeGRangesFromDataFrame(test[mouse_or_human=="mouse", .(chr,strand,',
		exonStartEndColDefText,
		',spliceID)], keep.extra.columns=T)')
	eval(parse(text=mouse_gr_expr))

	if(whichExon=="current"){
		assign("human_gr", human_gr, envir=.GlobalEnv)
		assign("mouse_gr", mouse_gr, envir=.GlobalEnv)
		aln_nmatch_prop<-list()
		assign("aln_nmatch_prop", aln_nmatch_prop, envir=.GlobalEnv)
	}

	# get sequences
	human_seqs<-Views(Hsapiens,human_gr)
	mouse_seqs<-Views(Mmusculus,mouse_gr)

	# getting pairwiseAlingment score matrix. maybe think about changing gapExtension
	exonAlnInfo_mat<-sapply(human_seqs,function(x){
		pairwiseAlignment(pattern=mouse_seqs,subject=x,type="global")})
	
	# assign to global enviroment, nmatch / length value for each comparison
	maxlen_mat<-outer(width(mouse_gr),width(human_gr),FUN=pmax)
	# in outer, the order is important!! mouse_gr then human_gr for consistency
	nmatch_prop_mat<-sapply(exonAlnInfo_mat,nmatch) / maxlen_mat
	colnames(nmatch_prop_mat)<-human_gr$spliceID
	rownames(nmatch_prop_mat)<-mouse_gr$spliceID
	nmatch_prop_list<-list(nmatch_prop_mat)
	names(nmatch_prop_list)<-whichExon
	aln_nmatch_prop<-append(aln_nmatch_prop, nmatch_prop_list)
	assign("aln_nmatch_prop", aln_nmatch_prop, envir=.GlobalEnv)

	return(exonAlnInfo_mat)
}

extractAlnDetails<-function(alnInfo_matrixList,functionToApply){
		sapply(names(alnInfo_matrixList),function(x){
			matList<-sapply(alnInfo_matrixList[[x]],functionToApply)
			colnames(matList)<-human_gr$spliceID
			rownames(matList)<-mouse_gr$spliceID
			return(matList)}, simplify=F)
}

test<-uniqSig_SE[commonID=="35622"]
test<-uniqSig_SE[commonID=="37776"] #sirt1
test<-uniqSig_SE[commonID=="44864"] #ubap2l
test<-uniqSig_SE[commonID=="36496"] #gpbp1

alnInfo_matList<-sapply(c("current","upstream","downstream"),function(x){
			get_alnInfo(test,x)}, simplify=F)
aln_nmatch_prop

#extractAlnDetails(alnInfo_matList,pid)
#extractAlnDetails(alnInfo_matList,score)

# get percentage identity of the sequences
alnPID<-extractAlnDetails(alnInfo_matList,pid)
## get a list of matrices that are truth tables for PID>70% for each exon
#sapply(alnPID,function(x)x>70,simplify=F)
# reduce to get truth table for PID>70% for all three exons
aln_over70p_TT<-Reduce(f=function(x,y){x&y}, x=sapply(alnPID,function(x)x>70,simplify=F))
aln_over70p_TT
Reduce(f=function(x,y){x+y}, x=aln_nmatch_prop) * aln_over70p_TT
alnScore_homMat<-Reduce(f=function(x,y){x+y}, x=extractAlnDetails(alnInfo_matList,score)) * aln_over70p_TT

Reduce(f=function(x,y){x+y}, x=extractAlnDetails(alnInfo_matList,score))
Reduce(f=function(x,y){x+y}, x=aln_nmatch_prop)

# make a list of spliceIDs, then write commonSpliceIDs (homologous splice event). 
alnScore_mat_forExtr<-alnScore_homMat
alnScore_mat_forExtr[alnScore_mat_forExtr==0]=-Inf
while( sum(alnScore_mat_forExtr!=-Inf) > 0){
	max_ind<-which(alnScore_mat_forExtr==max(alnScore_mat_forExtr),arr.ind=T)
	human_spliceID<-colnames(alnScore_mat_forExtr)[max_ind[1,2]]
	mouse_spliceID<-rownames(alnScore_mat_forExtr)[max_ind[1,1]]
	test[spliceID==mouse_spliceID,commonSpliceID:=human_spliceID]
	alnScore_mat_forExtr<-alnScore_mat_forExtr[-max_ind[1,1],-max_ind[1,2]]
}
test[is.na(commonSpliceID),commonSpliceID:=spliceID,by=.I]



####################################################################################3
library(UpSetR)
upset(fromList(SE_list),nintersects=NA,nsets=10,empty.intersections = "on",order.by="freq")
upset(fromList(SE_list),nsets=10,empty.intersections = "on",order.by="degree")
upset(fromList(SE_list),nsets=10,order.by="degree")
dev.off()





