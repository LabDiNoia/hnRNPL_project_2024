library(rtracklayer)
library(data.table)

# step0: make common IDs for mouse and human genes
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
commonSpliceID_lookup_dt<-data.table()
counter<-0
nhomologs<-length(commonIDs_homolog_diffSpliced)
for(cur_ID in commonIDs_homolog_diffSpliced){
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
		}
		cur_SE_dt[is.na(commonSpliceID),commonSpliceID:=spliceID,by=.I]
		commonSpliceID_lookup_dt<-rbind(commonSpliceID_lookup_dt,cur_SE_dt[,.(spliceID,commonSpliceID)])
	}
	counter<-counter+1
	cat("Done with",counter ,"of",nhomologs ,":", cur_ID,"(",cur_SE_dt[1,geneSymbol],"). ",sep=" ")
	cat("Took",Sys.time()-start_time,"seconds","\n")
	rm(cur_SE_dt)
}

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

test<-sapply(commonIDs_homolog_diffSpliced[1:20],get_commonSpliceID_lookup_dt,simplify=F)
test<-rbindlist(test)


saveRDS(commonSpliceID_lookup_dt,"commonSpliceID_lookup_dt_FDR_lt_pt2.rds")
commonSpliceID_lookup_dt<-readRDS("commonSpliceID_lookup_dt_FDR_lt_pt2.rds")

# merge commonSpliceID_lookup_dt and uniqSig_SE
temp<-copy(uniqSig_SE)
temp<-merge(temp,commonSpliceID_lookup_dt,by="spliceID",all.x=T)
temp[is.na(commonSpliceID),commonSpliceID:=spliceID,by=.I]
nrow(temp) #49655
length(unique(temp[,commonSpliceID])) # 48757

# merge uniqSig_SE with commonSpliceIDs and allSig_SE
merged_temp<-merge(allSig_SE, temp, by=colnames(allSig_SE)[ 1:match("downstreamEE",colnames(allSig_SE)) ])

# hits in at least 8 / 9 datasets
merged_temp[,commonSpliceID,by=dataset][,.N,by=commonSpliceID][N>=8]

# make list for upset plot
SE_list<-sapply(datasets,function(x){merged_temp[dataset==x,unique(commonSpliceID)]})

####################################################################################3
library(UpSetR)
#upset(fromList(SE_list),nintersects=NA,nsets=10,empty.intersections = "on",order.by="freq")
#upset(fromList(SE_list),nsets=10,empty.intersections = "on",order.by="degree")
upset(fromList(SE_list),nsets=10,order.by="degree")
dev.off()



