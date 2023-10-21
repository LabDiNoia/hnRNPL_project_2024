# convenience functions to compile lists of genes affected
library(data.table)

folders<-c("../HEK293T/","../bjFibro/","../hepg2/","../k562/","../kcytes/","../lncap/","../thymocytes/","../fetalLiv/","../actb_d1/")
datasets<-basename(folders)

#cur_dataset<-datasets[1]
dir.create("./R_rds_files/genes_with_splicing_changes/")

for(cur_dataset in datasets){
  # get file names
  cur_rMATS_folder_path<-paste0("../", cur_dataset, "/rna_seq/hnrnpl/results/rMATS/")
  splicing_files<-list.files(cur_rMATS_folder_path, pattern="*.MATS.JCEC.txt", full.names = T)
  names(splicing_files)<-sub(".MATS.JCEC.txt", "", basename(splicing_files))
  
  # keep only info relevant for now
  # filter for FDR<0.1 & IncLevelDifference>abs(0.1), keep only gene ID/name, FDR and IncLevelDifference 
  allSplice_dt<-lapply(splicing_files,
                            function(file_path){
                              require(data.table)
                              temp_dt<-fread(file_path)
                              temp_dt[,splicing_type:=sub(".MATS.JCEC.txt", "", basename(file_path))]
                              temp_dt<-temp_dt[FDR<0.1 & abs(IncLevelDifference)>0.1]
                              temp_dt[,.(splicing_type,GeneID,geneSymbol,FDR,IncLevelDifference)]
                            })
  allSplice_dt<-rbindlist(allSplice_dt)
  allSplice_dt[,gene_id_noVer:=sub("\\.[0-9]*","",GeneID)]
  
  #output genelist to a file
  cur_out_file<-paste0("./R_rds_files/genes_with_splicing_changes/",cur_dataset,"_allDiffSplicedGenes.txt")
  sink(cur_out_file)
  cat(allSplice_dt[,unique(gene_id_noVer)], sep="\n")
  sink()
  
  cat("Printed", length(allSplice_dt[,unique(gene_id_noVer)]), "geneIDs to", cur_out_file,"\n")
}

# summary - splicing numbers (gene and event counts)
splEventNum_summary_dt<-data.table()
splGenesNum_summary_dt<-data.table()

for(cur_dataset in datasets){
  # get file names
  cur_rMATS_folder_path<-paste0("../", cur_dataset, "/rna_seq/hnrnpl/results/rMATS/")
  splicing_files<-list.files(cur_rMATS_folder_path, pattern="*.MATS.JCEC.txt", full.names = T)
  names(splicing_files)<-sub(".MATS.JCEC.txt", "", basename(splicing_files))
  
  # keep only info relevant for now
  # filter for FDR<0.1 & IncLevelDifference>abs(0.1), keep only gene ID/name, FDR and IncLevelDifference 
  allSplice_dt<-lapply(splicing_files,
                       function(file_path){
                         require(data.table)
                         temp_dt<-fread(file_path)
                         temp_dt[,splicing_type:=sub(".MATS.JCEC.txt", "", basename(file_path))]
                         temp_dt<-temp_dt[FDR<0.1 & abs(IncLevelDifference)>0.1]
                         temp_dt[,.(splicing_type,GeneID,geneSymbol,FDR,IncLevelDifference)]
                       })
  allSplice_dt<-rbindlist(allSplice_dt)
  
  event_num_dt<-allSplice_dt[,.N,by=splicing_type]
  event_num_dt[,dataset:=cur_dataset]
  splEventNum_summary_dt<-rbind(splEventNum_summary_dt,event_num_dt)
  
  gene_num_dt<-unique(allSplice_dt[,.(GeneID,splicing_type)])[,.N,by=splicing_type]
  gene_num_dt[,dataset:=cur_dataset]
  splGenesNum_summary_dt<-rbind(splGenesNum_summary_dt,gene_num_dt)
  
}

write.table(dcast(splEventNum_summary_dt, formula = splicing_type~dataset, value.var = "N"), 
            "./R_rds_files/spliceEvents_count_summary.txt", sep="\t", row.names=F, quote=F)

write.table(dcast(splGenesNum_summary_dt, formula = splicing_type~dataset, value.var = "N"), 
            "./R_rds_files/spliceGenes_count_summary.txt", sep="\t", row.names=F, quote=F)
