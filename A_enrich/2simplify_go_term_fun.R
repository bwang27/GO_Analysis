

shorten_go_terms.fun <-function(go_cat, go_sig_data,  p_value_name, go_tern_name="uniq_go_ids", pval_cut=0.01, spe="hs", overlap_cut=0.75,overlap_min=0.2){
  #delete redundant GOs in a list based on shared genes
  #if shared gene number of GO > cutoff (75%), delete less significant GO
  #go2gene_data, relationship of geneset and its associated genes
	#go2gene_file<-paste(root_dir,"data/go/",spe,"/",spe,".",go_cat,".go2gene.tbl",sep="")
	go2gene_file<- go_cat2_go2id_f(go_cat, spe)
	go2gene_data<-read.table(go2gene_file,sep="\t",header=T,comment="")
	names(go2gene_data)[1]<-"GO_ID"
	
	if_select_go <- go_sig_data[[p_value_name]]<pval_cut & !is.na(go_sig_data[[p_value_name]])
	if(sum(if_select_go,na.rm=T)<=3){
		return (as.numeric(if_select_go))
	}
	ana_GOs <- as.character(go_sig_data[[go_tern_name]][if_select_go])
	ana_GOs_pval <- go_sig_data[[p_value_name]][if_select_go]
	names(ana_GOs_pval) <- ana_GOs
	if_select_go_tmp <- rep(T,length(ana_GOs))
	
  go2geneids_list <- NULL
  for(ana_1go in ana_GOs){
  	genes <- as.character(unique(go2gene_data$Gene_ID[go2gene_data[,"GO_ID"]==ana_1go]))
	  if(length(genes)>0)  {
	     go2geneids_list <- c(go2geneids_list, list(genes))
	     names(go2geneids_list)[length(go2geneids_list)] <- ana_1go
	  }
  }
  go2geneids_list_sizes <- sapply(go2geneids_list,length)
  
	judge_compare2go<-function(go1,go2,overlap_cut,overlap_min){
		#go1 <- ana_GOs[i]
		#go2 <- ana_GOs[j]
		common_genes_num <- length( intersect(go2geneids_list[[go1]],go2geneids_list[[go2]]) )
		overlap_in_go1<-common_genes_num/go2geneids_list_sizes[go1]
		overlap_in_go2<-common_genes_num/go2geneids_list_sizes[go2]
		if(max(c(overlap_in_go1,overlap_in_go2),na.rm=T)>overlap_cut & min(c(overlap_in_go1,overlap_in_go2),na.rm=T)>overlap_min){
			pval_magn_diff <- log10(ana_GOs_pval[go1]/ana_GOs_pval[go2])
			if(ana_GOs_pval[go1]==ana_GOs_pval[go2]){
				return( ifelse(go2geneids_list_sizes[go1] < go2geneids_list_sizes[go2], go2, go1) )
			}else if(abs(pval_magn_diff)< 0.5){ #similar significance
				return( ifelse(go2geneids_list_sizes[go1] < go2geneids_list_sizes[go2], go2, go1) ) #filter larger one
			}else{
				return( ifelse(pval_magn_diff < 0, go2, go1) ) #filter less significant one
			}
		}else{
			return(0)
		}
	}
  
  removed_gos <- NULL
  retained_gos <- NULL
  for(i in 1:(length(ana_GOs)-1)){
  	for (j in (i+1):length(ana_GOs)){
  		go1 <- ana_GOs[i]; go2 <- ana_GOs[j]; 
  		if(any(ana_GOs[c(i,j)] %in% removed_gos)){next}
  		removed_go <- judge_compare2go(go1,go2, overlap_cut,overlap_min)
  		if(removed_go!=0){
  			removed_gos <- c(removed_gos,removed_go)
  			retained_gos <- c(retained_gos, setdiff(ana_GOs[c(i,j)], removed_go))
  		}
  	}
  }
  #print(removed_gos)
  #print(retained_gos)
  #if_select_go[go_sig_data[[go_tern_name]] %in% removed_gos] <- FALSE
  if_select_go <- as.numeric(if_select_go)
  if_select_go[ match(removed_gos,go_sig_data[[go_tern_name]]) ] <- retained_gos
  return(if_select_go)
}

ano_gids_withGO <- function(GoIds, spe="hs", go_cat=NULL, GoDescs=NULL){ #given a list of GO ids, anotate gene ids with these go ids
	#spe="hs" or "mm"
	go_cats <- c("bp","cc","mf")
	#GoIds<-c("GO:0042742","GO:0006875","GO:0006941")
	if(!is.null(go_cat)){go_cats<-go_cat}
	go2gene_data<- NULL
	go_anot_data<- NULL
	for(go_cat in go_cats){
		go2gene_file<-paste(root_dir,"data/go/",spe,"/",spe,".",go_cat,".go2gene.tbl",sep="")
		go2gene_data<-rbind(go2gene_data, read.table(go2gene_file,sep="\t",header=T,comment=""))
		if(is.null(GoDescs)){
			go2_anot_file <- sub("go2gene.tbl","go2gene.num.tbl",go2gene_file)
			tmp_D <- read.table(file=go2_anot_file, sep="\t", header=T,quote="",comment.char = "")
			names(tmp_D)[1]<-"GO_ID" ##GO_ID  GO_term Num_Child       Num_Parent      Num_Genes
			go_anot_data <- rbind(go_anot_data, tmp_D)
		}
	}
	names(go2gene_data)[1]<-"GO_ID"  #GO_ID  Gene_ID
	go2gene_data2 <- go2gene_data[go2gene_data$GO_ID %in% GoIds,]
	if(is.null(GoDescs)){
		go2gene_data2["GoDescs"] <- go_anot_data$GO_term[match(go2gene_data2$GO_ID, go_anot_data$GO_ID)]
	}else{
		go2gene_data2["GoDescs"] <- GoDescs[match(go2gene_data2$GoIds, GoIds)]
	}
	#for every gene id, combine all required GO anotation
	Gene2GoDesc <- tapply(paste(go2gene_data2$GO_ID,go2gene_data2$GoDescs,sep=":"), go2gene_data2$Gene_ID, paste, collapse="||")
}


#overlap_cut <- 0.75 
#overlap_min <- 0.2
#pval_cut <- 0.05
#desribe_str<-paste("#overlap_cut =",overlap_cut,", overlap_min =",overlap_min,", pval_cut =",pval_cut)
#
#go_cats <- c("bp","cc","mf")
#


