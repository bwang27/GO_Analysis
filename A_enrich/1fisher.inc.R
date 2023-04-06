doFisher41GoExtype <- function(num_vec, num_total_type1, num.pop){
	num_vec[is.na(num_vec)]=0
	num_total_go <- as.numeric(num_vec[1])
	num_go_and_type1 <- as.numeric(num_vec[2])
	if(num_go_and_type1==0) {return(1)}else{
		fisher_matrix <- matrix(c(num_go_and_type1, num_total_type1-num_go_and_type1, num_total_go-num_go_and_type1, num.pop-num_total_go-num_total_type1+num_go_and_type1),ncol=2)
		return(fisher.test(fisher_matrix,alternative="g")$p.value)
	}
}

cal_strsplit_num <-function(str_in, sep_str){
	str_vec <- unlist(strsplit(str_in, sep_str))
	length(str_vec)
}

do_1go_ana <- function(genetb_data2,go_cat){
	############################################################
	#GO associated genes
	GO_tbl <- read.table(file=go2geneid_file, sep="\t",header=T)
	colnames(GO_tbl)[c(1,2)] <- c("GO", "Gene") ; nrow(GO_tbl)
	go_anot_data <- read.table(file=go2ano_file, sep="\t", header=T,quote="",comment.char = "")
	names(go_anot_data)=change_values(names(go_anot_data), c("STANDARD_NAME","DESCRIPTION_BRIEF"), c("GO_ID","GO_term") )
	if("MEMBERS_SYMBOLIZED" %in% names(go_anot_data) & !("Num_Genes" %in% names(go_anot_data)) ){
		go_anot_data["Num_Genes"] <- sapply(as.character(go_anot_data$MEMBERS_SYMBOLIZED),cal_strsplit_num,",")
	}
	if("Category" %in% names(go_anot_data)){
		GO_tbl<-GO_tbl[ GO_tbl$GO %in% go_anot_data$GO_ID[go_anot_data$Category==go_cats2[go_cat]], ]
	}
	go_genes <- unique(GO_tbl[,2]) ; length(go_genes)
	
	if(!is.null(gene_groups)){
		genetb_data2 <- genetb_data2[genetb_data2[[gene_grp_name]] %in% gene_groups,]
	}
	sel_genes <- intersect(genetb_data2$gene_id, go_genes) ; 
	print (paste("final selected gene number =",length(sel_genes)))
	genetb_data2 <- genetb_data2[genetb_data2$gene_id %in% sel_genes,]; nrow(genetb_data2)
	genetb_data2 <- genetb_data2[!duplicated(genetb_data2$gene_id),]
	gene_grp_types <- sort(unique(genetb_data2[[gene_grp_name]]))
	gene_grp_types <- setdiff(gene_grp_types,discard_grps)
	
	print (paste(c("analyzed gene groups =",gene_grp_types), collapse=" "))

	##link GO -> gene id -> up/dn numbers
	GO_tbl <- GO_tbl[GO_tbl$Gene %in% sel_genes,]; nrow(GO_tbl)
	uniq_go_ids <- as.character(unique(GO_tbl$GO))
	
	exp_ch_in_go_tb <- genetb_data2[[gene_grp_name]][match(GO_tbl$Gene,genetb_data2$gene_id) ]
	
	print("calculate gene numbers associated with go terms")
	Num_Exp <- table(GO_tbl$GO)[uniq_go_ids]
	out_data <- transform(cbind(uniq_go_ids,Num_Exp))
	for (gene_grp_type in gene_grp_types){
		out_data <- cbind(out_data, table(GO_tbl$GO[exp_ch_in_go_tb == gene_grp_type])[uniq_go_ids])
		names(out_data)[ncol(out_data)] <- paste("Num_",gene_grp_type,sep="")
	}
	print (paste("finally analyzed GO term:",nrow(out_data)))
	
	############################################################
	num.pop <- length(sel_genes) #all expressed gene (have go) number
	out_str <- paste("#Numbers: total expressed genes with GO = ",num.pop,sep="")
	
	
	print ("do fisher's exact test for every selected GO")
	for (gene_grp_type in gene_grp_types){
		num_total_type1 <- sum(genetb_data2[[gene_grp_name]] == gene_grp_type)
		print (paste("number of genes to be",gene_grp_type,"=", num_total_type1))
		out_str <- paste(out_str, ", ",gene_grp_type,"=",num_total_type1,sep="")
		tmpdata <- apply(out_data[,c("Num_Exp",paste("Num_",gene_grp_type,sep=""))], 1, doFisher41GoExtype, num_total_type1,  num.pop)
		out_data <- cbind(out_data,tmpdata)
		names(out_data)[ncol(out_data)] <- paste("P_",gene_grp_type,sep="")
	}
	write(out_str,file=out_summary_file, append=T)
	
	P_Min=apply(out_data[grep("P_",names(out_data),value=T)], 1, min, na.rm=T)
	P_Max=apply(out_data[grep("P_",names(out_data),value=T)], 1, max, na.rm=T)
	P_magnitude_diff <- round(log10(P_Max)-log10(P_Min),1)
	out_data <- cbind(out_data, P_Min, P_magnitude_diff)
	out_data <- out_data[sort.list(out_data$P_Min),]
	
	############################
	##add go anotation
	print ("add anotation to every GO")
	add_ano_headers=intersect(c("GO_term","Num_Child","Num_Parent","Num_Genes"), names(go_anot_data))
	out_data <- cbind(out_data, go_anot_data[match(out_data$uniq_go_ids, go_anot_data$GO_ID), add_ano_headers ])
	
	#######for each significant gene set, generate randomized pseudo-gene set with same gene number and similar expression level, calculate distribution P value for randomized datasets, calculate Q value)
	#GO term -> gene ids -> expression distribution -> select similar expression pattern gene set
	
	if(ifgen_randomized_data){
		ifsel_rows=out_data$P_Min<Ra_pmin_cut
		if("Num_Child" %in% names(out_data) ){ ifsel_rows=ifsel_rows & out_data$Num_Child<go_childnum_cut; }
		sel_rowids<- (1:nrow(out_data))[ ifsel_rows ];
		out_data=do_randomSampling(out_data, genetb_data2, gene_grp_name,gene_grp_types, sel_rowids, GO_tbl)
	}#end if(ifgen_randomized_data)
	
	out_data
	############################
}

do_randomSampling<-function(out_data, genetb_data2, gene_grp_name,gene_grp_types, sel_rowids, GO_tbl, num.pop=nrow(genetb_data2), sample_name=NULL ){
		gene_grp_nums<- table(genetb_data2[[gene_grp_name]])
		Q_names <- paste("Q_",gene_grp_types,sep="")
		#Q5q_names <- paste("Q5q_",gene_grp_types,sep="") #0.05 quantile of randomized p value
		if(!is.null(sample_name)){
			Q_names<-paste(Q_names,sample_name,sep='.')
		}
		P_names <- sub("^Q_","P_",Q_names)
		print(P_names)
		if_P_logged=F
		ori_P_tb=out_data[P_names]
		if(any( abs(out_data[P_names])>1 ) ){
			if_P_logged=T
			ori_P_tb=10^(-abs(ori_P_tb)) #change to un-logged
		}
		print(paste("#random sampling for", length(sel_rowids), "GOs"))
		add_Q_names=setdiff(Q_names,names(out_data))
		if(length(add_Q_names)>0){
			out_data[add_Q_names]<-rep(NA, nrow(out_data)*length(add_Q_names))
		}
		
		#out_data[Q5q_names]<-rep(NA, nrow(out_data)*length(Q_names))
		#out_data$gex_bin_median<-NA
		GO_ID_name=intersect( c('uniq_go_ids','GO_ID'), names(out_data) )[1]
		if(is.na(GO_ID_name)){print("Error: herder uniq_go_ids or GO_ID not found in GO table!"); exit;}
		for( rowi in sel_rowids ){
			if( all(!is.na(out_data[rowi,Q_names])) ){
				next
			}
			GO_ID1<-as.character(out_data[rowi,GO_ID_name])
			if_inGset<- genetb_data2$gene_id %in% GO_tbl$Gene[GO_tbl$GO==GO_ID1]; table(if_inGset)
			num_total_go <- sum(if_inGset)
			print(paste("randomize for row",rowi, "num.pop=",num.pop, "num_total_go=",num_total_go))
			#out_data$gex_bin_median[rowi] <- median(genetb_data2$gex_bin[if_inGset])
			gex_bin_nums <- table(genetb_data2$gex_bin[if_inGset])[as.character(1:RaGex_binnum)]
			gex_bin_nums[is.na(gex_bin_nums)]<-0
			names(gex_bin_nums) <- as.character(1:RaGex_binnum)
			Ra_Pvals<- sapply(1:Ra_times, function(Ra_i){#random sampling one time
				sample_rowids<- sample( (1:nrow(genetb_data2)), sum(if_inGset), prob=gex_bin_nums[as.character(genetb_data2$gex_bin)] )
				if_in_Ra_gset<-rep(FALSE,nrow(genetb_data2))
				if_in_Ra_gset[sample_rowids]<-TRUE
				unlist(sapply(gene_grp_types,function(gene_grp_type){
					num_total_type1<-gene_grp_nums[gene_grp_type]
					num_go_and_type1<- sum(genetb_data2[,gene_grp_name]==gene_grp_type & if_in_Ra_gset)
					doFisher41GoExtype(c(num_total_go, num_go_and_type1), num_total_type1,num.pop)
				}))
			})
			Ra_Pvals<-matrix(unlist(Ra_Pvals),nrow=length(gene_grp_types))
			out_data[rowi,Q_names]<-unlist(sapply(1:length(gene_grp_types), function(type_i){
				sum(Ra_Pvals[type_i,]<=ori_P_tb[rowi,P_names[type_i]],na.rm=T)/Ra_times
			}))
			#out_data[rowi,Q5q_names]<-unlist(sapply(1:length(gene_grp_types), function(type_i){
			#	quantile(Ra_Pvals[type_i,], probs=0.05, na.rm=T)
			#}))
			print("done")
		}#end for(rowi in sel_rowids)
		out_data	
}


do_all_go_ana <- function(genetb_data1,organism,ana_summary_str,outfilebase){
	write(paste("#",ana_summary_str,sep=""), file=out_summary_file)
	for(go_cat in go_cats){
		print("")
		print (paste("GO group =",go_cat))
		go2geneid_file <<- go_cat2_go2id_f(go_cat, organism)
		go2ano_file <<- go_cat2_go2ano_f(go_cat, organism)
		############################################################
		write(paste("\n\n#go=",go_cat,sep=""), file=out_summary_file, append=T)
		out_data <- do_1go_ana(genetb_data1,go_cat)
		write.table(out_data[out_data$P_Min<out_pval_cut,], file=out_summary_file, col.names=T, row.names=FALSE, sep="\t", quote=FALSE, append=T)
	}
}

go_cat2_go2id_f <- function(go_cat, organism){
	#paste(go_dir,organism,"/",organism,".",go_cat,".go2gene.tbl", sep="")
	if(go_cat %in% c("bp","cc","mf") ){ #from NCBI GO
		paste(go_dir,organism,".go2gene.tbl", sep="")  ###
	}else if( go_cat %in% c("c3.tft","c3.mir") ){ #from msigdb
		paste(go_dir,"1Msigdb/",organism,"/",go_cat,".",gset_version,".gid2gset.tb", sep="")
	}
}
go_cat2_go2ano_f <- function(go_cat, organism){
	#paste(go_dir,organism,"/",organism,".",go_cat,".go2gene.num.tbl", sep="")	###
	if(go_cat %in% c("bp","cc","mf") ){ #from NCBI GO
		paste(go_dir,organism,".go2gene.num.tbl", sep="")	###
	}else if( go_cat %in% c("c3.tft","c3.mir") ){ #from msigdb
		paste(go_dir,"2Msigdb_ano/","msigdb_",gset_version,".ano.tbl", sep="")
	}
}

pre_process_data <- function(organism,genetb_data, gene_grp_name=NULL, gene_grp_names=NULL, comb_grp_l=NULL, gene_groups=NULL, filter_names=NULL, filter_val_list=NULL,
	auto_setup=NULL, rm_redu_gid_sort_var=NULL, if_sort_decr=F, if_only_remove_old=FALSE, only_top_reg_num=NULL, top2reg_names=NULL, sel_top_sort_var=NULL,
	ifgen_randomized_data=FALSE, RanGex_vars=NULL, RaGex_binnum=10){
	
	if(gene_id_var!="gene_id"){
		genetb_data$gene_id<-genetb_data[,gene_id_var]
	}

	if("gene_grp_names" %in% ls() ){ #combine several column as gene group 
		if( length(gene_grp_names)>1){
			print(paste(c('combine several columns (', gene_grp_names, ')TO',gene_grp_name), collapse=' '))
			print(genetb_data[1,])
			genetb_data[gene_grp_name] <- apply(genetb_data[gene_grp_names],1,paste,collapse="_")
			print(table(genetb_data[gene_grp_name]))
	}}
	
	if("comb_grp_l" %in% ls()){
		if(!is.null(comb_grp_l)){
			print ("use comb_grp_l to combine values in one column to one group:")
			for (new_grp_name in names(comb_grp_l)){
				genetb_data[[gene_grp_name]] [genetb_data[[gene_grp_name]] %in% comb_grp_l[[new_grp_name]]] <- new_grp_name
				print (c(comb_grp_l[[new_grp_name]], "->", new_grp_name))
			}
			print(table(genetb_data[gene_grp_name]))
		}
	}
	if(!is.null(gene_grp_name)){
		genetb_data <- genetb_data[!is.na(genetb_data[,gene_grp_name]), ]
	}
	
	if(!is.null(gene_groups)){
		genetb_data <- genetb_data[genetb_data[,gene_grp_name] %in% gene_groups,]
		print(table(genetb_data[gene_grp_name]))
	}
	
		##update refgex_d #forLVH_Seq only
		#refgex_d[,gene_grp_name] <- genetb_data[gene_grp_name]
	
	if(!is.null(filter_names) ){
		if( length(filter_names)>0){
			print("use filter:")
			if(!is.null(filter_unSelVal_list) ){ #not selected values
				names(filter_unSelVal_list) <- filter_names
				for(filter_name in filter_names){
					genetb_data <- genetb_data[!(genetb_data[[filter_name]] %in% filter_unSelVal_list[[filter_name]]), ]
					print (paste(c(filter_name,"!=", filter_unSelVal_list[[filter_name]])))
					print(paste("row# after filtring: ", nrow(genetb_data) ))
				}
			}else{ #use filter_val_list
				names(filter_val_list) <- filter_names
				for(filter_name in filter_names){
					genetb_data <- genetb_data[genetb_data[[filter_name]] %in% filter_val_list[[filter_name]],]
					print (paste(c(filter_name,"=", filter_val_list[[filter_name]])))
					print(paste("row# after filtring: ", nrow(genetb_data) ))
				}
			}
			
			
			if(!is.null(gene_grp_name)){ print(table(genetb_data[gene_grp_name])) }
			filterstr <<- paste(paste(filter_names,unlist(lapply(filter_val_list,paste,collapse="")), sep="="), collapse="; ")
		}else{
			filterstr <<- NULL
		}
	}else{
		filterstr <<- NULL
	}
	
	if(!is.null(auto_setup)){
		print ("use auto_setup")
		if(auto_setup=="pA_gexch"){ #for a table of polyA type and gene  expression change (also used in cis elements study)
		genetb_data <- read.table(genetb_file,sep="\t",header=T)
		genetb_data <- genetb_data[genetb_data$pA_type %in% c("3S","3F"),]
		genetb_data <- genetb_data[genetb_data$gex_regu_type %in% c("up","dn"),] ###
		genetb_data["pA_gexreg_type"] <- paste(genetb_data$pA_type, genetb_data$gex_regu_type, sep="_")
		gene_grp_name <- "pA_gexreg_type"
		gene_groups <- NULL #compare all together
		ana_summary_str <- paste(study_name,expe_name,gene_grp_name)
	}}
	
	#add gene_id if no but has gene symbol
	if(!("gene_id" %in% names(genetb_data))  ){
		if( "gene_symbol" %in% names(genetb_data) ){
			genetb_data$gene_id <- geneidmapping(genetb_data$gene_symbol, inspe=organism, outspe=organism)
		}else if( "Gene_sym" %in% names(genetb_data) ){
			genetb_data$gene_id <- geneidmapping(genetb_data$Gene_sym, inspe=organism, outspe=organism)
		}
	}
	
	if("rm_redu_gid_sort_var" %in% ls()){
		if(!is.null(rm_redu_gid_sort_var)){
			print (paste(c("rm_redu_gid_sort_var = ", rm_redu_gid_sort_var, "; sort decreasing=", if_sort_decr, "; sum_sort_var_meth=",sum_sort_var_meth), collapse=" "))
			if(!sum_sort_var_meth %in% c("NA","na") ){
				genetb_data <- genetb_data[order( apply(abs(genetb_data[rm_redu_gid_sort_var]),1,sum_sort_var_meth,na.rm=T) , decreasing=if_sort_decr), ]
			}else{
				genetb_data <- genetb_data[order( genetb_data[,rm_redu_gid_sort_var], decreasing=if_sort_decr), ]
			}
			
			genetb_data<- genetb_data[!duplicated(genetb_data$gene_id) & !is.na(genetb_data$gene_id),]
			print (paste( "row# after remove duplicated gene IDs: ", nrow(genetb_data) ))
			print(table(genetb_data[gene_grp_name]))
		}
	}
	
	if("only_top_reg_num" %in% ls()){
		if(!is.null(only_top_reg_num)){
			print(paste("only_top_reg_num",only_top_reg_num,"; sort by",sel_top_sort_var))
			genetb_data<-genetb_data[!is.na(genetb_data[,sel_top_sort_var]), ]
			for( val_type in c("low","high") ){
				if(val_type %in% names(top2reg_names) ){
					genetb_data<-genetb_data[order(genetb_data[,sel_top_sort_var], decreasing=(val_type=="high") ), ]
					old_type_rows <- (1:nrow(genetb_data))[ genetb_data[,gene_grp_name]==top2reg_names[[val_type]] ]
					old_num <- length(old_type_rows)
					if(old_num>only_top_reg_num){
						genetb_data[old_type_rows[(only_top_reg_num+1):old_num], gene_grp_name] <- top2reg_names[["mid"]]
					}else if(old_num<only_top_reg_num & !if_only_remove_old){
						genetb_data[1 : only_top_reg_num , gene_grp_name] <- top2reg_names[[val_type]]
					}
				}
			}
			print(table(genetb_data[,gene_grp_name]))
		}
	}
	
	#(generate randomized pseudo-gene set with different gene number and expression level, calculate distribution of ratio and P value for randomized datasets as refeerence )
	if(ifgen_randomized_data){
		genetb_data$gex_mean<- rowMeans(genetb_data[RanGex_vars],na.rm=T)
		genetb_data<-genetb_data[order(genetb_data$gex_mean),]
		genetb_data$gex_bin<-sep_val2block(1:nrow(genetb_data), RaGex_binnum); table(genetb_data$gex_bin)
	}
	return(genetb_data)
}
