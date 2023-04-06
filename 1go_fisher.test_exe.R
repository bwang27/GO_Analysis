
#02/14/2012 add function to only select fixed number of top-regulated genes (setting only_top_reg_num)
#04/10/2012, add function to do randomization, control gene set size and gene expression level (set ifgen_randomized_data<-T; RanGex_vars )
#05/14/2012 this script is used to accept command line parameters, 
#08/20/2012 add function (mode=update_Qval) to  do randomization (control gene expression) after combine GO results of different experiment (combine mode)
#04/08/2015 changed the code to deal with msigdb data as well

string2list<-function(s){
	l<-list()
	for(v1 in unlist(strsplit(s,";"))){ #eg: "REGU:UP DN;XXXX:XX XXX
		v2<-unlist(strsplit(v1, ":")) #PA:P0 P1 P2 P3
		l[[v2[1]]] <- unlist(strsplit(v2[2], " "))
	}
	return(l)
}


##Parameter
args<- commandArgs(TRUE)
args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}
args<-NULL

if(length(grep("liw8|wl314/",getwd()))==1){
	root_dir="/HPCTMP_NOBKUP/wl314/"
}else if(!is.na(args_v["root_dir"])){
	root_dir=args_v["root_dir"]
}else{
	root_dir="/Home/cocochen/"
}
#setwd(paste(root_dir,"/analyze/go/A_enrich",sep=""))
project_root="/home/Research_Archive/ProcessedData/RNAseq/202302_HEK293_strains_RNAseq"
#setwd(paste(project_root,"/analyze/go/A_enrich",sep=""))
setwd("/home/Research_Archive/ProcessedData/RNAseq/202302_HEK293_strains_RNAseq")

go_dir<- paste("/drive2/wli/analyze/Pipeline/GO_analysis/build_go2gid/2go2gene/",sep="")  #

source("/drive2/wli/analyze/Pipeline/GO_analysis/A_enrich/1fisher.inc.R")
source("/drive2/wli/analyze/Pipeline/GO_analysis/A_enrich/2simplify_go_term_fun.R")
source("/drive2/wli/analyze/Pipeline/SharedCodes/Rfunc.inc.R")
go_cats <- c("bp","cc","mf")
go_cats2 <- c("biological_process","cellular_component","molecular_function"); names(go_cats2)<-go_cats
  
organism <- "hs" ###hs  mm rn
out_pval_cut <- 1.1 #in the raw output file: require p_min < out_pval_cut
low_pval_cut <- 0.01 #in the summary out file: require SigGrp_P<low_pval_cut
low_Qval_cut <- 0.05 #in the summary out file: require q value <low_Qval_cut  if  ifgen_randomized_data=TRUE
high_pval_cut <- 0.1 #in the summary out file: require all othGrp_Ps>high_pval_cut -Inf
min_gene_num <- 8 #minimal gene number required in an significant GO term
SigGrp_Ps<-NULL #NULL; c("P_nc","P_co") for neighbor analysis
go_childnum_cut <- 500 #when simplify the results, only keep go terms with child number < this cutoff
discard_grps<-c("0","nc_nc","NC","nc")




#default parameters:
gene_grp_name<-NULL; gene_grp_names=NULL
comb_grp_l=NULL; gene_groups=NULL
filter_names=NULL; filter_val_list=NULL; filter_unSelVal_list=NULL
auto_setup=NULL
rm_redu_gid_sort_var=NULL; if_sort_decr=F; sum_sort_var_meth="max"
if_only_remove_old=FALSE; only_top_reg_num=NULL; top2reg_names=NULL; sel_top_sort_var=NULL
##parameter for randomization
ifgen_randomized_data=FALSE; RanGex_vars=NULL; Ra_times<-500; RaGex_binnum<-10; Ra_pmin_cut<-0.05 #only do randomization if P_min<Ra_pmin_cut
gene_id_var<-'gene_id'
run_mode="fisher"
go_f_ext<-""
GOnum4genelist<-30 #extract GO number for gene list
gene_grp_var_name<-"";
gene_grp_name_all<-NULL

##parameter for combine
out_prefix<-""
comb_2P <- NULL  #NULL, c("P_UP","P_DN"), c("P_up","P_dn"), c("P_le","P_mo"); combine two p values for each sample using the more significant one, convert the second P to negative. 
discard_p <- c("P_0","P_co","P_nc_nc","P_NC")  #NULL, P_nc, P_0,  c("P_0","P_co") 
##parameters for GeneList
unsel_gene_group <- "NC" # "NC","nc", "unkn", "0" gene group must NOT be all == unsel_gene_group; set "" to disable this
discard_headers <- c("pAset","rbp","DEV.NEG","DEV.POS","BURGE","gene_symbol.1","pA_id.1","if_C_or_A","pA_type1","pAu_regu_type","test_rpm")
discard_header_pattern=NULL
change_headers_fr=NULL

##########################samples


if(length(args_v)>0){
	if(!is.na(args_v["go_dir"])){ go_dir<-args_v["go_dir"] }
	if(!is.na(args_v["go_cats"])){ go_cats<- unlist(strsplit(args_v["go_cats"]," ")); go_cats2=go_cats  }
	gset_version<- ifelse(is.na(args_v["gset_version"]), "v5.0", args_v["gset_version"])
	
	
	organism<- ifelse(is.na(args_v["organism"]), "mm", args_v["organism"])
	if(!is.na(args_v["run_mode"])){ run_mode<-args_v["run_mode"] } #fisher; combine; GeneList (extract gene for top GOs)
	if(!is.na(args_v["study_name"])){ study_name<-args_v["study_name"] }
#	out_dir=ifelse(is.na(args_v["out_dir"]), paste0("1GOFisher/",study_name,"/"), args_v["out_dir"])
	out_dir=ifelse(is.na(args_v["out_dir"]), paste0(study_name,"/","1GOFisher/"), args_v["out_dir"])
	if(!is.na(args_v["expe_name"])){ expe_name<-args_v["expe_name"] }
	if(!is.na(args_v["genetb_file"])){ genetb_file<-args_v["genetb_file"] }
	if(!is.na(args_v["gene_grp_name"])){ gene_grp_name<-args_v["gene_grp_name"] }
	if(!is.na(args_v["gene_grp_names"])){ gene_grp_names<- unlist(strsplit(args_v["gene_grp_names"]," "))  }
	if(!is.na(args_v["gene_groups"])){ gene_groups<- unlist(strsplit(args_v["gene_groups"]," "))  }
	if(!is.na(args_v["comb_grp_l"])){ comb_grp_l<- string2list(args_v["comb_grp_l"]) }
	if(!is.na(args_v["ifgen_randomized_data"])){ if(args_v["ifgen_randomized_data"] %in% c('TRUE','T','1')){ifgen_randomized_data<-TRUE} }
	if(!is.na(args_v["RanGex_vars"])){ RanGex_vars<- unlist(strsplit(args_v["RanGex_vars"]," "))  }
	if(!is.na(args_v["rm_redu_gid_sort_var"])){ rm_redu_gid_sort_var<- unlist(strsplit(args_v["rm_redu_gid_sort_var"]," ")) }
	if(!is.na(args_v["if_sort_decr"])){ if_sort_decr<- args_v["if_sort_decr"] %in% c("T","1","True") }
	if(!is.na(args_v["sum_sort_var_meth"])){ sum_sort_var_meth<-args_v["sum_sort_var_meth"] }
	if(!is.na(args_v["filter_names"])){ filter_names<- unlist(strsplit(args_v["filter_names"]," "))  }
	if(!is.na(args_v["filter_val_list"])){ filter_val_list<- string2list(args_v["filter_val_list"]) }
	if(!is.na(args_v["filter_unSelVal_list"])){ filter_unSelVal_list<- string2list(args_v["filter_unSelVal_list"]) }
	if(!is.na(args_v["gene_id_var"])){ gene_id_var<-args_v["gene_id_var"] } #eg: Gene_id
	if(!is.na(args_v["Ra_times"])){ Ra_times<-as.numeric(args_v["Ra_times"]) } 
	if(!is.na(args_v["RaGex_binnum"])){ RaGex_binnum<-as.numeric(args_v["RaGex_binnum"]) } 
	if(!is.na(args_v["Ra_pmin_cut"])){ Ra_pmin_cut<-as.numeric(args_v["Ra_pmin_cut"]) } 
	if(!is.na(args_v["low_pval_cut"])){ low_pval_cut<-as.numeric(args_v["low_pval_cut"]) } 
	if(!is.na(args_v["low_Qval_cut"])){ low_Qval_cut<-as.numeric(args_v["low_Qval_cut"]) } 
	if(!is.na(args_v["only_top_reg_num"])){ only_top_reg_num<-as.numeric(args_v["only_top_reg_num"]) } 
	if(!is.na(args_v["top2reg_names"])){ top2reg_names<- string2list(args_v["top2reg_names"]) }
	if(!is.na(args_v["sel_top_sort_var"])){ sel_top_sort_var<-args_v["sel_top_sort_var"] }
	if(!is.na(args_v["if_only_remove_old"])){ if_only_remove_old<-TRUE }
	
	#for combine mode
	if(!is.na(args_v["exp_type"])){ exp_type<-args_v["exp_type"] } 
	if(!is.na(args_v["sample_names"])){ sample_names<- unlist(strsplit(args_v["sample_names"]," "))  }
	if(!is.na(args_v["comb_2P"])){ comb_2P<- unlist(strsplit(args_v["comb_2P"]," "))  }
	if_comb_2P_signs=T
	if(!is.na(args_v["if_comb_2P_signs"])){ if_comb_2P_signs<- args_v["if_comb_2P_signs"] %in% c("T","1","True")  }
	comb_2P_SS_method=ifelse(is.na(args_v["comb_2P_SS_method"]), "moreSig", args_v["comb_2P_SS_method"]) #moreSig or sum
	if(!is.na(args_v["go_f_ext"])){ go_f_ext<-args_v["go_f_ext"] } 
	if(!is.na(args_v["discard_p"])){ discard_p<- unlist(strsplit(args_v["discard_p"]," "))  }
	if(!is.na(args_v["out_prefix"])){ out_prefix<-args_v["out_prefix"] } 
	combine_meth=ifelse(is.na(args_v["combine_meth"]), "row_min_p", args_v["combine_meth"]) #can set to topNinEachP
	eachP_GO_num=ifelse(is.na(args_v["eachP_GO_num"]), 10, as.numeric(args_v["eachP_GO_num"]))
	max_row_num=ifelse(is.na(args_v["max_row_num"]), 60, as.numeric(args_v["max_row_num"]))
	comb_sample_l=NULL
	if(!is.na(args_v["comb_sample_l"])){ comb_sample_l<- string2list(args_v["comb_sample_l"]) }
	comb_sample_SS_method=ifelse(is.na(args_v["comb_sample_SS_method"]), "average", args_v["comb_sample_SS_method"])
	if_cluster_cb_GO_tb_rows=F; if_cluster_cb_GO_tb_cols=F;
	if(!is.na(args_v["if_cluster_cb_GO_tb_rows"])){ if_cluster_cb_GO_tb_rows<- args_v["if_cluster_cb_GO_tb_rows"] %in% c("T","1","True")  }
	if(!is.na(args_v["if_cluster_cb_GO_tb_cols"])){ if_cluster_cb_GO_tb_cols<- args_v["if_cluster_cb_GO_tb_cols"] %in% c("T","1","True")  }

	#for run_mode=="GeneList"
	if(!is.na(args_v["GOnum4genelist"])){ GOnum4genelist<-as.numeric(args_v["GOnum4genelist"]) } 
	if(!is.na(args_v["unsel_gene_group"])){ unsel_gene_group<-args_v["unsel_gene_group"] } 
	if(!is.na(args_v["gene_grp_name_all"])){ gene_grp_name_all<- unlist(strsplit(args_v["gene_grp_name_all"]," "))  } #or set gene_grp_var_name and sample_names
	if(!is.na(args_v["gene_grp_var_name"])){ gene_grp_var_name<- args_v["gene_grp_var_name"]  }
	if(!is.na(args_v["discard_headers"])){ discard_headers<- unlist(strsplit(args_v["discard_headers"]," "))  }
	if(!is.na(args_v["discard_header_pattern"])){ discard_header_pattern<- args_v["discard_header_pattern"]  }
	if(!is.na(args_v["change_headers_fr"])){ change_headers_fr<- unlist(strsplit(args_v["change_headers_fr"]," "))  }
	if(!is.na(args_v["change_headers_to"])){ change_headers_to<- unlist(strsplit(args_v["change_headers_to"]," "))  }
	
}else{quit()}

system(paste("mkdir -p ",out_dir,"/",sep=""))
out_img_f=paste(out_dir,"/run_mode.",run_mode,".",gene_grp_name,".img",sep="")
save.image(out_img_f)

if(run_mode=="fisher"){
	genetb_data <- read.table(genetb_file,sep="\t",header=T,comment.char="",quote="",stringsAsFactors=F)
	
	##################pre_process_data
	genetb_data <- pre_process_data(organism,genetb_data, gene_grp_name, gene_grp_names=gene_grp_names, comb_grp_l=comb_grp_l, gene_groups=gene_groups, 
		filter_names=filter_names, filter_val_list=filter_val_list,
		auto_setup=auto_setup, rm_redu_gid_sort_var=rm_redu_gid_sort_var, if_sort_decr=if_sort_decr, 
		if_only_remove_old=if_only_remove_old, only_top_reg_num=only_top_reg_num, top2reg_names=top2reg_names, sel_top_sort_var=sel_top_sort_var,
		ifgen_randomized_data=ifgen_randomized_data, RanGex_vars=RanGex_vars, RaGex_binnum=RaGex_binnum)
	
	#########################################################################
	# note:
	# genetb_data should has two columes: gene_id, gene_grp_name ('no', 'up', 'dn', 'spli_changed' ... values for gene_grp_name)
	# should set gene_grp_name, study_name, expe_name,   values before run function do_all_go_ana
	##optional variables:
	#rm_redu_gid_sort_var #use abs(rm_redu_gid_sort_var) to sort genetb_data (decreasing=if_sort_decr), then remove duplicated gene_id
	#comb_grp_l #combine values in one column (gene_grp_name); eg: comb_grp_l <- list(regu=c("in","ex"))
	#use filter_names and filter_val_list to filter rows in input table. eg. use filter_names=c("pA_type"), filter_val_list=list(c("3F","3M")) to only study 3F and 3M polyA sites
	#use gene_grp_names to combine combine several column as gene group 
	
	###############1,  do fisher's exact test, output raw data################
	
	gene_groups_str <- paste(gene_groups,collapse="_VS_"); gene_groups_str
	ana_summary_str <- paste(study_name,expe_name,gene_grp_name, substr(filterstr,1,50)); print(ana_summary_str)
	outfilebase <- paste(out_dir,sep="")
	out_summary_file <- paste(outfilebase,expe_name,".",gene_grp_name,".",gene_groups_str,".GO_FisherT.tbl",sep=""); out_summary_file
	mkdir_if_not_exist(out_summary_file)
	do_all_go_ana(genetb_data, organism, ana_summary_str, outfilebase)
	
	###########################################################################
	####2, simplify raw results (remove redundance)############################
	if(all(go_cats2 %in% c("bp","cc","mf")) ){
		SigGrp_Ps<-NULL ###
		
		simpl_out_f <- paste(outfilebase,expe_name,".",gene_grp_name,".",gene_groups_str,".GO_FisherT.sel_sigGO.tbl",sep="")
		print (simpl_out_f)
		goIn_D <- read.table(out_summary_file, header=T, sep="\t", comment.char="#", quote="", strip.white=T, stringsAsFactors=F) #
		headerRowIds <- (1:nrow(goIn_D))[goIn_D[,1]==names(goIn_D)[1]]
		goIn_D["go_cat"] <- c(rep("bp",headerRowIds[1]-1), rep("cc",headerRowIds[2]-headerRowIds[1]), rep("mf",nrow(goIn_D)-headerRowIds[2]+1))
		goIn_D <- goIn_D[goIn_D[,1]!=names(goIn_D)[1],]
		
		numericVars <- grep("^P_|^Num_", names(goIn_D), value=T)
		temp <- sapply(numericVars, function(v){goIn_D[[v]]<<- as.numeric(goIn_D[[v]])})
		if(is.null(SigGrp_Ps) | length(SigGrp_Ps)<1){
			SigGrp_Ps <- setdiff(grep("^P_", numericVars, value=T), c("P_Min","P_magnitude_diff","P_nc","P_NC","P_0","P_nc_nc"))
		}
		out_tb <- NULL
		for(SigGrp_P in SigGrp_Ps){
			othGrp_Ps <- setdiff(grep("^P_", numericVars, value=T), c("P_Min","P_magnitude_diff",SigGrp_P))
			geneNum_name<- sub("P_","Num_",SigGrp_P)
			goIn_D1 <- goIn_D[ goIn_D[,geneNum_name]>=min_gene_num & goIn_D$Num_Child<go_childnum_cut & goIn_D[,SigGrp_P]<low_pval_cut ,] #
			if(length(othGrp_Ps)>0){
				goIn_D1 <- goIn_D[ rowSums(goIn_D[othGrp_Ps]<high_pval_cut)==0 ,]
			}
			#add filter column (for selected GO terms)
			goIn_D1$filter<-as.character(unlist(sapply(go_cats, function(go_cat){
				item_num <- sum(goIn_D1$go_cat==go_cat)
				if(item_num==0){
					NULL
				}else if(item_num<=5){
					rep("1",item_num)
				}else{
					tmp_D <- goIn_D1[goIn_D1$go_cat==go_cat,]
					shorten_go_terms.fun(go_cat, tmp_D,  SigGrp_P, pval_cut=low_pval_cut, spe=organism)
				}
			})))
			goIn_D1$sig_grp <- rep(sub("P_","",SigGrp_P),nrow(goIn_D1))
			goIn_D1<- goIn_D1[goIn_D1$filter==1, ]
			goIn_D1<- goIn_D1[order(goIn_D1$go_cat, goIn_D1[[SigGrp_P]]),]
			out_tb<- rbind(out_tb,goIn_D1)
		}
		#output simplified out file
		out_tb<-out_tb[ c( setdiff(names(out_tb),c("P_Min","P_magnitude_diff","P_nc","P_0","filter","GO_term")), "GO_term") ]
		write.table(out_tb, file=simpl_out_f, col.names=T, row.names=F, sep="\t", quote=F)
	}

}


if(run_mode=="combine"){
	########## combine go term analysis for multi analysis######
	# ##mm_3pS.HURKD_3pS.miR133KD_3pS.dsRBD_OE_3pS
	# #Gene expression #1:6 c(1,7,8)
	# sample_names <- paste("Gene_ReguType",test_samples[1:6],ref_samples[1:6],sep="_")
	# exp_type <- "Gex" ###
	# go_f_names <- paste(exp_type,'.',sample_names,".DN_VS_NC_VS_UP.GO_FisherT.tbl",sep="")
	# 
	# ##eg: Carol
	# #APA
	#	sample_names <- paste(test_samples,ref_samples,sep="_")
	#	exp_type <- "pAu.APA_Type_" ###
	#	go_f_ext<- ".REGU_VS_NC.GO_FisherT.tbl"
	#	#gene expression
	#	exp_type <- "Gex.Gene_ReguType_"
	#	go_f_ext<- ".DN_VS_NC_VS_UP.GO_FisherT.tbl"
	go_f_names <- paste(exp_type,sample_names,go_f_ext,sep="")
	
	##########3, RUN:   ################################
	comtb_comm_headers <-c("GO_term","Num_Child","Num_Parent","Num_Genes")
	go_fs <- ifelse( go_f_names %in% grep("GOFisher",go_f_names,value=T) , go_f_names, paste(out_dir, go_f_names, sep="") ); go_fs
	names(go_fs) <- sample_names
	out_root <- paste(out_dir,"/",out_prefix, exp_type, sep=""); print(out_root)
	
	d <- list()
	all_go_l <- list()
	for(sample_name in sample_names){ #read raw GO results data
		go_f <- go_fs[sample_name]
		print(go_f)
		goIn_D <-read.table(go_f, header=T, sep="\t", comment.char="#", quote="", strip.white=T, stringsAsFactors=F) #
		headerRowIds <- (1:nrow(goIn_D))[goIn_D[,1]==names(goIn_D)[1]]
		print( paste(c("headerRowIds=",headerRowIds),collapse=" ") )
		#goIn_D["go_cat"] <- c(rep("bp",headerRowIds[1]-1), rep("cc",headerRowIds[2]-headerRowIds[1]), rep("mf",nrow(goIn_D)-headerRowIds[2]+1))
		go_cats_v=rep(go_cats[1],nrow(goIn_D))
		if(length(headerRowIds)>0){
			for(rowId_i in 1:length(headerRowIds)){
				headerRowId=headerRowIds[rowId_i]
				go_cats_v[headerRowId:length(go_cats_v)]=go_cats[rowId_i+1]
			}
		}
		goIn_D["go_cat"]=go_cats_v
		comtb_comm_headers=intersect(comtb_comm_headers,names(goIn_D))
		goIn_D <- goIn_D[goIn_D[,1]!=names(goIn_D)[1],]
		numericVars <- grep("^P_|^Num_", names(goIn_D), value=T)
		temp <- sapply(numericVars, function(v){goIn_D[[v]]<<- as.numeric(goIn_D[[v]])})
		SigGrp_Ps <- setdiff(grep("^P_", numericVars, value=T), c("P_Min","P_magnitude_diff",discard_p))
		for(go_cat in go_cats){
			all_go_l[[go_cat]] <- unique(c(all_go_l[[go_cat]], goIn_D$uniq_go_ids[goIn_D$go_cat==go_cat]))
		}
		d[[sample_name]] <- goIn_D
	}
	
	####combine p values , 
	for(go_cat in go_cats){
		cb_tb <- data.frame(GO_ID=all_go_l[[go_cat]])
		for(sample_name in sample_names){
			print(paste(go_cat,sample_name))
			tmptb <- d[[sample_name]]
			tmptb <- tmptb[tmptb$go_cat==go_cat,]
			
			SigGrp_Ps <- setdiff(grep("^P_", names(tmptb), value=T), c("P_Min","P_magnitude_diff",discard_p))
			comb_p_name <- SigGrp_Ps
			tmptb[SigGrp_Ps][tmptb[SigGrp_Ps]<10^(-99)] <- 10^(-99) #pseudo p value when P=0
			if(!is.null(comb_2P) ){ #combine two p values
				if(all(comb_2P %in% SigGrp_Ps)){
					tmptb[[comb_2P[2]]] <- -tmptb[[comb_2P[2]]]
					if(if_comb_2P_signs){
						comb_p_name <- paste(comb_2P,collapse="_")
						if(comb_2P_SS_method=="moreSig"){
							tmptb[comb_p_name] <- apply(abs(tmptb[comb_2P]),1, min)
							ifminus <- abs(tmptb[[comb_2P[2]]]) < tmptb[[comb_2P[1]]]
							tmptb[ifminus,comb_p_name] <- -tmptb[ifminus,comb_p_name]					
						}else if(comb_2P_SS_method=="sum"){
							ss_tb=-log10(abs(tmptb[,comb_2P])) * sign(tmptb[,comb_2P])
							ss_sum=rowSums(ss_tb,na.rm=T)
							tmptb[,comb_p_name]=(10^( -abs(ss_sum)) ) * ifelse(ss_sum==0,1, sign(ss_sum))
						}
					}
				}
			}
			
			gnum_names<-c("Num_Exp", sub("^P_","Num_",SigGrp_Ps))
			Q_names<-sub("^P_","Q_",SigGrp_Ps)
			Q5q_names<-sub("^P_","Q5q_",SigGrp_Ps)
			update_names<-c(gnum_names,comb_p_name)
			if(ifgen_randomized_data){
				update_names<-c(update_names,Q_names)
			}
			print(update_names)
			if(sample_name==sample_names[1]){
				cb_tb[c(comtb_comm_headers, paste(update_names,sample_name,sep="."))] <- tmptb[match(cb_tb$GO_ID,tmptb$uniq_go_ids), c(comtb_comm_headers,update_names)]
			}else{
				cb_tb[ paste(update_names,sample_name,sep=".")] <- tmptb[match(cb_tb$GO_ID,tmptb$uniq_go_ids), update_names]
			}
		}
		
		print(names(cb_tb))
		#pnames_tmp=grep("P_",names(cb_tb),value=T)
		#ifsel=rowSums(cb_tb[,pnames_tmp]==0,na.rm=T)>0
		#print(cb_tb[ifsel,c("GO_ID",pnames_tmp)])
		
		if(!is.null(comb_sample_l) ){ #combine P from different sample (eg. replicates)
			for(gene_group in gene_groups){
				for(sample1 in names(comb_sample_l)){
					cbP_names=paste("P_",gene_group,".",comb_sample_l[[sample1]],sep="")
					new_P_name=paste("P_",gene_group,".",sample1,sep="")
					print (paste(c("combine",comb_sample_SS_method,cbP_names, "--->", new_P_name), collapse=" "))
					ss_tb=-log10(abs(cb_tb[,cbP_names])) * sign(cb_tb[,cbP_names])
					if(comb_sample_SS_method=="average"){
						ss_mean=rowMeans(ss_tb,na.rm=T)
						cb_tb[,new_P_name]=(10^( -abs(ss_mean)) ) * ifelse(ss_mean==0,1, sign(ss_mean))
					}else if(comb_sample_SS_method=="min"){ #select least significant SS
						cb_tb[,new_P_name]=apply(cb_tb[,cbP_names],1,function(v){
							v[which.max(abs(v))]
						})
					}
					remove_P_names=setdiff(cbP_names,new_P_name)
					cb_tb=cb_tb[setdiff(names(cb_tb), remove_P_names)]
				}
			}
		}
		
		#save.image("tmp.img")
		all_pvalnames <- grep("^P_",names(cb_tb), value=T); print(all_pvalnames)
		all_qvalnames <- grep("^Q_",names(cb_tb), value=T)
		all_q5qvalnames <- grep("^Q5q_",names(cb_tb), value=T)
		all_Expnum_names <- grep("^Num_Exp",names(cb_tb), value=T)
		all_reguNum_names <- setdiff(grep("^Num_",names(cb_tb), value=T), c(comtb_comm_headers,all_Expnum_names)); all_reguNum_names
		
		table(cb_tb[all_pvalnames]==0)
		
		cb_tb[all_pvalnames][is.na(cb_tb[all_pvalnames])] <- 1
		if(ifgen_randomized_data){
			cb_tb[c(all_qvalnames,all_q5qvalnames)]<- as.numeric(unlist(cb_tb[c(all_qvalnames,all_q5qvalnames)]))
		}

		cb_tb$min_p <- apply(abs(cb_tb[all_pvalnames]),1, min, na.rm=T)
		cb_tb1 <- cb_tb[cb_tb$min_p<low_pval_cut & !is.na(cb_tb$min_p), ]
		
		if(length(all_qvalnames)>0){
			min_q<-apply(cb_tb1[,all_qvalnames],1,min, na.rm=T)
			cb_tb1 <- cb_tb1[min_q<low_Qval_cut & !is.na(min_q) & abs(min_q)!=Inf, ]
		}
		if("Num_Child" %in% names(cb_tb1)){
			cb_tb1<-cb_tb1[ !is.na(cb_tb1$Num_Child) & cb_tb1$Num_Child<go_childnum_cut, ]
		}
		if(length(all_pvalnames)==length(all_reguNum_names)){
			cb_tb1<-cb_tb1[ rowSums(cb_tb1[,all_reguNum_names]>=min_gene_num & cb_tb1[,all_pvalnames]<low_pval_cut)>0,  ]
		}else{
			cb_tb1<-cb_tb1[ apply(cb_tb1[all_reguNum_names],1,max)>=min_gene_num,  ]
		}
		
		cb_tb1 <- cb_tb1[order(cb_tb1$min_p),]
		
		if(combine_meth=="topNinEachP"){#top GO in each p value
			cb_tb1$filter="0"
			cb_tb2=cb_tb1
			cb_tb2$ori_rowid=1:nrow(cb_tb2)
			all_selRowIDs=NULL
			for(pval1 in all_pvalnames){
				cb_tb2[,pval1]=abs(cb_tb2[,pval1])
				cb_tb2=cb_tb2[order(cb_tb2[,pval1]),]
				if(nrow(cb_tb2)>eachP_GO_num*10){
					cb_tb3=cb_tb2[1:(eachP_GO_num*10),]
				}else{
					cb_tb3=cb_tb2
				}
				if(go_cat %in% c("bp","cc","mf")){
					cb_tb3$filter=shorten_go_terms.fun(go_cat, cb_tb3,  pval1, go_tern_name="GO_ID",  pval_cut=low_pval_cut, spe=organism)
				}else{
					cb_tb3$filter=1
				}
				selRowIDs=cb_tb3$ori_rowid[!is.na(cb_tb3$filter) & cb_tb3$filter=="1"]
				if(length(selRowIDs)>eachP_GO_num){selRowIDs=selRowIDs[1:eachP_GO_num]}
				all_selRowIDs=c(all_selRowIDs,selRowIDs)
				print( pval1 )
				print(selRowIDs)
			}
			all_selRowIDs=all_selRowIDs[!duplicated(all_selRowIDs)]
			cb_tb1=cb_tb1[all_selRowIDs,]
		}else{
				if(go_cat %in% c("bp","cc","mf")){
					cb_tb1$filter <- shorten_go_terms.fun(go_cat, cb_tb1,  "min_p", go_tern_name="GO_ID",  pval_cut=low_pval_cut, spe=organism)
				}else{
					cb_tb1$filter=1
				}
			cb_tb1<- cb_tb1[cb_tb1$filter=="1", ] #remove redundant GO
		}
		if(nrow(cb_tb1)>max_row_num){ #filter by row number
			cb_tb1 <- cb_tb1[1:max_row_num, ]
		}
		
		cb_tb1[c(all_pvalnames,all_q5qvalnames)] <- round( -log10(abs(cb_tb1[c(all_pvalnames,all_q5qvalnames)])) * sign(cb_tb1[c(all_pvalnames,all_q5qvalnames)]), 1)
		cb_tb1[all_pvalnames][is.na(cb_tb1[all_pvalnames])] <- 0
		if(if_cluster_cb_GO_tb_rows){
				tb=cb_tb1[,all_pvalnames]
				hc_str=try( hclust(dist(tb), "ave") )
				if(length(grep("Error",hc_str))==0){
					hc <- hclust(dist(tb), "ave")
					cb_tb1=cb_tb1[hc$order, ]
				}
		}
		if(if_cluster_cb_GO_tb_cols){
				tb=cb_tb1[,all_pvalnames]
				hc_str=try( hclust(dist(t(tb)), "ave") )
				if(length(grep("Error",hc_str))==0){
					hc <- hclust(dist(tb), "ave")
					all_pvalnames=all_pvalnames[hc$order]
				}
		}
		cb_tb1<-cb_tb1[ c("GO_ID",comtb_comm_headers[-1],all_Expnum_names,all_reguNum_names,all_pvalnames,all_qvalnames,all_q5qvalnames,"GO_term") ]
		
		out_f <- paste(out_root,go_cat, "txt", sep=".")
		write.table(cb_tb1, file=out_f, col.names=T, row.names=F, sep="\t", quote=F)
	}

}

if(run_mode=="update_Qval"){
	sel_GOrows <- 1:GOnum4genelist
	genetb_data <- read.table(genetb_file,sep="\t",header=T,comment.char="#",quote="",stringsAsFactors=F)
	genetb_data <- pre_process_data(organism,genetb_data, gene_grp_name, gene_grp_names=gene_grp_names, comb_grp_l=comb_grp_l, gene_groups=gene_groups, filter_names=filter_names, filter_val_list=filter_val_list,
		auto_setup=auto_setup, rm_redu_gid_sort_var=rm_redu_gid_sort_var, if_sort_decr=if_sort_decr, if_only_remove_old=if_only_remove_old, only_top_reg_num=only_top_reg_num, top2reg_names=top2reg_names,
		ifgen_randomized_data=ifgen_randomized_data, RanGex_vars=RanGex_vars, RaGex_binnum=RaGex_binnum)
	gene_grp_types <- sort(unique(genetb_data[[gene_grp_name]]))
	gene_grp_types <- setdiff(gene_grp_types,c(discard_grps,unsel_gene_group))
	
	out_root <- paste(out_dir,"/",out_prefix, exp_type, sep=""); out_root
	for(go_cat in go_cats){ #for each top go term, find regulated gene list 
		#top GO data
		topGO_f <- paste(out_root,go_cat, "txt", sep=".")
		print(paste("open",topGO_f))
		topGo_d <- read.table(topGO_f, header=T, sep="\t", stringsAsFactors=F, comment.char="", quote="")
		names(topGo_d)<-change_values(names(topGo_d), "uniq_go_ids","GO_ID")
		sel_GOrows2<-intersect(sel_GOrows, 1:nrow(topGo_d))
		
		#GO annotation file
		go2geneid_file <- go_cat2_go2id_f(go_cat, organism)
		GO_tbl <- read.table(file=go2geneid_file, sep="\t")
		colnames(GO_tbl) <- c("GO", "Gene") ; nrow(GO_tbl)
		print(genetb_data[1,])
		print( c("total rows of genetb_data:", nrow(genetb_data)) )
		topGo_d=do_randomSampling(topGo_d, genetb_data, gene_grp_name,gene_grp_types, sel_GOrows2, GO_tbl, sample_name=sample_names)
		
		#write topGo_d, replace original one
		write.table(topGo_d,file=topGO_f, col.names=T, row.names=F, sep="\t", quote=F)
	}	
}


if(run_mode=="GeneList"){
	################4, for selected GOs(usually top ones), extract coresponding genes expression or other information########
	#######need to use output of step 3, and need genetb_data in step1, remove unnesessary rows first
	sel_GOrows <- 1:GOnum4genelist
	
	
	
	##mm_3pS.HURKD_3pS.miR133KD_3pS.dsRBD_OE_3pS
	##gex
	#gene_grp_name_all <- paste("Gene_ReguType",test_samples,ref_samples,sep="_"); gene_grp_name_all
	if(is.null(gene_grp_name_all)){
		gene_grp_name_all<-paste(gene_grp_var_name,sample_names,sep="_")
	}
	
	###############RUN
	genetb_data <- read.table(genetb_file,sep="\t",header=T,comment.char="#",quote="",stringsAsFactors=F)
	##################pre_process_data
	genetb_data <- pre_process_data(organism,genetb_data, gene_grp_name, gene_grp_names=gene_grp_names, comb_grp_l=comb_grp_l, gene_groups=gene_groups, filter_names=filter_names, filter_val_list=filter_val_list,
		auto_setup=auto_setup, rm_redu_gid_sort_var=rm_redu_gid_sort_var, if_sort_decr=if_sort_decr, if_only_remove_old=if_only_remove_old, only_top_reg_num=only_top_reg_num, top2reg_names=top2reg_names,
		ifgen_randomized_data=ifgen_randomized_data, RanGex_vars=RanGex_vars, RaGex_binnum=RaGex_binnum)

	out_root <- paste(out_dir,"/",out_prefix, exp_type, sep=""); out_root
	#!!!need genetb_data, make sure genetb_data is new opened, (certain steps (eg. step 1) remove rows in genetb_data)
	
	#add gene description:
	if( !("gene_desc" %in% names(genetb_data)) & !("Gene_desc" %in% names(genetb_data)) ){
		genetb_data["gene_desc"] <- gid2_gene_desc(genetb_data$gene_id, spe=organism)[,"gene_desc"]
	}
	
	for(go_cat in go_cats){ #for each top go term, find regulated gene list 
		topGO_f <- paste(out_root,go_cat, "txt", sep=".")
		print(paste("open",topGO_f))
		topGo_d <- read.table(topGO_f, header=T, sep="\t", stringsAsFactors=F, comment.char="", quote="")
		names(topGo_d)<-change_values(names(topGo_d), "uniq_go_ids","GO_ID")
		out_f <- paste(topGO_f,".top.gene",sep="") ###
		out_f2 <- paste(out_f,"detail",sep=".") ###
		sel_GOrows2<-intersect(sel_GOrows, 1:nrow(topGo_d))
		selGOs <- topGo_d$GO_ID[sel_GOrows2]
		selGO_descs <- topGo_d$GO_term[sel_GOrows2]; print(selGO_descs)
		names(selGO_descs) <- selGOs
		##get gene ids associated with selGOs
		go2geneid_file <- go_cat2_go2id_f(go_cat, organism)
		GO_tbl <- read.table(file=go2geneid_file, sep="\t")
		colnames(GO_tbl) <- c("GO", "Gene") ; nrow(GO_tbl)
		sel_gids <- GO_tbl$Gene[GO_tbl$GO %in% selGOs]
		
		#sel_genetb <- genetb_data[!is.na(genetb_data$gene_id) & genetb_data$gene_id %in% sel_gids, ]
		sel_genetb <- genetb_data
		if("sort_vars" %in% ls()){
			print (paste(c("sort_vars=",sort_vars),collapse=" " ))
			if(sort_vars %in% names(sel_genetb)){
				sel_genetb <- sel_genetb[order( apply(abs(sel_genetb[sort_vars]),1,max,na.rm=T), decreasing=T ), ]
			}
		}
		
		print(gene_grp_name_all)
		if(unsel_gene_group != ""){
			print (paste(c("unsel_gene_group=",unsel_gene_group),collapse=" " ))
			sel_reg_genetb <- sel_genetb[rowSums(!is.na(sel_genetb[gene_grp_name_all]) & sel_genetb[gene_grp_name_all]!=unsel_gene_group &  sel_genetb[gene_grp_name_all]!="NA" & sel_genetb[gene_grp_name_all]!="na" & sel_genetb[gene_grp_name_all]!="unkn")>0, ]
		}else{
			sel_reg_genetb=sel_genetb
		}
		
		#anotate genes with GO
		write("top regulated:", out_f2)
		i<-0
		for(selGO in selGOs){
			i=i+1;
			if_go <- as.numeric(sel_reg_genetb$gene_id %in% GO_tbl$Gene[GO_tbl$GO==selGO])
			sel_reg_genetb[selGO_descs[selGO]] <- if_go
			write("", out_f2,append=T)
			write.table(topGo_d[topGo_d$GO_ID==selGO,c("GO_ID","GO_term")], file=out_f2,col.names=F, row.names=F, sep="|", quote=F, append=T)
			go1_tb <- data.frame(GO_rank=rep(i,sum(if_go)), sel_reg_genetb[if_go==1, 1:ncol(genetb_data)])
			#names(go1_tb)[1]<-"GO_rank"
			go1_tb <- go1_tb[ setdiff(names(go1_tb),discard_headers)  ]
			if(!is.null(discard_header_pattern)){go1_tb <- go1_tb[ , -grep(discard_header_pattern,names(go1_tb))  ]}
			if(!is.null(change_headers_fr)){ names(go1_tb)<-change_values(names(go1_tb), change_headers_fr,change_headers_to)  }
			write.table(go1_tb, file=out_f2, col.names=T, row.names=F, sep="\t", quote=F, append=T)
		}
		sel_reg_genetb <- sel_reg_genetb[ setdiff(names(sel_reg_genetb),discard_headers)  ]
		if(!is.null(discard_header_pattern)){sel_reg_genetb <- sel_reg_genetb[ , -grep(discard_header_pattern,names(sel_reg_genetb))  ]}
		if(!is.null(change_headers_fr)){ names(sel_reg_genetb)<-change_values(names(sel_reg_genetb), change_headers_fr,change_headers_to)  }
		write.table(sel_reg_genetb,file=out_f, col.names=T, row.names=F, sep="\t", quote=F)
	}
	
	
}

