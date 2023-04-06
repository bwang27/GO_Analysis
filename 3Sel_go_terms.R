##select GO terms based on interest

root_dir="../../../"
source("../../Rfunc.inc.R")
go_dir<- paste("2go2gene/",sep="")  #
organism="hs"
tax_id=tax_ids[organism]
latin_name=spe_latin_names[as.character(tax_id)]

go2geneid_file=paste(go_dir,organism,".go2gene.tbl", sep="")
go2ano_file=paste(go_dir,organism,".go2gene.num.tbl", sep="")

go2gene_d <- read.table(file=go2geneid_file, sep="\t",header=T)
go_anot_data <- read.table(file=go2ano_file, sep="\t", header=T,quote="",comment.char = "")


##1, build gene set for interested GO terms
GOgroup_l=list(neurology="nervous|neuron|brain", 
		ophthalmology="eye development|\\svision|ophthalmology", 
		muscle="muscle")
all_gene_d=data.frame(gene_id=unique(go2gene_d$Gene_ID))
all_gene_d[c("gene_symbol","gene_desc")]=gid2_gene_desc(all_gene_d$gene_id, spe=organism)[c("gene_symbol","gene_desc")]
for(GOgroup_name in names(GOgroup_l)){
	search_patt=GOgroup_l[[GOgroup_name]]
	out_go2gene_f=paste("3GO_set/",organism,".",GOgroup_name,".gene.txt",sep="")
	out_goAno_f=paste("3GO_set/",organism,".",GOgroup_name,".ano.txt",sep="")
	mkdir_if_not_exist(out_go2gene_f)

	write(paste("#search pattern=",search_patt), file=out_goAno_f)
	go_anot_d2=go_anot_data[grep(search_patt, go_anot_data$GO_term, value=F),]; dim(go_anot_d2)
	go_anot_d2$GO_term #manually check

	gids=unique(go2gene_d$Gene_ID [go2gene_d$GO_ID %in% go_anot_d2$GO_ID ]); length(gids)
	write(paste("#found GO terms=",nrow(go_anot_d2)), file=out_goAno_f, append=T)
	write(paste("#found genes=", length(gids)), file=out_goAno_f, append=T)
	write.table(go_anot_d2, file=out_goAno_f, row.names=F, col.names=T, sep="\t", append=T, quote=F)

	out_go2gene_d=data.frame(gene_id=gids)
	out_go2gene_d[c("gene_symbol","gene_desc")]=gid2_gene_desc(gids, spe=organism)[c("gene_symbol","gene_desc")]
	write.table(out_go2gene_d, file=out_go2gene_f, row.names=F, col.names=T, sep="\t", quote=F)

	#update all_gene_d
	all_gene_d[,GOgroup_name]=ifelse(all_gene_d$gene_id %in% out_go2gene_d$gene_id, "Y","N")
}
#output all_gene_d
out_allgene_with_go_f=paste("3GO_set/",organism,".gene.interested_GO.txt",sep="")
write.table(all_gene_d, out_allgene_with_go_f, row.names=F, col.names=T, sep="\t", quote=F)



##2, build a table with all genes and all non-generic GOs associated
update_headers=setdiff(names(go_anot_data), "GO_ID")
go2gene_d2=go2gene_d
go2gene_d2[update_headers]=go_anot_data[match(go2gene_d2$GO_ID, go_anot_data$GO_ID),update_headers]
ifsel=go2gene_d2$Category=="biological_process" & go2gene_d2$Num_Genes<=2500
go2gene_d2=go2gene_d2[ifsel,]; nrow(go2gene_d2)


raw_gid2go_f=paste("../../../data/ncbi/gene/go/gene2go",sep="")
raw_gid2go_d=read.table(raw_gid2go_f, sep="\t", header=F,quote="",comment.char = "#")
names(raw_gid2go_d)=c("tax_id","Gene_ID","GO_ID","Evidence","Qualifier","GO_term","PubMed","Category")
raw_gid2go_d=raw_gid2go_d[raw_gid2go_d$tax_id==tax_id,]
raw_gid2go_d$Num_Genes=go_anot_data$Num_Genes[match(raw_gid2go_d$GO_ID,go_anot_data$GO_ID)]
ifsel=raw_gid2go_d$Category=="Process" & raw_gid2go_d$Num_Genes<=2500; sum(ifsel)


#convert to each gene one row table
input_d=go2gene_d2; input_d_name="gene2allGO"
input_d=raw_gid2go_d[ifsel,]; input_d_name="gene2endGO"

gene2multiGO_d=data.frame(gene_id=unique(input_d$Gene_ID))
gene2multiGO_d$GO_anos=tapply(input_d$GO_term, input_d$Gene_ID, paste, collapse=";")[as.character(gene2multiGO_d$gene_id)]
gene2multiGO_d[c("gene_symbol","gene_desc")]=gid2_gene_desc(gene2multiGO_d$gene_id, spe=organism)[c("gene_symbol","gene_desc")]
gene2multiGO_d=gene2multiGO_d[c("gene_symbol","gene_id","GO_anos","gene_desc")]
out_gene2multiGO_f=paste("3GO_set/",organism,".",input_d_name,".txt",sep="")
write.table(gene2multiGO_d, file=out_gene2multiGO_f, row.names=F, col.names=T, sep="\t",  quote=F)


