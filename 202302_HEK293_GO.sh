
#<<SAMPLE_SETTINGS


#mkdir -p /home/Research_Archive/RawData/RNAseq/public_data/202302_HEK293_strains_RNAseq
#cd /home/Research_Archive/RawData/RNAseq/public_data/202302_HEK293_strains_RNAseq

#samples="SRR20828725 SRR20828726 SRR20828727 SRR20828728 SRR20828729 SRR20828730"
##re-download the following:
#samples="SRR17653682 SRR17653685 SRR17653695"
#for sample in $samples; do
# echo $sample
# fastq-dump -I --split-files --gzip $sample
#done

src_dir=/mnt/research_archive_rawdata_azure/RawData/RNAseq/2023/202302_HEK293_strains_RNAseq/19205-75-03162023 
#mkdir /home/Research_Archive/ProcessedData/RNAseq/202302_HEK293_strains_RNAseq
#mkdir -p /home/Research_Archive/ProcessedData/RNAseq/202302_HEK293_strains_RNAseq/fastq
cd /home/Research_Archive/ProcessedData/RNAseq/202302_HEK293_strains_RNAseq/fastq
#ll $src_dir/*gz


##this pipeline start with fastq files, do mapping (using STAR), gene expression quantification and differential gene expression analysis (using DESeq2)
#programs needed to install in linux:
# perl, R
# STAR
# samtools
# bedtools
# R library: DESeq2, ggplot2, edgeR, BiocParallel, fmsb, gplots, openxlsx, DEXSeq, knitr, pander, rmarkdown.


root_dir=/drive2/wli/
study_name=202302_HEK293_strains_RNAseq #project name
pipeline_root=/drive2/wli/analyze/Pipeline/ #this is root directory of the pipeline/scripts (where the setGlobalVars.sh is located)
project_root=/home/Research_Archive/ProcessedData/RNAseq/202302_HEK293_strains_RNAseq/ #project folder (all the input and output files related to this project)
mapping_data_root=/home/Research_Archive/ProcessedData/RNAseq/202302_HEK293_strains_RNAseq/ #use this directory to save large files (eg. fastq files and mapped bam files). expect a fastq folder with all fastq files under this directory
##create STAR index for hg19 genome +  E1A,E1B sequences 
#cd $root_dir/analyze/dep_seq/1mapping (can't write in wencheng's directory)
geno=hg19 #hg19 mm9
out_folder=/home/Research_Archive/ProcessedData/RNAseq/202302_HEK293_strains_RNAseq/STAR_index/hg19_ad5   ### HELA_HEK293/STAR_index/hg19_ad5
#mkdir -p $out_folder
geno_f=/drive2/wli/analyze/Pipeline/ReferenceDB/ucsc/genomes/$geno/$geno.fa
junc_f=/drive2/wli/analyze/club/0remove_club_redundance/8club_junc_info/1109/$geno.1109.STAR.junc
addi_fa_dir=/drive2/wli/analyze/summary/GeneTherapy/NanoporeSequencing/VirusGenomeSequences
#addi_fa="$addi_fa_dir./Ad5/Human_Ad5_first5k.fa $addi_fa_dir/Helper/pHelper.fa $addi_fa_dir/RepCap/pAAV9-WT.fa"
addi_fa="$addi_fa_dir/Ad5/Human_Ad5_first5k.fa"
#star=/drive2/bwang/software/STAR-2.7.9a/bin/Linux_x86_64
 
#STAR --runMode genomeGenerate --genomeDir $out_folder --genomeFastaFiles $geno_f $addi_fa --runThreadN 5 --sjdbFileChrStartEnd $junc_f &

#wait

source $pipeline_root/setGlobalVars.sh
if_PE=1; #whether this is pair end sequencing ()
if_stranded=1 #0 #1; #whether the RNA library is stranded (contain gene strand information)
if_read1_antiSense=1 #0 #1; #whether the first read is antisense to gene. This is usually true for TruSeq Stranded mRNA Library (second read is sense strand)
if_do_splicing_analysis=1; #whether do splicing analysis

nCPUs=5
nNodes=1
DEXSeq_nCPUs=3
DEG_method=DESeq2 #DESeq2 fisher
Gblock_anaMethod=DEXSeq #DEXSeq FET
rnum_idname=refseqid
GexFoldchange_cut=2; GexP_cut=0.05; GexHigher_avg_RPKM_cutoff=1

geno=hg19
samples=(VPC2.0.1 VPC2.0.2 VPC2.0.3 PTCSUS.2 PTCSUS.3 AAV293.1 AAV293.2 AAV293.3 STATBAXKO.1 STATBAXKO.2 STATBAXKO.3 HEK293.1 HEK293.2 HEK293.3 Expi293.1 Expi293.2 Expi293.3 HEK293F.1 HEK293F.2 HEK293F.3)
sample_repl_l="VPC2.0:VPC2.0.1 VPC2.0.2 VPC2.0.3;PTCSUS:PTCSUS.2 PTCSUS.3;AAV293:AAV293.1 AAV293.2 AAV293.3;STATBAXKO:STATBAXKO.1 STATBAXKO.2 STATBAXKO.3;HEK293:HEK293.1 HEK293.2 HEK293.3;Expi293:Expi293.1 Expi293.2 Expi293.3;HEK293F:HEK293F.1 HEK293F.2 HEK293F.3;"
test_samples=(VPC2.0 Expi293 Expi293 VPC2.0 VPC2.0 HEK293F PTCSUS AAV293 STATBAXKO) #Hep3B.HHT HepG2.HHT.10uM HepG2.HHT.100nM HepG2.HHT.500nM HepG2.HHT.1uM Huh6.HHT Huh7.HHT PLC.HHT)
ref_samples=(Expi293 HEK293F HEK293 HEK293 HEK293F HEK293 HEK293 HEK293 HEK293) #Hep3B.ctrl HepG2.ctrl HepG2.ctrl HepG2.ctrl HepG2.ctrl Huh6.ctrl Huh7.ctrl PLC.ctrl)

#ucsc genome browser track setting
col_indexs=(0 0 0 1 1 1 2 2 2 3 3 3 4 4 0 0 0 1 1 1); #4 4 0 0 1 1 2 2 3 3 4 4 0 0 4 4 0 0 4 4 0 0 4 4);
colors=(000,000,125 000,000,255 125,000,150 125,000,020 255,000,000)
priority_start=41

##PSI analysis setting:
PSIuse_p_type="Pfisher"; PSI_P_cut="0.001"; PSI_cal_rnumMin=20; delta_PSI_cut=10; PSIana_setting=delta_PSI${delta_PSI_cut}_${PSIuse_p_type}$PSI_P_cut

####R.1, mapping
##using STAR to map

cd $mapping_data_root
#mkdir STAR_map
fq_dir=$mapping_data_root/fastq/ #the directory of the fastq files 
#genomeDir=$pipeline_root/ReferenceDB/map_index/STAR/$geno #the directory of mapping index file
genomeDir=$mapping_data_root/STAR_index/hg19_ad5 #the directory of mapping index file
chr_sizes_f=$genomeDir/chrNameLength.txt
##create additional refflat format gene annotation file
addi_fa_dir=/drive2/wli/analyze/summary/GeneTherapy/NanoporeSequencing/VirusGenomeSequences
addi_anno_f="$addi_fa_dir/Ad5/Human_Ad5_first5k.refflat"
combined_gene_anno_f=$project_root/combined_gene_refflat_ano.txt
#i=0

function generate_sample_pairs {
  sample_pairs=()
  for i in `seq 1 ${#test_samples[*]}`; do
    indx=$(($i-1));
    sample_pairs[$indx]="${test_samples[$indx]}_${ref_samples[$indx]}";
  done
  echo ${sample_pairs[*]}
}


####7, GO analysis for gene expression
organism=hs
geno=hg19
study_name=202302_HEK293_strains_RNAseq
cd $root_dir/analyze/go/A_enrich
#test_samples=(p5403.1000 p7671.30 p7671.300 p7744.300 p7750.30 p7808.30 p7814.300)
#ref_samples=(DMSO DMSO DMSO DMSO DMSO DMSO DMSO)
#build full sample_pairs array
generate_sample_pairs

#genetb_file="/drive2/wli/analyze/dep_seq/2gene_readsnum/8gexch_p/$study_name/refseqcds.DESeq2.format_tb/allGene.txt"
genetb_file="/home/Research_Archive/ProcessedData/RNAseq/202302_HEK293_strains_RNAseq/05.Gene_DE/refseqcds.DESeq2.format_tb/allGene.txt"

for sample_pair in ${sample_pairs[*]}; do
  Rscript /drive2/bwang/analysis/RNASeq/202302_HEK293_strains_RNAseq/GO/1go_fisher.test_exe.R -run_mode fisher -organism $organism -study_name $study_name \
  -expe_name Gex  \
  -gene_grp_name GexType_$sample_pair -gene_groups "dn nc up" -gene_id_var "gene_id" \
  -sel_top_sort_var adjSLog10P_$sample_pair -only_top_reg_num 300 -top2reg_names "low:dn;mid:nc;high:up;" \
  -genetb_file "$genetb_file" &
  p=$(($p+1)); echo $p; if [ $p -ge 8 ]; then  p=0 ;  wait; fi
done

#combine datasets #
#test_samples=(p7671.300 p7744.300 p7750.30 p7808.30 p7814.300)
#ref_samples=(DMSO DMSO DMSO DMSO DMSO)
#build full sample_pairs array
generate_sample_pairs

Rscript /drive2/bwang/analysis/RNASeq/202302_HEK293_strains_RNAseq/GO/1go_fisher.test_exe.R -run_mode combine -organism $organism -low_Qval_cut 0.1  -study_name $study_name -exp_type "Gex.GexType_" \
  -go_f_ext ".dn_VS_nc_VS_up.GO_FisherT.tbl" \
  -out_prefix top10each -combine_meth topNinEachP -eachP_GO_num 10 -low_pval_cut 0.05 -comb_2P "P_up P_dn" -comb_2P_SS_method sum -discard_p "P_nc" -max_row_num 80 \
  -sample_names "${sample_pairs[*]}" &

  -out_prefix top10each -combine_meth topNinEachP -eachP_GO_num 10 -low_pval_cut 0.05 -comb_2P "P_up P_dn" -comb_2P_SS_method sum -discard_p "P_nc" -max_row_num 80 \
  -out_prefix mostSigP -combine_meth row_min_p -low_pval_cut 0.05 -comb_2P "P_up P_dn" -comb_2P_SS_method sum -discard_p "P_nc"  \
  -out_prefix mostSigP.up.dn -combine_meth row_min_p -low_pval_cut 0.05 -discard_p "P_nc"  \

#extract genes for top GOs
Rscript /drive2/bwang/analysis/RNASeq/202302_HEK293_strains_RNAseq/GO/1go_fisher.test_exe.R -run_mode GeneList -GOnum4genelist 10 -unsel_gene_group "" \
  -organism $organism  -study_name $study_name -genetb_file "$genetb_file" \
  -exp_type "Gex.GexType_"  -gene_id_var "gene_id" \
  -gene_grp_var_name "GexType" -sample_names "${sample_pairs[*]}" \
  -discard_header_pattern "num_"



