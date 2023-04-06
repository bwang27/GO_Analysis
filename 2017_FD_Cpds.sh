
####1, use STAR to map 
####2, analyze RNA-seq data
####3, create bigwig files for visualization:
####5, cut genes to blocks (exon and intron part) and do expression analysis
####5.5, study all exon splicing using PSI
####6, gene expression analysis
####7, GO analysis for gene expression
####10, do cis elements analysis

while pgrep -xc "perl|R|Rscript"; do  sleep 1m; echo "."; done; echo "finished old work!";

function generate_sample_pairs {
  sample_pairs=()
  for i in `seq 1 ${#test_samples[*]}`; do
    indx=$(($i-1));
    sample_pairs[$indx]="${test_samples[$indx]}_${ref_samples[$indx]}";
  done
  echo ${sample_pairs[*]}
}

root_dir=/drive2/wli/
study_name=2017_FD_Cpds

processedData_dir=/home/nas02/ProcessedData/RNAseq/$study_name/
mkdir -p $processedData_dir/fastq/
src_dir=/home/nas02/RawData/RNAseq/2017/2017_FD_Cpds/Samples.2017-08-25_26213/20171007-H73154.FASTQs

cd $processedData_dir/fastq/

ln $src_dir/FDPCDMSO-01_L1.D707_505_1.clipped.fastq.gz  DMSO.01.1.fq.gz
ln $src_dir/FDPCDMSO-01_L1.D707_505_2.clipped.fastq.gz  DMSO.01.2.fq.gz
ln $src_dir/FDPCDMSO-02_L1.D708_505_1.clipped.fastq.gz  DMSO.02.1.fq.gz
ln $src_dir/FDPCDMSO-02_L1.D708_505_2.clipped.fastq.gz  DMSO.02.2.fq.gz
ln $src_dir/FDPCp5403-1000-03_L1.D709_505_1.clipped.fastq.gz  p5403.1000_1.1.fq.gz
ln $src_dir/FDPCp5403-1000-03_L1.D709_505_2.clipped.fastq.gz  p5403.1000_1.2.fq.gz
ln $src_dir/FDPCp5403-1000-04_L1.D710_505_1.clipped.fastq.gz  p5403.1000_2.1.fq.gz
ln $src_dir/FDPCp5403-1000-04_L1.D710_505_2.clipped.fastq.gz  p5403.1000_2.2.fq.gz
ln $src_dir/FDPCp7671-300-07_L1.D701_506_1.clipped.fastq.gz  p7671.300_1.1.fq.gz
ln $src_dir/FDPCp7671-300-07_L1.D701_506_2.clipped.fastq.gz  p7671.300_1.2.fq.gz
ln $src_dir/FDPCp7671-300-08_L1.D702_506_1.clipped.fastq.gz  p7671.300_2.1.fq.gz
ln $src_dir/FDPCp7671-300-08_L1.D702_506_2.clipped.fastq.gz  p7671.300_2.2.fq.gz
ln $src_dir/FDPCp7671-30-05_L1.D711_505_1.clipped.fastq.gz  p7671.30_1.1.fq.gz
ln $src_dir/FDPCp7671-30-05_L1.D711_505_2.clipped.fastq.gz  p7671.30_1.2.fq.gz
ln $src_dir/FDPCp7671-30-06_L1.D712_505_1.clipped.fastq.gz  p7671.30_2.1.fq.gz
ln $src_dir/FDPCp7671-30-06_L1.D712_505_2.clipped.fastq.gz  p7671.30_2.2.fq.gz
ln $src_dir/FDPCp7744-300-09_L1.D703_506_1.clipped.fastq.gz  p7744.300_1.1.fq.gz
ln $src_dir/FDPCp7744-300-09_L1.D703_506_2.clipped.fastq.gz  p7744.300_1.2.fq.gz
ln $src_dir/FDPCp7744-300-10_L1.D704_506_1.clipped.fastq.gz  p7744.300_2.1.fq.gz
ln $src_dir/FDPCp7744-300-10_L1.D704_506_2.clipped.fastq.gz  p7744.300_2.2.fq.gz
ln $src_dir/FDPCp7750-30-11_L1.D705_506_1.clipped.fastq.gz  p7750.30_1.1.fq.gz
ln $src_dir/FDPCp7750-30-11_L1.D705_506_2.clipped.fastq.gz  p7750.30_1.2.fq.gz
ln $src_dir/FDPCp7750-30-12_L1.D706_506_1.clipped.fastq.gz  p7750.30_2.1.fq.gz
ln $src_dir/FDPCp7750-30-12_L1.D706_506_2.clipped.fastq.gz  p7750.30_2.2.fq.gz
ln $src_dir/FDPCp7808-30-13_L1.D707_506_1.clipped.fastq.gz  p7808.30_1.1.fq.gz
ln $src_dir/FDPCp7808-30-13_L1.D707_506_2.clipped.fastq.gz  p7808.30_1.2.fq.gz
ln $src_dir/FDPCp7808-30-14_L1.D708_506_1.clipped.fastq.gz  p7808.30_2.1.fq.gz
ln $src_dir/FDPCp7808-30-14_L1.D708_506_2.clipped.fastq.gz  p7808.30_2.2.fq.gz
ln $src_dir/FDPCp7814-300-15_L1.D709_506_1.clipped.fastq.gz  p7814.300_1.1.fq.gz
ln $src_dir/FDPCp7814-300-15_L1.D709_506_2.clipped.fastq.gz  p7814.300_1.2.fq.gz
ln $src_dir/FDPCp7814-300-16_L1.D710_506_1.clipped.fastq.gz  p7814.300_2.1.fq.gz
ln $src_dir/FDPCp7814-300-16_L1.D710_506_2.clipped.fastq.gz  p7814.300_2.2.fq.gz



####1, mapping
root_dir=/drive2/wli/


##using STAR to map
study_name=2017_FD_Cpds
geno=hg19
samples="DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2"

processedData_dir=/home/nas02/ProcessedData/RNAseq/$study_name/
cd $processedData_dir
mkdir STAR_map
fq_dir=fastq/
genomeDir=$root_dir/analyze/dep_seq/1mapping/7map_index/STAR/$geno.noJunc #genome index with no pre-specified junction information
genomeDir=$root_dir/analyze/dep_seq/1mapping/7map_index/STAR/$geno

##unzip
p=0
for sample1 in $samples; do
  fq1=$fq_dir/${sample1}.1.fq.gz
  fq2=$fq_dir/${sample1}.2.fq.gz
  gunzip $fq1 &
  gunzip $fq2 &
  p=$(($p+2)); echo $p; if [ $p -ge 8 ]; then  p=0 ;  wait; fi
done
wait


samples="DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2"
for sample1 in $samples; do
  fq1=$fq_dir/${sample1}.1.fq
  fq2=$fq_dir/${sample1}.2.fq
  out_dir=STAR_map/$sample1/
  mkdir $out_dir
  STAR --runMode alignReads -f --genomeDir $genomeDir --readFilesIn $fq1 $fq2 --runThreadN 5 --outSAMtype BAM SortedByCoordinate \
    --genomeLoad LoadAndKeep --limitBAMsortRAM  10000000000 \
    --outSAMstrandField intronMotif   --outSAMattributes NH HI NM MD jM jI XS --outFileNamePrefix $out_dir
done


# /home/wli/soft/biosoft/seq_map/STAR/STAR-STAR_2.4.2a/bin/Linux_x86_64


# --genomeLoad NoSharedMemory
STAR  --genomeDir $genomeDir --genomeLoad Remove  #release memory

cd STAR_map
grep "Number of input reads" */Log.final.out >Input.read.num.log.out





##create index for IGV visualization
cd $processedData_dir/STAR_map
p=0
for sample1 in $samples; do
  ln $sample1/Aligned.sortedByCoord.out.bam $sample1.bam
  echo $sample1
  samtools index $sample1.bam &
  p=$(($p+1)); echo $p; if [ $p -ge 7 ]; then  p=0 ;  wait; fi
done
wait





####2, analyze RNA-seq data

####2.1  find uniquely mapped
study_name=2017_FD_Cpds
samples="DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2"
cd $root_dir/analyze/dep_seq/2gene_readsnum

p=0
for sample in $samples; do
  perl 6cal_geno_junc_readsnum_frSAM.pl -n $study_name -s "$sample"  -d "$processedData_dir/STAR_map/" \
    -e ".bam" -t STAR -q 10 -m 5 -u 1 -p SE &
  p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait


cd $root_dir/analyze/dep_seq/2gene_readsnum/6out/
perl 1cal_readnum.pl -n $study_name


####2.2 calculate junction read map to gene info (as well as gene expression based on CDS region of refseq transcripts)
geno=hg19
study_name=2017_FD_Cpds
cd $root_dir/analyze/dep_seq/2gene_readsnum
samples="DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2"
p=0
for sample in $samples; do
  perl 14cal_gene_cds_readsnum.pl -s $study_name -g $geno -i refseqid -p refseqcds_SE. -e CDS -r 1 -d 1 -a "$sample" \
    -D 6out/$study_name/ -j 1 -E "JunReadNum GenoMapNum"  &
    p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait

cd $root_dir/analyze/dep_seq/2gene_readsnum/14gene_Rnum
perl 1cal_readnum.pl -n $study_name -p refseqcds_SE




####2.4 combine junction read mapping info and do splicing analysis using DEXSeq
geno=hg19
study_name=2017_FD_Cpds
cd $root_dir/analyze/dep_seq/2gene_readsnum
samples="DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2"
Rscript 26comb_junc_map_info.R -study_name $study_name -geno $geno -samples "$samples" -indir 14gene_Rnum/$study_name/refseqcds_SE.




####3, create bigwig files for visualization:
study_name=2017_FD_Cpds
geno=hg19
PROG_DIR=$root_dir/soft/biosoft/genome/UCSC
genofasta=$root_dir/data/ucsc/genomes/$geno/$geno.fa
fext='.bam.unique' #string after sample name in input sam file name

function sam2bw {
 totalReadNum=`wc -l $sample$fext | sed s/[[:blank:]].*//`
 echo "step2, for file $sample$fext, TotalReadNum=$totalReadNum\n"
 echo step3, ____ $sample, convert to bam ____
 samtools view -S -b -T $genofasta -o $sample$fext.bam  $sample$fext 
 echo step4,  ____ $sample, sort bam file ____
 #samtools sort $sample$fext.bam  $sample$fext.bam.sorted ##!!!no .bam in output name
 echo step5, ____ $sample, convert to bedgraph ____
 genomeCoverageBed -bg -split -ibam $sample$fext.bam -g $PROG_DIR/$geno.chrom.sizes > $sample.bedgraph
 #genomeCoverageBed -bg -split -ibam $sample$fext.bam.2.bam -strand + -g $PROG_DIR/$geno.chrom.sizes > $sample.p.bedgraph #-split 
 #genomeCoverageBed -bg -split -ibam $sample$fext.bam.2.bam -strand '-' -g $PROG_DIR/$geno.chrom.sizes > $sample.m.bedgraph #-split 
 echo step6, ____ $sample, normolize bedgraph counts ____
 norm_bedgraph.pl -t $totalReadNum -i "$sample.bedgraph" 
 echo step8, ____ $sample, convert to bw ____
 $PROG_DIR/bedGraphToBigWig $sample.bedgraph.normolized $PROG_DIR/$geno.chrom.sizes $sample.bw 
 #$PROG_DIR/bedGraphToBigWig $sample.p.bedgraph.normolized $PROG_DIR/$geno.chrom.sizes $sample.p.bw 
 #$PROG_DIR/bedGraphToBigWig $sample.m.bedgraph.normolized $PROG_DIR/$geno.chrom.sizes $sample.m.bw 
}


samples="DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2"
p=0
cd $processedData_dir/STAR_map/
for sample in  $samples; do 
  sam2bw  &
  p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait

cd $processedData_dir/STAR_map/
rm *.bedgraph
rm *.bedgraph.normolized
rm *.bam.unique
rm *.unique.bam.2




####create track strings
study_name=2017_FD_Cpds
colors=(000,000,125 000,000,255 125,000,150 125,000,020 255,000,000)
samples=(DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2) 
col_indexs=(0 0 1 1 2 2 2 2 3 3 2 2 3 3 4 4); priority_start=41

 for samp_i in `seq 1 ${#samples[*]}`; do
  indx=$(($samp_i-1))
  priority=$(( $priority_start+$indx ))
  echo track type=bigWig visibility=2 alwaysZero=on color=${colors[ ${col_indexs[$indx]} ]}   graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=$study_name \
   priority=$priority name=\"${samples[$indx]}\" description=\"RNA-seq ${samples[$indx]}\" \
   bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/$study_name/${samples[$indx]}.bw
 done

track type=bigWig visibility=2 alwaysZero=on color=000,000,125 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=41 name="DMSO.01" description="RNA-seq DMSO.01" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/DMSO.01.bw
track type=bigWig visibility=2 alwaysZero=on color=000,000,125 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=42 name="DMSO.02" description="RNA-seq DMSO.02" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/DMSO.02.bw
track type=bigWig visibility=2 alwaysZero=on color=000,000,255 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=43 name="p5403.1000_1" description="RNA-seq p5403.1000_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p5403.1000_1.bw
track type=bigWig visibility=2 alwaysZero=on color=000,000,255 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=44 name="p5403.1000_2" description="RNA-seq p5403.1000_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p5403.1000_2.bw
track type=bigWig visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=45 name="p7671.30_1" description="RNA-seq p7671.30_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7671.30_1.bw
track type=bigWig visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=46 name="p7671.30_2" description="RNA-seq p7671.30_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7671.30_2.bw
track type=bigWig visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=47 name="p7671.300_1" description="RNA-seq p7671.300_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7671.300_1.bw
track type=bigWig visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=48 name="p7671.300_2" description="RNA-seq p7671.300_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7671.300_2.bw
track type=bigWig visibility=2 alwaysZero=on color=125,000,020 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=49 name="p7744.300_1" description="RNA-seq p7744.300_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7744.300_1.bw
track type=bigWig visibility=2 alwaysZero=on color=125,000,020 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=50 name="p7744.300_2" description="RNA-seq p7744.300_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7744.300_2.bw
track type=bigWig visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=51 name="p7750.30_1" description="RNA-seq p7750.30_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7750.30_1.bw
track type=bigWig visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=52 name="p7750.30_2" description="RNA-seq p7750.30_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7750.30_2.bw
track type=bigWig visibility=2 alwaysZero=on color=125,000,020 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=53 name="p7808.30_1" description="RNA-seq p7808.30_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7808.30_1.bw
track type=bigWig visibility=2 alwaysZero=on color=125,000,020 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=54 name="p7808.30_2" description="RNA-seq p7808.30_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7808.30_2.bw
track type=bigWig visibility=2 alwaysZero=on color=255,000,000 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=55 name="p7814.300_1" description="RNA-seq p7814.300_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7814.300_1.bw
track type=bigWig visibility=2 alwaysZero=on color=255,000,000 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=56 name="p7814.300_2" description="RNA-seq p7814.300_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7814.300_2.bw
##UCSC genome browser track string (bigbed format)
#track type=bigBed itemRgb=On visibility=3 colorByStrand='255,0,0 0,0,255' group=pA priority=1 name=GeneBlock description="2017_FD_Cpds Gene Blocks"  bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/geneBlock.flat.tbl.ano.bed.sorted.bb

track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=000,000,125 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=41 name="DMSO.01" description="RNA-seq DMSO.01" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/DMSO.01.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=000,000,125 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=42 name="DMSO.02" description="RNA-seq DMSO.02" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/DMSO.02.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=000,000,255 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=43 name="p5403.1000_1" description="RNA-seq p5403.1000_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p5403.1000_1.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=000,000,255 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=44 name="p5403.1000_2" description="RNA-seq p5403.1000_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p5403.1000_2.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=45 name="p7671.30_1" description="RNA-seq p7671.30_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7671.30_1.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=46 name="p7671.30_2" description="RNA-seq p7671.30_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7671.30_2.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=47 name="p7671.300_1" description="RNA-seq p7671.300_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7671.300_1.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=48 name="p7671.300_2" description="RNA-seq p7671.300_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7671.300_2.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=125,000,020 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=49 name="p7744.300_1" description="RNA-seq p7744.300_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7744.300_1.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=125,000,020 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=50 name="p7744.300_2" description="RNA-seq p7744.300_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7744.300_2.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=51 name="p7750.30_1" description="RNA-seq p7750.30_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7750.30_1.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=52 name="p7750.30_2" description="RNA-seq p7750.30_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7750.30_2.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=125,000,020 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=53 name="p7808.30_1" description="RNA-seq p7808.30_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7808.30_1.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=125,000,020 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=54 name="p7808.30_2" description="RNA-seq p7808.30_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7808.30_2.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=255,000,000 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=55 name="p7814.300_1" description="RNA-seq p7814.300_1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7814.300_1.bw
track type=bigWig autoScale=off yLineOnOff=on viewLimits=0:3 visibility=2 alwaysZero=on color=255,000,000 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=2017_FD_Cpds priority=56 name="p7814.300_2" description="RNA-seq p7814.300_2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/2017_FD_Cpds/p7814.300_2.bw

#autoScale=off yLineOnOff=on viewLimits=0:6

####5, cut genes to blocks (exon and intron part) and do expression analysis


####5.1 cut gene to blocks based on splice site
study_name=2017_FD_Cpds
geno=hg19
cd $root_dir/analyze/dep_seq/2gene_readsnum
#perl 29cut_gene3blocks_basedOnSplicing.pl -s $study_name -g $geno -R 4 -S 2 -F 0.03 -j 14gene_Rnum/$study_name/refseqcds_SE.combine.junc2gene.tbl
perl 29cut_gene3blocks_basedOnSplicing.pl -s $study_name -g $geno -R 4 -S 2 -F 0.2 -j 14gene_Rnum/$study_name/refseqcds_SE.combine.junc2gene.tbl

##create bed and bigbed file
perl 29b_Gblock2bigBed.pl -s $study_name -g $geno 


####5.2, extract splice site flanking region sequence
cd $root_dir/analyze/club
study_name=2017_FD_Cpds
geno=hg19

#extract 100 bp upstream and downstream of splice site (this input file has more junctions, some of them are from gene model only, not from reads)
perl 23extract_gene_subRegion_seq.pl -d $study_name -g $geno -b 1 -s ../dep_seq/2gene_readsnum/29geneBlocks/$study_name/AllJunc.tbl \
  -i "contig strand juncpos5" -e juncpos5 -n ss5_PM100 -c m99p100 &
 
perl 23extract_gene_subRegion_seq.pl -d $study_name -g $geno -b 1 -s ../dep_seq/2gene_readsnum/29geneBlocks/$study_name/AllJunc.tbl \
  -i "contig strand juncpos3" -e juncpos3 -n ss3_PM100 -c m100p99 &
wait

####5.3, calculat4e gene block expression (exonic and intronic parts) 
geno=hg19
study_name=2017_FD_Cpds
cd $root_dir/analyze/dep_seq/2gene_readsnum
samples="DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2"
p=0
for sample in $samples; do
  perl 14cal_gene_cds_readsnum.pl -s $study_name -g $geno -i Gblock_id -p Gblock. -e transcript -r 1 -d 1 -a "$sample" \
    -D 6out/$study_name/ -j 0 -E "JunReadNum GenoMapNum"  \
    -A 29geneBlocks/$study_name/geneBlock.flat.tbl  &
    p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait
cd $root_dir/analyze/dep_seq/2gene_readsnum/14gene_Rnum
perl 1cal_readnum.pl -n $study_name -p Gblock




####5.4 do splicing analysis using DEXSeq 
cd $root_dir/analyze/dep_seq/2gene_readsnum

study_name=2017_FD_Cpds
geno=hg19


sample_repl_l="DMSO:DMSO.01 DMSO.02;p5403.1000:p5403.1000_1 p5403.1000_2;p7671.30:p7671.30_1 p7671.30_2;p7671.300:p7671.300_1 p7671.300_2;p7744.300:p7744.300_1 p7744.300_2;p7750.30:p7750.30_1 p7750.30_2;p7808.30:p7808.30_1 p7808.30_2;p7814.300:p7814.300_1 p7814.300_2;"
test_samples="p5403.1000 p7671.30 p7671.300 p7744.300 p7750.30 p7808.30 p7814.300"
ref_samples="DMSO DMSO DMSO DMSO DMSO DMSO DMSO"
samples="DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2"
ana_unit=geneBlock 
#note: need relatively large memory, do not use too many workers if memory is small
Rscript 27ana_junc_map_info.R -study_name $study_name -sample_repl_l "$sample_repl_l" \
   -test_samples "$test_samples" -ref_samples "$ref_samples" -all_samples "$samples" \
   -workers 2 \
   -gene_id_header refseqid \
   -in_cbf 14gene_Rnum/$study_name/Gblock.ReadNum.tbl \
   -juncAno_f  29geneBlocks/$study_name/geneBlock.flat.tbl.ano \
   -as_ano_root 29geneBlocks/$study_name/ASano. \
   -geneAno_f ../../club/3add_gene_ano/$geno.refflat.desc \
   -ana_unit $ana_unit \
   -what2do all -use_p_type Padj -P_cut 0.001 -Log2ratio_cut 1
 

  -what2do add_as_ano
  -what2do all -use_p_type Padj -P_cut 0.001 -Log2ratio_cut 1
   -what2do cal_reguType_output -use_p_type Padj -P_cut 0.001 -Log2ratio_cut 1
   -what2do cal_reguType_output -use_p_type Padj -P_cut 0.1 
   -what2do add_ss_liftover -second_spl_name 2016_RO247_RO067_mk -second_spl_update_names "reg_len region_ano Log2R_MK22473uM_MK2DMSO Log2R_MK20671uM_MK2DMSO P_MK22473uM_MK2DMSO P_MK20671uM_MK2DMSO" -second_spl_f 27ana_juncreads/2016_RO247_RO067_mk/geneBlock/combine_d.tbl \
      -ss5liftover_f 29geneBlocks/$study_name/liftover.hg19ToRheMac2/AllJunc.tbl.juncpos5.bed.out -ss3liftover_f 29geneBlocks/$study_name/liftover.hg19ToRheMac2/AllJunc.tbl.juncpos3.bed.out
   -what2do add_ss_liftover -second_spl_name 2016_RO247_RO067_mm -second_spl_update_names "reg_len region_ano Log2R_C2C122473uM_C2C12DMSO Log2R_C2C120671uM_C2C12DMSO P_C2C122473uM_C2C12DMSO P_C2C120671uM_C2C12DMSO" -second_spl_f 27ana_juncreads/2016_RO247_RO067_mm/geneBlock/combine_d.tbl \
      -ss5liftover_f 29geneBlocks/$study_name/liftover.hg19ToMm9/AllJunc.tbl.juncpos5.bed.out -ss3liftover_f 29geneBlocks/$study_name/liftover.hg19ToMm9/AllJunc.tbl.juncpos3.bed.out






####5.5, study all exon splicing using PSI

cd $root_dir/analyze/AS/5AS_deepSeq
study_name=2017_FD_Cpds
geno=hg19


sample_repl_l="DMSO:DMSO.01 DMSO.02;p5403.1000:p5403.1000_1 p5403.1000_2;p7671.30:p7671.30_1 p7671.30_2;p7671.300:p7671.300_1 p7671.300_2;p7744.300:p7744.300_1 p7744.300_2;p7750.30:p7750.30_1 p7750.30_2;p7808.30:p7808.30_1 p7808.30_2;p7814.300:p7814.300_1 p7814.300_2;"
test_samples="p5403.1000 p7671.30 p7671.300 p7744.300 p7750.30 p7808.30 p7814.300"
ref_samples="DMSO DMSO DMSO DMSO DMSO DMSO DMSO"
samples="DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2"

Rscript 4ana_exon_PSI.R -study_name $study_name -sample_repl_l "$sample_repl_l" \
   -test_samples "$test_samples" -ref_samples "$ref_samples" -all_samples "$samples" \
   -basal_PSI_samples "DMSO" \
   -read_num_f ../../dep_seq/2gene_readsnum/14gene_Rnum/$study_name/refseqcds_SE.combine.junc2gene.tbl \
   -workers 8 \
   -use_p_type "Pfisher" -P_cut "0.001" -PSI_cal_rnumMin 5 -delta_PSI_cut 30 \
   -what2do calReguType

   -what2do all
   -what2do cal_PSI
   -what2do calReguType
   -use_p_type "Pfisher" -P_cut "0.01" -PSI_cal_rnumMin 2 -delta_PSI_cut 5 \
   -use_p_type "Pfisher" -P_cut "0.001" -PSI_cal_rnumMin 2 -delta_PSI_cut 10 \

##study all annotated exons
Rscript 4ana_exon_PSI.R -study_name $study_name -sample_repl_l "$sample_repl_l" \
   -test_samples "$test_samples" -ref_samples "$ref_samples" -all_samples "$samples" \
   -basal_PSI_samples "DMSO" \
   -read_num_f ../../dep_seq/2gene_readsnum/14gene_Rnum/$study_name/refseqcds_SE.combine.junc2gene.tbl \
   -workers 6 \
   -use_p_type "Pfisher" -P_cut "0.001" -PSI_cal_rnumMin 5 -delta_PSI_cut 10 \
   -study_exon_type  annotated_exon -what2do all



##using a consolidated gene blocks file
root_dir=/drive2/wli/
pipeline_root=/drive2/wli/analyze/Pipeline/ #this is root directory of the pipeline/scripts (where the setGlobalVars.sh is located)

cd $pipeline_root/RNAseq_Splicing
study_name=2017_FD_Cpds #project name
geno=hg19

Rscript 06.ana_exon_PSI.R -study_name $study_name -sample_repl_l "$sample_repl_l" \
   -test_samples "$test_samples" -ref_samples "$ref_samples" -all_samples "$samples" \
   -basal_PSI_samples "DMSO" \
   -gblock_f "/drive2/wli/analyze/dep_seq/othProject/compound_splicing/1Splicing_comb_tb/FD_Cpds_2015to2019/Gblock.Combine.PSI.P0.001.Ch10.txt" \
   -read_num_f /drive2/wli/analyze/dep_seq/2gene_readsnum/14gene_Rnum/$study_name/refseqcds_SE.combine.junc2gene.tbl \
   -out_root /drive2/wli/analyze/AS/5AS_deepSeq/04exon_PSI/$study_name/consolidated.FD_Cpds_2015to2019/ \
   -workers 3 \
   -use_p_type "Pfisher" -PSI_cal_rnumMin 20 \
   -what2do all -P_cut "0.0000001"  -delta_PSI_cut 30 &






####6.1  extract pair-end reads mapping info
cd $root_dir/analyze/dep_seq/2gene_readsnum
samples="DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2"
study_name=2017_FD_Cpds

p=0
for sample in $samples; do
  perl 10cal_PEreadsCov_Gene_pA_Rnum.pl -s $study_name -S "$sample"  -i "$processedData_dir/STAR_map/" -e ".bam" \
    -q 10  -d 0 -m "" &
  p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait


####6.2 calculate gene expression using 14cal_gene_cds_readsnum.pl
study_name=2017_FD_Cpds
geno=hg19
cd $root_dir/analyze/dep_seq/2gene_readsnum
samples=(DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2)
p=0
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1))
  perl 14cal_gene_cds_readsnum.pl -s $study_name -g $geno -i refseqid -p refseqcds_PE. \
    -e CDS -r 1 -d 1 -a "${samples[indx]}" \
    -D 10PE_reasCov/$study_name/ -E reads.cov.tbl  &
  p=$(($p+1));  if [ "$p" -ge 5 ]; then p=0; wait;  fi
done
wait

cd $root_dir/analyze/dep_seq/2gene_readsnum/14gene_Rnum
perl 1cal_readnum.pl -n $study_name -p refseqcds_PE

####6.3 calculate gene expression regulation (log2 ratio of IP RPKM / IgG RPKM)
study_name=2017_FD_Cpds
geno=hg19
#cd $root_dir/analyze/dep_seq/2gene_readsnum
pipeline_root=/drive2/wli/analyze/Pipeline/ #this is root directory of the pipeline/scripts (where the setGlobalVars.sh is located)
project_root=/drive2/wli/analyze/dep_seq/othProject/2017_FD_Cpds/ #project folder (all the input and output files related to this project)
cd $pipeline_root/RNAseq

gene_rnum_f=/drive2/wli/analyze/dep_seq/2gene_readsnum/14gene_Rnum/2017_FD_Cpds/refseqcds_PE.ReadNum.tbl
rnum_idname=refseqid
sample_repl_l="DMSO:DMSO.01 DMSO.02;p5403.1000:p5403.1000_1 p5403.1000_2;p7671.30:p7671.30_1 p7671.30_2;p7671.300:p7671.300_1 p7671.300_2;p7744.300:p7744.300_1 p7744.300_2;p7750.30:p7750.30_1 p7750.30_2;p7808.30:p7808.30_1 p7808.30_2;p7814.300:p7814.300_1 p7814.300_2;"
test_samples="p5403.1000 p7671.30 p7671.300 p7744.300 p7750.30 p7808.30 p7814.300"
ref_samples="DMSO DMSO DMSO DMSO DMSO DMSO DMSO"
samples="DMSO.01 DMSO.02 p5403.1000_1 p5403.1000_2 p7671.30_1 p7671.30_2 p7671.300_1 p7671.300_2 p7744.300_1 p7744.300_2 p7750.30_1 p7750.30_2 p7808.30_1 p7808.30_2 p7814.300_1 p7814.300_2"
out_prefix=""


/bin/Rscript 05.cal_gexch.template.R -p_cut 0.05 -pval_fun "DESeq2" -foldchange_cut 1.5 -IfCal_RPM 0 \
  -gene_model refseqcds  -geno $geno -study_name $study_name -gene_rnum_f $gene_rnum_f \
  -output_root  $project_root/05.Gene_DE/ \
  -nCores 8 \
  -all_sample_name "$samples" \
  -gene_len_var cds_len -test_samples "$test_samples" -ref_samples "$ref_samples" \
  -g_ano_f "/drive2/wli/analyze/club/3add_gene_ano/hg19.refflat.desc" -rnum_idname "$rnum_idname" \
  -comb_sample_l "$sample_repl_l" \
  -out_prefix "$out_prefix" \
  -what2do  all





####10, do cis elements analysis
cd $root_dir/analyze/motif_ana/B_count_num
study_name=2017_FD_Cpds
geno=hg19

##4.3.1 compare splice site sequence
cd $root_dir/analyze/motif_ana/B_count_num
study_name=2017_FD_Cpds
geno=hg19
test_samples=(p5403.1000 p7671.30 p7671.300 p7744.300 p7750.30 p7808.30 p7814.300)
ref_samples=(DMSO DMSO DMSO DMSO DMSO DMSO DMSO)
generate_sample_pairs

##for PSI data
reguCutOff_str="delta_PSI10_Pfisher0.001" #delta_PSI10_Pfisher0.001 delta_PSI5_Pfisher0.01 
inf=$root_dir/analyze/AS/5AS_deepSeq/04exon_PSI/$study_name/all_exon/formatted_tb/exons.$reguCutOff_str.allReguType.txt
for sample_pair in ${sample_pairs[*]}; do
   perl 6count_num_fr_seq_and_idtables_fisher.pl -d $study_name -v $geno  -e ReguType_$sample_pair -g "UP DN NC" \
     -n 20 \
     -i $inf  \
     -p "PSI.$reguCutOff_str/ss5_m50m5." -h "contig strand end_pos" -s ss5 -q ../../club/23gene_subRegSeq/$study_name/ss5_PM100.fa -l "4 5 6" -H "ss5_m50m5 50|45" -f "if_gt_ag_ss_exon reg_len" -t "1 >50"  &
done

     -p "PSI.$reguCutOff_str/ss5_m50m5." -h "contig strand end_pos" -s ss5 -q ../../club/23gene_subRegSeq/$study_name/ss5_PM100.fa -l "4 5 6" -H "ss5_m50m5 50|45" -f "if_gt_ag_ss_exon reg_len" -t "1 >50"  &
     -p "PSI.$reguCutOff_str/ss5_m4p6." -h "contig strand end_pos" -s ss5 -q ../../club/23gene_subRegSeq/$study_name/ss5_PM100.fa -l "4 6 8 10" -H "ss5_m4p6 96|10" -f if_gt_ag_ss_exon -t "1"  &
     -p "PSI.$reguCutOff_str/ss5_m4m1." -h "contig strand end_pos" -s ss5 -q ../../club/23gene_subRegSeq/$study_name/ss5_PM100.fa -l "4" -H "ss5_m4m1 96|4" -f if_gt_ag_ss_exon -t "1"  &
     -p "PSI.$reguCutOff_str/ss5_p1p6." -h "contig strand end_pos" -s ss5 -q ../../club/23gene_subRegSeq/$study_name/ss5_PM100.fa -l "6" -H "ss5_p1p6 100|6" -f if_gt_ag_ss_exon -t "1"  &
     -p "PSI.$reguCutOff_str/ss5_p7p100." -h "contig strand end_pos" -s ss5 -q ../../club/23gene_subRegSeq/$study_name/ss5_PM100.fa -l "4 5 6" -H "ss5_p7p100 106|94" -f if_gt_ag_ss_exon -t "1"  &
     -p "PSI.$reguCutOff_str/ss3_m40m3." -h "contig strand start_pos" -s ss3 -q ../../club/23gene_subRegSeq/$study_name/ss3_PM100.fa -l "4 5 6" -H "ss3_m40m3 60|38" -f if_gt_ag_ss_exon -t "1"  &
     -p "PSI.$reguCutOff_str/ss3_p1p4." -h "contig strand start_pos" -s ss3 -q ../../club/23gene_subRegSeq/$study_name/ss3_PM100.fa -l "4" -H "ss3_p1p4 100|4" -f if_gt_ag_ss_exon -t "1"  &
     -p "PSI.$reguCutOff_str/ss3_m100m41." -h "contig strand start_pos" -s ss3 -q ../../club/23gene_subRegSeq/$study_name/ss3_PM100.fa -l "4 5 6" -H "ss3_m100m41 0|60" -f if_gt_ag_ss_exon -t "1"  &
     -p "PSI.$reguCutOff_str/ss3_p5p50." -h "contig strand start_pos" -s ss3 -q ../../club/23gene_subRegSeq/$study_name/ss3_PM100.fa -l "4 5 6" -H "ss3_p5p50 105|45" -f "if_gt_ag_ss_exon reg_len" -t "1 >50"  &


#extract results and combine
cd $root_dir/analyze/motif_ana/B_count_num
filenIDs=${sample_pairs[*]}
regions=(ss5_m4m1.ss5 ss5_p1p6.ss5 ss5_m4p6.ss5 ss3_p1p4.ss3 ss3_m40m3.ss3 ss3_m100m41.ss3 ss5_p7p100.ss5)
merlenStrs=(4 6 46810 4 456 456 456); infstr2="_if_gt_ag_ss_exon1"
comb_sum_f=6cisnum_fisher/$study_name/PSI.$reguCutOff_str/summary1.tbl

regions=(ss3_p5p50.ss3 ss5_m50m5.ss5)
merlenStrs=(456 456); infstr2="_if_gt_ag_ss_exon_reg_len1_>50"
comb_sum_f=6cisnum_fisher/$study_name/PSI.$reguCutOff_str/summary2.tbl

echo "" >$comb_sum_f
for compare_2grp in UP.VS.NC DN.VS.NC; do  #
  for regioni in `seq 1 ${#regions[*]}`; do
     indx=$(($regioni-1));
     region1=${regions[indx]}
     merlenStr=${merlenStrs[indx]}
     echo "" >>$comb_sum_f
     echo "regulation cutoff=$reguCutOff_str; comparing two groups=$compare_2grp; region=$region1" >>$comb_sum_f
     out_f1=6cisnum_fisher/$study_name/PSI.$reguCutOff_str/tmp.txt
     Rscript 13extract_combine_results.R -infstr1 6cisnum_fisher/$study_name/PSI.$reguCutOff_str/$region1.ReguType_ \
      -search_patterns ".$compare_2grp" -extractNum 20 -grpNum 3 \
      -infvar1 "$filenIDs" -infstr2 "$infstr2/summary.mer" -infvar2 "$merlenStr" -infstr3 ".top20.tbl" \
      -filenIDs "$filenIDs" -out_f $out_f1
      sed "/#/d" $out_f1 >>$comb_sum_f
    done
done




####7, GO analysis for gene expression
organism=hs
geno=hg19
study_name=2017_FD_Cpds
cd $root_dir/analyze/go/A_enrich
test_samples=(p5403.1000 p7671.30 p7671.300 p7744.300 p7750.30 p7808.30 p7814.300)
ref_samples=(DMSO DMSO DMSO DMSO DMSO DMSO DMSO)
#build full sample_pairs array
generate_sample_pairs

genetb_file="/drive2/wli/analyze/dep_seq/2gene_readsnum/8gexch_p/$study_name/refseqcds.DESeq2.format_tb/allGene.txt"
for sample_pair in ${sample_pairs[*]}; do
  Rscript 1go_fisher.test_exe.R -run_mode fisher -organism $organism -study_name $study_name \
  -expe_name Gex  \
  -gene_grp_name GexType_$sample_pair -gene_groups "dn nc up" -gene_id_var "gene_id" \
  -sel_top_sort_var adjSLog10P_$sample_pair -only_top_reg_num 300 -top2reg_names "low:dn;mid:nc;high:up;" \
  -genetb_file "$genetb_file" &
  p=$(($p+1)); echo $p; if [ $p -ge 8 ]; then  p=0 ;  wait; fi
done

#combine datasets # 
test_samples=(p7671.300 p7744.300 p7750.30 p7808.30 p7814.300)
ref_samples=(DMSO DMSO DMSO DMSO DMSO)
#build full sample_pairs array
generate_sample_pairs

Rscript 1go_fisher.test_exe.R -run_mode combine -organism $organism -low_Qval_cut 0.1  -study_name $study_name -exp_type "Gex.GexType_" \
  -go_f_ext ".dn_VS_nc_VS_up.GO_FisherT.tbl" \
  -out_prefix top10each -combine_meth topNinEachP -eachP_GO_num 10 -low_pval_cut 0.05 -comb_2P "P_up P_dn" -comb_2P_SS_method sum -discard_p "P_nc" -max_row_num 80 \
  -sample_names "${sample_pairs[*]}" &

  -out_prefix top10each -combine_meth topNinEachP -eachP_GO_num 10 -low_pval_cut 0.05 -comb_2P "P_up P_dn" -comb_2P_SS_method sum -discard_p "P_nc" -max_row_num 80 \
  -out_prefix mostSigP -combine_meth row_min_p -low_pval_cut 0.05 -comb_2P "P_up P_dn" -comb_2P_SS_method sum -discard_p "P_nc"  \
  -out_prefix mostSigP.up.dn -combine_meth row_min_p -low_pval_cut 0.05 -discard_p "P_nc"  \

#extract genes for top GOs
Rscript 1go_fisher.test_exe.R -run_mode GeneList -GOnum4genelist 10 -unsel_gene_group "" \
  -organism $organism  -study_name $study_name -genetb_file "$genetb_file" \
  -exp_type "Gex.GexType_"  -gene_id_var "gene_id" \
  -gene_grp_var_name "GexType" -sample_names "${sample_pairs[*]}" \
  -discard_header_pattern "num_"


