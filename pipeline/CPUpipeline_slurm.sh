#!/bin/bash

set -euo pipefail

# CPU pipeline
# May 2018 tjongh@JGM GT
# Backbone: ChIA-PET Utilities (credit to CheeHongWong), macs2, samtools, bedtools, kent, juicer, and powerful BASH :-)
# The Jackson Laboratory for Genomic Medicine

# REQUIRES 2 arguments
configfile=$1 #configuration general to a set of libraries
fastq_list=$2 #fastq_file with RUN id list
#configfile contains parameters of this pipeline
############################# Example of config file ##################################### 
# #Define directories & filenames
# R1suf="_R1_001.fastq.gz" # GT fastq file sufffix
# R2suf="_R2_001.fastq.gz" # read2 suffix
# genome='hg19' #genome version
# 
# #Chiapet params
# extbp=500
# selfbp=8000
# 
# #Peak calling params
# macsq=0.05
# organism='hs' #used in macs: human='hs'  mouse='mm'
# peakext=250
# 
# #Chiasig params
# PET=3
# 
# 
# #cluster job setting parameters:
# NTHREAD=16 #cpu usage per job
# hrs=12 #number of hours allocation
# GB=64 #memory in GB
# seed=12639
# 
# #Not important for running, metadata only
# 
# # The type of sequencing run:
# #    "miseq" - around 30 million reads
# #    "nextseq" - around 300 million reads
# #    "hiseq" - around 300 million reads
# #    "pooled" - around 1 billion reads
# run_type="nextseq"
# 
# # The factor for which the IP was performed
# ip_factor="Pol2"
# 
# # Cell type
# cell_type="cmy"
##########################################################################################


# fastq_list has 2 columns:
# prefix with full path of fastq files & RUN_ID

#$$$$$$$$$$$$$$$ Example of fastq_list $$$$$$$$$$$$$$$$$$$$$$$$$# 
# /projects/wei-lab/proj-in-situ-chia-pet/silvia/WT2_GT18-06298_AGGCAGAA-AGAGGATA_S2 wt2
# /projects/wei-lab/proj-in-situ-chia-pet/silvia/MUT1_GT18-06300_TAAGGCGA-ATAGAGAG_S4 mut1
##########################################################################################

die () {
 echo >&2 "$@"
 exit 1
}
[ $# -gt 1 ] || die "config file and fastq_list in this order required!"



maindir=$PWD

source "$configfile"

: "${R1suf:?R1suf must be defined in config}"
: "${R2suf:?R2suf must be defined in config}"
: "${genome:?genome must be defined in config}"
: "${extbp:?extbp must be defined in config}"
: "${selfbp:?selfbp must be defined in config}"
: "${macsq:?macsq must be defined in config}"
: "${organism:?organism must be defined in config}"
: "${peakext:?peakext must be defined in config}"
: "${NTHREAD:?NTHREAD must be defined in config}"
: "${memMb:?memMb must be defined in config}"
: "${minMapScore:?minMapScore must be defined in config}"
: "${read_id_run_id:?read_id_run_id must be defined in config}"
: "${sif:?sif must be defined in config}"
: "${conda_activate:?conda_activate must be defined in config}"
: "${trim_env:?trim_env must be defined in config}"
: "${python_env_activate:?python_env_activate must be defined in config}"

#tools and refs
: "${fasta:?fasta must be defined in config}"
: "${cpuprog:?cpuprog must be defined in config}"
: "${annot:?annot must be defined in config}"
: "${annotPerl:?annotPerl must be defined in config}"
: "${extractsummary:?extractsummary must be defined in config}"
: "${Rdumpitx_bedpe:?Rdumpitx_bedpe must be defined in config}"
: "${addpath:?addpath must be defined in config}"

#defined names
pairlabel='singlelinker.paired'
singlabel='singlelinker.single'
nonelabel='none'

#------------------------------------

function odir {
outputDir=$1
if [ ! -d $outputDir ]; then
    echo "Create $outputDir"
    mkdir $outputDir
fi
}


function pbsjob {
    fqfile=$1
    RUN=$2
    read1=${fqfile}${R1suf}
    read2=${fqfile}${R2suf}
    pbs="$RUN.insituChiapet.pbs"
    peakfile="${RUN}_peaks.narrowPeak"
    cis_prefix="${RUN}.e${extbp}.clusters.cis.chiasig"
    trans_prefix="${RUN}.e${extbp}.clusters.trans.chiasig"
    chiasigoutput="${cis_prefix}.BE2.sigf.interactions"

    juicer="singularity run ${sif} juicer_tools"
    Samtools="singularity run ${sif} samtools"
    Bedtools="singularity run ${sif} bedtools"
    MACS2="singularity run ${sif} macs2"
    Pigz="singularity run ${sif} pigz"
    
    Trim_galore="singularity run ${sif} trim_galore"
    trim_cutsite="singularity run ${sif} trimcutsite"
    Bioawk="singularity run ${sif} bioawk"
    
    BamCoverage="singularity run ${sif} bamCoverage"
    : "${BedGraphToBigWig:?BedGraphToBigWig must be defined in config}"
    : "${BedSort:?BedSort must be defined in config}"
    
    Java="singularity run ${sif} java"
    Rscript="singularity run ${sif} Rscript"
    Python="singularity run ${sif} python"
    Perl="singularity run ${sif} perl"

#Create directory for each RUN
odir $RUN
cd $RUN

cat << ENDHERE > $pbs
#!/bin/bash 
#SBATCH -N 1 # number of nodes
#SBATCH -n ${NTHREAD} # number of cores
#SBATCH --mem ${memMb} # memory pool for all cores
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

cd \$SLURM_SUBMIT_DIR
echo "Workdir: \$SLURM_SUBMIT_DIR"

# 清理可能冲突的环境模块
module purge

# 使用 Singularity 运行所需工具
module load apps/singularity/3.11.4


LOGFILE=${RUN}.insituChiapet_processing.log

export PATH=$addpath:\$PATH


# perform read ID conversion

parallel -j ${NTHREAD} "
   zcat {} | awk -v run_id='${read_id_run_id}' '
   NR%4==1 {
     printf \"@%s:0100:000%d 1:N:0:GATTACCG+TTGGACTG\\\n\", 
            run_id,  int(NR/4)+1
   }
   NR%4==3 { print \"+\" }
   NR%4!=1 && NR%4!=3 { print \\\$0 }
 ' | gzip > ${RUN}_IDconverted_R1.fastq.gz
 echo \"[info] \$(date) finish converting fastq_1 header of ${RUN} for {}\"
" ::: $read1 &

parallel -j ${NTHREAD} "
   zcat {} | awk -v run_id='${read_id_run_id}' '
   NR%4==1 {
     printf \"@%s:0100:000%d 2:N:0:GATTACCG+TTGGACTG\\\n\", 
            run_id, int(NR/4)+1
   }
   NR%4==3 { print \"+\" }
   NR%4!=1 && NR%4!=3 { print \\\$0 }
 ' | gzip > ${RUN}_IDconverted_R2.fastq.gz
 echo \"[info] \$(date) finish converting fastq_2 header of ${RUN} for {}\"
" ::: $read2

wait

source ${conda_activate}
conda activate ${trim_env}

# FastQC
fastqc ${RUN}_IDconverted_R1.fastq.gz -t ${NTHREAD}
fastqc ${RUN}_IDconverted_R2.fastq.gz -t ${NTHREAD}

# Trim Galore 
#trim_galore --paired --fastqc --stranded_illumina -q 20 --stringency 3 --gzip --length 20 --output_dir . --basename ${RUN}_IDconverted ${RUN}_IDconverted_R1.fastq.gz ${RUN}_IDconverted_R2.fastq.gz

 ${Trim_galore} -q 20 -j ${NTHREAD} --fastqc --stranded_illumina --stringency 3 --gzip --length 20 --output_dir . --paired ${RUN}_IDconverted_R1.fastq.gz ${RUN}_IDconverted_R2.fastq.gz

## Rename
#mv ${RUN}_IDconverted_val_1.fq.gz ${RUN}_IDconverted_R1.trimmed.fq.gz
#mv ${RUN}_IDconverted_val_2.fq.gz ${RUN}_IDconverted_R2.trimmed.fq.gz

echo "trim_galore completed!"
conda deactivate

#rm ${RUN}_IDconverted_R1.fastq.gz 
#rm ${RUN}_IDconverted_R2.fastq.gz


# perform linker detection and generation of different category of fastq files
echo "Linker detection on: ${RUN} " 2>\${LOGFILE}

$cpuprog stag -W -T 18 -t ${NTHREAD} -O ${RUN} ${RUN}_IDconverted_R1_val_1.fq.gz ${RUN}_IDconverted_R2_val_2.fq.gz  2>>\${LOGFILE}

echo "--- linker detection completed ----" >>\${LOGFILE}

#Get the stat
$cpuprog stat -s -p -T 18 -t ${NTHREAD} ${RUN}.cpu 2>>\${LOGFILE} 1>${RUN}.stat
echo "--- statistics done  ----" >>\${LOGFILE}

date  >>\${LOGFILE}

echo "--- pigziiping   ----" >>\${LOGFILE}
${Pigz} -p ${NTHREAD} ${RUN}.singlelinker.paired.fastq 2>>\${LOGFILE}
${Pigz} -p ${NTHREAD} ${RUN}.singlelinker.single.fastq 2>>\${LOGFILE}
${Pigz} -p ${NTHREAD} ${RUN}.none.fastq 2>>\${LOGFILE}
${Pigz} -p ${NTHREAD} ${RUN}.conflict.fastq 2>>\${LOGFILE}
${Pigz} -p ${NTHREAD} ${RUN}.tied.fastq 2>>\${LOGFILE}

#------------------------------------------------------------------------------

#Pair
#-- perform hybrid bwa men and bwa all mapping, de-duplication, span computation, and tag clustering --#
# mapping
echo  START  ${RUN} cpu memaln .. | tee -a \${LOGFILE}
echo  Mapping paired tags .. | tee -a \${LOGFILE}
$cpuprog memaln -T $minMapScore -t ${NTHREAD} $fasta ${RUN}.$pairlabel.fastq.gz 2>> \${LOGFILE} 1>${RUN}.$pairlabel.sam

${Pigz} -p ${NTHREAD} ${RUN}.$pairlabel.sam | tee -a \${LOGFILE}
echo  ENDED pair mapping | tee -a \${LOGFILE}


#pairing
echo  STARTED ${RUN} cpu pair .. | tee -a \${LOGFILE}
echo  Pairing paired tags .. | tee -a \${LOGFILE}
$cpuprog pair -s $selfbp -S -t ${NTHREAD} ${RUN}.$pairlabel.sam.gz 1>${RUN}.$pairlabel.stat.xls 2>> \${LOGFILE}
echo  ENDED ${RUN} cpu pair .. | tee -a \${LOGFILE}

# span
echo  STARTED ${RUN} cpu span .. | tee -a \${LOGFILE}
echo  Computing span of paired tags .. | tee -a \${LOGFILE}
$cpuprog span -s $selfbp -g -t ${NTHREAD} ${RUN}.$pairlabel.UU.bam 2>> \${LOGFILE} 1>${RUN}.$pairlabel.UU.span.xls
echo  ENDED ${RUN} span pair .. | tee -a \${LOGFILE}

# deduplication
echo  STARTED ${RUN} cpu dedup .. | tee -a \${LOGFILE}
echo  De-duplicating paired tags UU .. | tee -a \${LOGFILE}
$cpuprog dedup -g -t ${NTHREAD} ${RUN}.$pairlabel.UU.bam 1>${RUN}.$pairlabel.UU.dedup.lc 2>> \${LOGFILE}
echo  ENDED ${RUN} cpu dedup .. | tee -a \${LOGFILE}

# deduplicated span
echo  STARTED ${RUN} cpu dedup span.. | tee -a \${LOGFILE}
echo  Computing span of paired tags UU nr .. | tee -a \${LOGFILE}
$cpuprog span -s $selfbp -t ${NTHREAD} ${RUN}.$pairlabel.UU.nr.bam 2>> \${LOGFILE} 1>${RUN}.$pairlabel.UU.nr.span.xls
echo  ENDED ${RUN} cpu dedup span.. | tee -a \${LOGFILE}

# cluster tags
echo  STARTED ${RUN} cpu clustering.. | tee -a \${LOGFILE}
$cpuprog cluster -s $selfbp -M -B 1000 -5 5,-20 -3 5,480 -t ${NTHREAD} -O ${RUN}.e$extbp -j -x -v 1 -g ${RUN}.$pairlabel.UU.nr.bam 2>&1 | tee -a \${LOGFILE}
echo  ENDED ${RUN} cpu clustering.. | tee -a \${LOGFILE}

# Heatmap for Juicebox
${juicer} pre -r  2500000,1000000,500000,250000,100000,50000,25000,10000 ${RUN}.e${extbp}.juice.gz $RUN.hic $genome  | tee -a \${LOGFILE}

#------------------------------------------------------------------------------
#None tag

# mapping
echo  STARTED ${RUN}.$nonelabel cpu memaln .. | tee -a \${LOGFILE}
$cpuprog memaln -T 15 -t ${NTHREAD} $fasta ${RUN}.$nonelabel.fastq.gz 2>> \${LOGFILE} 1>${RUN}.$nonelabel.sam
pigz -p ${NTHREAD} ${RUN}.$nonelabel.sam 2>&1  | tee -a \${LOGFILE}
echo  ENDED ${RUN} cpu memaln .. | tee -a \${LOGFILE}

# pairing
echo Pairing $nonelabel tags .. | tee -a \${LOGFILE}
$cpuprog pair -S -t ${NTHREAD} ${RUN}.$nonelabel.sam.gz 1>${RUN}.$nonelabel.stat.xls 2>> \${LOGFILE}
echo  ENDED ${RUN} cpu pair .. | tee -a \${LOGFILE}

# span
echo  STARTED ${RUN}.$nonelabel cpu span .. | tee -a \${LOGFILE}
$cpuprog span -g -t ${NTHREAD} ${RUN}.$nonelabel.UU.bam 2>> \${LOGFILE} 1>${RUN}.$nonelabel.UU.span.xls
echo  ENDED ${RUN}.$nonelabel span pair .. | tee -a \${LOGFILE}

# deduplication
echo  STARTED ${RUN}.$nonelabel cpu dedup .. | tee -a \${LOGFILE}
$cpuprog dedup -g -t ${NTHREAD} ${RUN}.$nonelabel.UU.bam 1>${RUN}.$nonelabel.UU.dedup.lc 2>> \${LOGFILE}


#------------------------------------------------------------------------------

#1tag

# mapping
echo STARTED ${RUN}.$singlabel cpu memaln .. | tee -a \${LOGFILE}
$cpuprog memaln -T 15 -t ${NTHREAD} $fasta ${RUN}.$singlabel.fastq.gz 2>> \${LOGFILE} 1>${RUN}.$singlabel.sam
${Pigz} -p ${NTHREAD} ${RUN}.$singlabel.sam 2>&1  | tee -a \${LOGFILE}

# pairing to get bam
echo Pairing $singlabel tags .. | tee -a \${LOGFILE}
$cpuprog pair -S -t ${NTHREAD} ${RUN}.$singlabel.sam.gz 1>${RUN}.$singlabel.stat.xls 2>> \${LOGFILE}

# deduplication
echo STARTED ${RUN}.$singlabel cpu dedup .. | tee -a \${LOGFILE}
$cpuprog dedup -g -t ${NTHREAD} ${RUN}.$singlabel.UxxU.bam 1>${RUN}.$singlabel.UxxU.dedup.lc 2>> \${LOGFILE}
echo  ENDED ${RUN}.$singlabel cpu dedup .. | tee -a \${LOGFILE}


#------------------------------------------------------------------------------

date  >>\${LOGFILE}

#Peak calling

echo Converting BAM to BED mapped-paired . | tee -a \${LOGFILE}
${Bedtools} bamtobed -i ${RUN}.$pairlabel.UU.nr.bam > ${RUN}.sl.p.uu.nr.bed 2>> \${LOGFILE}

echo Converting BAM to BED none | tee -a \${LOGFILE}
${Bedtools} bamtobed -i ${RUN}.$nonelabel.UU.nr.bam > ${RUN}.none.uu.nr.bed 2>> \${LOGFILE}

echo Converting BAM to BED 1tag | tee -a \${LOGFILE}
${Bedtools} bamtobed -i ${RUN}.$singlabel.UxxU.nr.bam > ${RUN}.1tag.nr.bed 2>> \${LOGFILE}


echo Peak calling using MACS2.0.. | tee -a \${LOGFILE}
${MACS2} callpeak --keep-dup all --nomodel --extsize $peakext -B --SPMR -t ${RUN}.sl.p.uu.nr.bed ${RUN}.none.uu.nr.bed ${RUN}.1tag.nr.bed  \
-f BED -g $organism -n ${RUN} --qvalue ${macsq} | tee -a \${LOGFILE}

echo ENDED ${RUN} cpu peak called .. | tee -a \${LOGFILE}

echo  Generating coverage density.. | tee -a \${LOGFILE}
${Samtools} view -H ${RUN}.$pairlabel.UU.nr.bam | grep '^@SQ' | cut -f 2-3 | sed s?SN:?? | sed s?LN:?? > ${RUN}.genome.length


source ${python_env_activate}
bedClip ${RUN}_treat_pileup.bdg ${RUN}.genome.length ${RUN}_treat_pileup.clip.bdg
sort -k1,1 -k2,2n ${RUN}_treat_pileup.clip.bdg > ${RUN}_treat_pileup.clip.sorted.bdg
conda deactivate

${BedGraphToBigWig} ${RUN}_treat_pileup.clip.sorted.bdg ${RUN}.genome.length ${RUN}.treat_pileup.NDP.bw
echo ENDED ${RUN} coverage density generated. | tee -a \${LOGFILE}

echo -e "Converting file formats..\n" >> \${LOGFILE}
${Samtools} sort -@ $NTHREAD -o ${RUN}.$pairlabel.nr.sorted.bam ${RUN}.$pairlabel.UU.nr.bam  
${Samtools} sort -@ $NTHREAD -o ${RUN}.$singlabel.nr.sorted.bam ${RUN}.$singlabel.UxxU.nr.bam  
${Samtools} sort -@ $NTHREAD -o ${RUN}.$nonelabel.nr.sorted.bam ${RUN}.$nonelabel.UU.nr.bam  
${Samtools} merge -@ $NTHREAD ${RUN}.forBASIC.bam ${RUN}.$pairlabel.nr.sorted.bam ${RUN}.$singlabel.nr.sorted.bam ${RUN}.$nonelabel.nr.sorted.bam

# Make bedgraph
${Bedtools} genomecov -ibam ${RUN}.forBASIC.bam -bg > ${RUN}.forBASIC.bedgraph

# Sort bedgraph
${BedSort} ${RUN}.forBASIC.bedgraph ${RUN}.forBASIC.sorted.bedgraph

# Make bigwig
${BedGraphToBigWig} ${RUN}.forBASIC.sorted.bedgraph  ${RUN}.genome.length ${RUN}.coverage.forBASIC.bw
echo -e  "Done converting file formats, ${RUN}.coverage.forBASIC.bw generated.\n" >> \${LOGFILE}

#------------------------------------------------------------------------------
date  >>\${LOGFILE}


###############get iPET > 1 ##################################

echo  Creating BE2 file  | tee -a \${LOGFILE}
zcat ${cis_prefix}.gz | awk '{ if ( \$7 > 1 ) print }' > ${cis_prefix}.BE2

perl -F'\t' -lane '
    # 在第6列后插入一列0.0，并在末尾追加一列0.0
    splice @F, 6, 0, "0.0";
    push @F, "0.0";
    print join("\t", @F);
' ${cis_prefix}.BE2 > ${cis_prefix}.BE2.4anno.tmp

# 在开头添加四行注释
perl -i -pe '
    print "#no chiasig\n#for anno\n#fake FDR\n#fake exp\n" if $. == 1;
' ${cis_prefix}.BE2.4anno.tmp

# 重命名最终文件
mv ${cis_prefix}.BE2.4anno.tmp ${cis_prefix}.BE2.4anno

#BE2 Annotation
echo "BE2 Annotating..."
source ${python_env_activate}
perl $annotPerl  ${cis_prefix}.BE2.4anno $annot $peakfile ${cis_prefix}.BE2.anno | tee -a \${LOGFILE}
conda deactivate

#dump out itx
echo  "Creating iPET>1 itx files"
${Rscript} ${Rdumpitx_bedpe} ${cis_prefix}.BE2.anno | tee -a \${LOGFILE}


###############get iPET > 2 ##################################

echo  Creating BE3 file  | tee -a \${LOGFILE}
zcat ${cis_prefix}.gz | awk '{ if ( \$7 > 2 ) print }' > ${cis_prefix}.BE3

perl -F'\t' -lane '
    # 在第6列后插入一列0.0，并在末尾追加一列0.0
    splice @F, 6, 0, "0.0";
    push @F, "0.0";
    print join("\t", @F);
' ${cis_prefix}.BE3 > ${cis_prefix}.BE3.4anno.tmp

# 在开头添加四行注释
perl -i -pe '
    print "#no chiasig\n#for anno\n#fake FDR\n#fake exp\n" if $. == 1;
' ${cis_prefix}.BE3.4anno.tmp

# 重命名最终文件
mv ${cis_prefix}.BE3.4anno.tmp ${cis_prefix}.BE3.4anno

#BE3 Annotation
echo "BE3 Annotating..."
source ${python_env_activate}
perl $annotPerl  ${cis_prefix}.BE3.4anno $annot $peakfile ${cis_prefix}.BE3.anno | tee -a \${LOGFILE}
conda deactivate

#dump out itx
echo  "Creating iPET>2 itx files"
${Rscript} ${Rdumpitx_bedpe} ${cis_prefix}.BE3.anno | tee -a \${LOGFILE}


##Significant interaction
##run chiasig
##$chiasigprog -c $PET -t $NTHREAD -p -m $selfbp ${cis_prefix}.BE2
#
#date  >>\${LOGFILE}
#
##Annotation
#echo "Annotating..."
# source the configured annotation environment if the optional ChiaSig block is enabled
#perl $annotPerl  $chiasigoutput $annot $peakfile $RUN | tee -a \${LOGFILE}
#conda deactivate

#dump out itx
#echo  Creating itx files  | tee -a \${LOGFILE}
#${Rscript} ${Rdumpitx} $RUN | tee -a \${LOGFILE}


#Get the summary
$extractsummary --conf $maindir/$configfile -r $RUN


#date  >>\${LOGFILE}
##generate matrix
#echo "Generating matrix"
#cat ${cis_prefix}.gz ${trans_prefix}.gz > $RUN.pooled.clusters.gz
#${Python} $matrixprog ${RUN}.pooled.clusters.gz $RUN $binsize
#
#date  >>\${LOGFILE}
#
##TIF CIF
#echo "Calculating normalized matrix and print TIF & CIF"
#${Rscript} $Rciftif $RUN $binsize

#Cleaning files
rm ${RUN}.*.fastq.gz 2>> \${LOGFILE}
rm ${RUN}.*.cluster 2>> \${LOGFILE}
rm ${RUN}.*.dedup 2>> \${LOGFILE}
rm ${RUN}.*.nr.bed 2>> \${LOGFILE}
rm ${RUN}_treat_pileup.clip.bdg ${RUN}_treat_pileup.bdg ${RUN}.genome.length 2>> \${LOGFILE}
rm ${RUN}.*dup.bam ${RUN}.*.self.*.bam ${RUN}.*.discordant.bam ${RUN}.*.nn.bam ${RUN}.*.xx.bam 2>> \${LOGFILE}
rm $RUN.*.sam.gz $RUN.*.intra.bam $RUN.*.inter.bam $RUN.*.failed.bam $RUN.*.nr.bam
rm $RUN.bedpe.anchors.bed $RUN.cis.sigf.anchors.intersect.bed $RUN.bedpe.anchors.gene.bed
rm ${RUN}.forBASIC.bedgraph*
rm ${RUN}_control_lambda.bdg
rm $RUN.pooled.clusters.gz
rm *anchors*.bed 2>> \${LOGFILE}
${Pigz} -p ${NTHREAD} $RUN.*.bedgraph
${Pigz} -p ${NTHREAD} $RUN.*.bdg #keep the clip.sorted.bdg

#------------------------------------------------------------------------------

ENDHERE

sbatch $pbs
cd $maindir 

}


#------- loop & submit
while IFS=$'\t' read -r fqfile RUN; do
   [[ -z "${fqfile}" || "${fqfile:0:1}" == "#" ]] && continue
   pbsjob "$fqfile" "$RUN"
done < "$fastq_list"
