#!/bin/bash
#SBATCH --nodes 1
#SBATCH --time=60:00:00
#SBATCH --mem=50G
#SBATCH --output=/home/moorej3/Job-Logs/jobid_%A.output
#SBATCH --error=/home/moorej3/Job-Logs/jobid_%A.error
#SBATCH --partition=5days

data=GM12878.CHiC
version=v2

#data=$1
#version=$2

echo $data
setDir= /Users/xuxi/Documents/cg/BENGI/Benchmark/Remove-Ambiguous-Pairs.Fixed-Ratio
train=$setDir/$data-Benchmark.$version.txt
totalTrain=$setDir/$data-Benchmark.*.txt
featureDir=/Users/xuxi/Documents/cg/BENGI/Scripts/Supervised-Methods/PEP/Feature-Matrices
outputDir=/Users/xuxi/Documents/cg/BENGI/Scripts/Supervised-Methods/PEP/Results

#enhancers=~/Lab/Target-Gene/Target-Finder/GM12878-Enhancers.bed
ccREs=/Users/xuxi/Documents/cg/BENGI/Benchmark/Annotations/hg19-cCREs.bed
scriptDir=/Users/xuxi/Documents/cg/BENGI/Scripts/Supervised-Methods/PEP
tss=/Users/xuxi/Documents/cg/BENGI/Benchmark/Annotations/TSS.2019.bed

genome=/Users/xuxi/Documents/cg/BENGI/Benchmark/Annotations/hg19.2bit
fimo=~/bin/meme/bin/fimo
motifs=/Users/xuxi/Documents/cg/BENGI/Benchmark/Annotations/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
outputDir=/Users/xuxi/Documents/cg/BENGI/Scripts/Supervised-Methods/PEP/Results
header=~/Lab/Target-Gene/PEP/HOCOMOCO-Motif-Header.txt

mkdir -p /tmp/moorej3/$SLURM_JOBID-$jid
cd /tmp/moorej3/$SLURM_JOBID-$jid

mkdir -p $outputDir

######## Creating Enhancer Feature Matrix ################
if [ ! -f "$featureDir/$data-Enhancer-Feature-Matrix.txt" ]
then
    echo -e "Generating enhancer feature matrix..."
    cat $totalTrain | awk '{print $1}' | sort -u  > cres
    awk 'FNR==NR {x[$1];next} ($5 in x)' cres $ccREs | \
    awk '{print $1 "\t" $2-4000 "\t" $3+4000 "\t" $5}' > enhancers
    /Users/xuxi/Documents/cg/BENGI/Benchmark/Annotations/twoBitToFa $genome enhancers.fa -bed=enhancers
    $fimo --text $motifs enhancers.fa  > fimo.out
    python $scriptDir/count.motifs.py $header enhancers fimo.out \
    > $featureDir/$data-Enhancer-Feature-Matrix.txt
fi

######## Creating TSS Feature Matrix ################
if [ ! -f "$featureDir/$data-TSS-Feature-Matrix.txt" ]
then
    echo -e "Generating tss feature matrix..."
    cat $totalTrain | awk '{print $2}' | sort -u  > genes
    awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | \
        awk '{print $1 "\t" $2-500 "\t" $3+500 "\t" $4}' > tss
    /Users/xuxi/Documents/cg/BENGI/Benchmark/Annotations/twoBitToFa $genome tss.fa -bed=tss
    $fimo --text $motifs tss.fa  > fimo.out
    python $scriptDir/count.motifs.py $header tss fimo.out \
    > $featureDir/$data-TSS-Feature-Matrix.txt
fi

######## Creating Distance Matrix ################
#if [ ! -f "$featureDir/$data-Distance.txt" ]
#then
#    echo -e "Generating distance matrix..."
#    cat $val $train | sort -u > pairs
#    python calculate.distance.py $tss $enhancers pairs > \
#        $featureDir/$data-Distance.txt 
#fi

######## Running Random Forest ################
echo -e "Running Model..."
cat $train | awk '{print $2}' | sort -u  > genes
awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | \
awk '{print $1 "\t" $2-500 "\t" $3+500 "\t" $4 "\t" $7 }' > tss

python $scriptDir/xgb.py $train $featureDir/$data-Enhancer-Feature-Matrix.txt \
    $featureDir/$data-TSS-Feature-Matrix.txt tss $data $outputDir $version 

rm -r /tmp/moorej3/$SLURM_JOBID-$jid
