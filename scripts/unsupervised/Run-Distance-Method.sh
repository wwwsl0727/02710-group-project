#!/bin/bash
setup() {
  # Fill in this helper function to do any setup if you need to.
  #
  # This function will be executed once at the beginning of the grading process.
  # Other functions may be executed repeatedly in arbitrary order.
  # Make use of this function to reduce any unnecessary overhead and to make your
  # code behave consistently.
  #
  # e.g., You should compile any code in this function.
  #
  # However, DO NOT use this function to solve any question or
  # generate any intermediate output.
  #
  # Examples:
  # mvn clean package
  # pip3 install -r requirements.txt
  #
  # Standard output format:
  # No standard output needed

  # convert the notebook into an executable Python script (DO NOT remove)
  jupyter nbconvert filtered.rank.distance.ipynb --to python

}

data=$1
version=$2
mode=$3

#cutoff=$4

#setDir=/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Benchmark/Remove-Ambiguous-Pairs.Fixed-Ratio
setDir=/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Benchmark/All-Pairs.Fixed-Ratio
train=$setDir/$data-Benchmark.$version.txt
outputDir=/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Unsupervised_result/data
ccres=/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Benchmark/Annotations/hg19-cCREs.bed
#outputDir=~/Lab/Target-Gene/Distance-Method/Results
#ccres=~/Lab/ENCODE/Encyclopedia/V4/Registry/V4-hg19/hg19-ccREs-Simple.bed
scriptDir=/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Scripts/Unsupervised-Methods
tss=/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Benchmark/Annotations/GENCODEv19-TSSs.bed
bedtools=/usr/local/bin/bedtools
#exp=~/Lab/Target-Gene/Benchmark/Characteristics/GM12878-TPM-Expression.txt

mkdir -p $outputDir


#if [ $mode == "normal" ]
#then
python $scriptDir/rank.distance.py $tss $ccres \
    $train $outputDir/$data-Distance.$version.txt $mode
#elif [ $mode == "filtered" ]
#then
#setup
#python $scriptDir/filtered.rank.distance.py $tss $ccres \
#    $train $outputDir/$data-Filtered-Distance.$version.txt
#fi
