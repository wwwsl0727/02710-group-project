#!/bin/bash
setup() {
  jupyter nbconvert analyze.distance.ipynb --to python

}
setup

data=$1
version=$2

inputDir=/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Unsupervised_learning
scriptDir=/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Scripts/Unsupervised-Methods
outputDir=inputDir


python $scriptDir/analyze.distance.py $inputDir/$data-Distance.$version.txt
#python $scriptDir/thredshold.py $inputDir/$data-Filtered-Distance.$version.txt $outputDir

#python $scriptDir/distance.acc.py $outputDir
