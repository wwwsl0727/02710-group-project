import xgboost, joblib
from operator import add
import random, subprocess,numpy, sys, sklearn, scipy, itertools, shutil, math
        
def Run_XGB(allLabels, allFeatures, valLabels, valFeatures, output,\
                      repetitions, i, header,outputDir,version, name):
    valOutput=output+"-validation-"+str(i)+"."+version+".txt"
    fiOutput=output+"-FI-"+str(i)+"."+version+".txt"
    output=open(valOutput, "w")
    output2=open(fiOutput, "w")
    oob=[]
    acc=[]
    o=[0]*len(valLabels)
    FI=[0]*len(valFeatures[0])
    y=numpy.array(allLabels)

    weight = float(len(y[y == 0]))/float(len(y[y == 1]))
    w1 = [1]*y.shape[0]
    k=0
    for entry in allLabels:
        if entry == 1:
            w1[k] = weight
        k+=1
    w1=numpy.array(w1)

    for x in range(0,repetitions):        
        xgb = xgboost.XGBClassifier(max_depth=10, learning_rate=0.1, \
                                    n_estimators=1000, nthread=50)
        xgb.fit(numpy.array(allFeatures), numpy.array(allLabels), sample_weight=w1)
        joblib.dump(xgb, outputDir+'/'+name+'.PEP.model.'+str(x)+'.'+str(i)+'.pkl')
        M=xgb.predict_proba(numpy.array(valFeatures))
        predictions=xgb.predict(numpy.array(valFeatures))

        i=0
        correct=0
        incorrect=0
        for x in predictions:
            if x == valLabels[i]:
                correct += 1
            else:
                incorrect += 1
            i += 1
        acc.append(correct/float(correct+incorrect))
        k=0
        for entry in M:
            o[k] += entry[1]
            k += 1
        FI=[x + y for x, y in zip(FI, xgb.feature_importances_)]
    print numpy.mean(acc), weight
    j=0
    for element in FI:
        print >> output2, header[j], "\t", element/float(repetitions)
        j+=1
    print "\n"
    i=0
    for entry in o:
        print >> output, valLabels[i], "\t", entry/float(repetitions)
        i += 1
    output.close()
    output2.close()
    shutil.move(valOutput, outputDir+"/"+valOutput)
    shutil.move(fiOutput, outputDir+"/"+fiOutput)

def Run_Model(trainFeat,trainLab, outputPrefix, cvList, cvGroups, header,\
    outputDir, version):
    i=1
    for group in cvList:
        print "Running cross-validation for "+ group
        trainMatrix=[]
        testMatrix=[]
        trainY=[]
        testY=[]
        j=0
        for entry in cvGroups:
            if entry == group:
                testMatrix.append(trainFeat[j])
                testY.append(trainLab[j])
            else:
                trainMatrix.append(trainFeat[j])
                trainY.append(trainLab[j])
            j+=1
        Run_XGB(trainY, trainMatrix, testY, testMatrix, outputPrefix, 1, \
            group, header, outputDir, version, outputPrefix)
        i+=1    
    return


def Create_Gene_Dict(tss):
    tssDict={}
    geneDict={}
    for line in tss:
        line=line.rstrip().split("\t")
        tssDict[line[3]]=line[4]
        if line[4] not in geneDict:
            geneDict[line[4]]=[line[3]]
        else:
            geneDict[line[4]].append(line[3])
    return tssDict, geneDict
    
def Process_ELS_Gene_Pairs(pairs):
    pairArray=[]
    cvList=[]
    for line in pairs:
        line=line.rstrip().split("\t")
        pairArray.append([line[0],line[1],int(line[2]),line[3]])
        if line[3] not in cvList:
            cvList.append(line[3])
    return pairArray, cvList

def Process_Peak_Matrix(matrix, mode, header):
    elementDict={}
    h=matrix.next().rstrip().split("\t")[1:]
    for entry in h:
        header.append(entry+"-"+mode)
    for line in matrix:
        line=line.rstrip().split("\t")
        elementDict[line[0]]=[float(i) for i in line[1:]]
    return elementDict, header
#geneDict key ENSG00000237738.1
#tssDict key ENST00000395864.3
def Create_Feature_Array(data, enhancerSignals, tssSignals, geneDict):
    labels=[]
    cvGroups=[]
    features=[]
    for pair in data:
        tssFeatures=[]
        flag = False
        for tss in geneDict[pair[1]]:
            if tss in tssSignals.keys():
                flag = True
            else:
                continue
            if len(tssFeatures) > 0:
                tssFeatures=map(add, tssFeatures, tssSignals[tss])
            else:
                tssFeatures=tssSignals[tss]
        if not flag:
            continue
        tssFeatures=[x / float(len(geneDict[pair[1]])) for x in tssFeatures]
        if pair[0] not in enhancerSignals.keys():
            continue
        features.append(enhancerSignals[pair[0]]+tssFeatures)
        labels.append(pair[2])
        cvGroups.append(pair[3])
    return features, labels, cvGroups

header=[]
trainingPairs=open(sys.argv[1])
trainingArray, cvList=Process_ELS_Gene_Pairs(trainingPairs)
trainingPairs.close()

enhancerMatrix=open(sys.argv[2])
enhancerSignals, header=Process_Peak_Matrix(enhancerMatrix, "enhancer",header)
#print(enhancerSignal['EH37E0111085'])
enhancerMatrix.close()

tssMatrix=open(sys.argv[3])
tssSignals, header=Process_Peak_Matrix(tssMatrix,"tss",header)
#print(tssSignals['ENST00000395864.3'])
tssMatrix.close()

tss=open(sys.argv[4])
tssDict, geneDict=Create_Gene_Dict(tss)
tss.close()

outputPrefix=sys.argv[5]
outputDir=sys.argv[6]
version=sys.argv[7]

print "processing feature matrices..."
trainFeat, trainLab, cvGroups, =Create_Feature_Array(trainingArray, \
    enhancerSignals, tssSignals, geneDict)
with open('trainFeat.txt', 'w') as filehandle1:
    filehandle1.writelines("%s\n" % tf for tf in trainFeat)

print "finish trainFeat"

with open('trainLab.txt', 'w') as filehandle2:
    filehandle2.writelines("%s\n" % label for label in trainLab)

print "finish label"

with open('cvGroups.txt', 'w') as filehandle3:
    filehandle3.writelines("%s\n" % group for group in cvGroups)

print "finish label"

print "finish trainFeat"

print "running models..."
'''
Run_Model(trainFeat, trainLab, outputPrefix, cvList, cvGroups, header, \
    outputDir, version)
'''
#python2 xgb.py /Users/xuxi/Documents/cg/BENGI/Benchmark/data/GM12878.HiC-Benchmark.v2.txt Enhancer-Feature-Matrix.txt tss-Feature-Matrix.txt tss GM12878.CHiC content/result v2

def reRun_Model(allLabels, allFeatures, valLabels, valFeatures, output, repetitions, i, header,outputDir,version, name)):
    path = "result/GM12878.CHiC.PEP.model.0.cv-"+str(i)+".pkl"
    xbg = joblib.load(path)
    predictions = gb_clf.predict(numpy.array(valFeatures))
    print("Classification Report")
    matrix = classification_report(y_val, predictions)
    print(matrix)
    auc = metrics.auc(matrix[0][recall], matrix[0][precision])
    print("auc for label 0")
    print(auc)
    print("auc for label 1")
    auc = metrics.auc(matrix[1][recall], matrix[1][precision])
    print(auc)
#auc for label 0
'''
Running cross-validation for cv-8
auc
0.1528135081539404
Running cross-validation for cv-1
auc
0.23054173460568497
Running cross-validation for cv-2
auc
0.1863796472846313
Running cross-validation for cv-7
auc
0.22108369908465314
Running cross-validation for cv-11
auc
0.2809053365837533
Running cross-validation for cv-4
auc
0.1902865740040688
Running cross-validation for cv-10
auc
0.18369012844620086
Running cross-validation for cv-0
auc
0.21773044587812285
Running cross-validation for cv-5
auc
0.18445869285304836
Running cross-validation for cv-3
auc
0.22727737897406425
Running cross-validation for cv-9
auc
0.19964667492571064
Running cross-validation for cv-6
auc
0.20517203171941512
'''
#auc for label1
'''
Running cross-validation for cv-8
auc
0.26869602181922697
Running cross-validation for cv-1
auc
0.6258281494526782
Running cross-validation for cv-2
auc
0.35063082140727275
Running cross-validation for cv-7
auc
0.3915276687667354
Running cross-validation for cv-11
auc
0.29727736532699567
Running cross-validation for cv-4
auc
0.32353422705459556
Running cross-validation for cv-10
auc
0.5296290043418359
Running cross-validation for cv-0
auc
0.7188039181892771
Running cross-validation for cv-5
auc
0.383982967712955
Running cross-validation for cv-3
auc
0.5847379815099122
Running cross-validation for cv-9
auc
0.3463025057412883
Running cross-validation for cv-6
auc
0.41958873491467286
'''

#accuracy
'''
Running cross-validation for cv-8
Classification Report
              precision    recall  f1-score   support

           0       0.81      0.95      0.87       150
           1       0.30      0.08      0.13        37

   micro avg       0.78      0.78      0.78       187
   macro avg       0.55      0.52      0.50       187
weighted avg       0.71      0.78      0.73       187

Running cross-validation for cv-1
Classification Report
              precision    recall  f1-score   support

           0       0.74      0.97      0.84        67
           1       0.78      0.23      0.36        30

   micro avg       0.74      0.74      0.74        97
   macro avg       0.76      0.60      0.60        97
weighted avg       0.75      0.74      0.69        97

Running cross-validation for cv-2
Classification Report
              precision    recall  f1-score   support

           0       0.78      0.99      0.87       155
           1       0.60      0.06      0.12        47

   micro avg       0.77      0.77      0.77       202
   macro avg       0.69      0.53      0.49       202
weighted avg       0.74      0.77      0.69       202

Running cross-validation for cv-7
Classification Report
              precision    recall  f1-score   support

           0       0.71      0.97      0.82       117
           1       0.33      0.04      0.07        49

   micro avg       0.69      0.69      0.69       166
   macro avg       0.52      0.50      0.44       166
weighted avg       0.60      0.69      0.60       166

Running cross-validation for cv-11
Classification Report
              precision    recall  f1-score   support

           0       0.74      0.98      0.85       114
           1       0.60      0.07      0.13        42

   micro avg       0.74      0.74      0.74       156
   macro avg       0.67      0.53      0.49       156
weighted avg       0.70      0.74      0.65       156

Running cross-validation for cv-4
Classification Report
              precision    recall  f1-score   support

           0       0.75      0.97      0.85        71
           1       0.33      0.04      0.07        24

   micro avg       0.74      0.74      0.74        95
   macro avg       0.54      0.51      0.46        95
weighted avg       0.64      0.74      0.65        95

Running cross-validation for cv-10
Classification Report
              precision    recall  f1-score   support

           0       0.72      0.99      0.84        72
           1       0.50      0.04      0.07        28

   micro avg       0.72      0.72      0.72       100
   macro avg       0.61      0.51      0.45       100
weighted avg       0.66      0.72      0.62       100

Running cross-validation for cv-0
Classification Report
              precision    recall  f1-score   support

           0       0.79      1.00      0.88        76
           1       1.00      0.31      0.47        29

   micro avg       0.81      0.81      0.81       105
   macro avg       0.90      0.66      0.68       105
weighted avg       0.85      0.81      0.77       105

Running cross-validation for cv-5
Classification Report
              precision    recall  f1-score   support

           0       0.76      0.99      0.86       142
           1       0.67      0.04      0.08        46

   micro avg       0.76      0.76      0.76       188
   macro avg       0.71      0.52      0.47       188
weighted avg       0.74      0.76      0.67       188

Running cross-validation for cv-3
Classification Report
              precision    recall  f1-score   support

           0       0.69      0.99      0.81        80
           1       0.83      0.12      0.21        41

   micro avg       0.69      0.69      0.69       121
   macro avg       0.76      0.55      0.51       121
weighted avg       0.74      0.69      0.61       121

Running cross-validation for cv-9
Classification Report
              precision    recall  f1-score   support

           0       0.76      0.96      0.85       140
           1       0.44      0.09      0.15        46

   micro avg       0.75      0.75      0.75       186
   macro avg       0.60      0.53      0.50       186
weighted avg       0.68      0.75      0.68       186

Running cross-validation for cv-6
Classification Report
              precision    recall  f1-score   support

           0       0.72      0.94      0.82       103
           1       0.40      0.10      0.15        42

   micro avg       0.70      0.70      0.70       145
   macro avg       0.56      0.52      0.48       145
weighted avg       0.63      0.70      0.62       145

'''