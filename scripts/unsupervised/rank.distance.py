import sys, numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
import pickle
def Create_TSS_Dict(tss):
    tss=open(tss)
    tssDict={}
    for line in tss:
        line=line.rstrip().split("\t")
        if line[6] in tssDict:
            tssDict[line[6]].append(int(line[1]))
        else:
            tssDict[line[6]]=[int(line[1])]
    tss.close()
    return tssDict

def Create_Enhancer_Dict(enhancers):
    enhancers=open(enhancers)
    enhancerDict={}
    for line in enhancers:
        line=line.rstrip().split("\t")
        enhancerDict[line[4]]=[int(line[1]),int(line[2])]
    enhancers.close()
    return enhancerDict


def threds_mode(filename,tssDict,enhancerDict):
    max = 200000
    threds = np.linspace(0, max, num=11, endpoint=True)
    print(threds)
    #threds=[200000]
    aucs = []
    for thred in threds:    
        y=[]
        pred=[]
        links=open(filename)
        for line in links:
            line=line.rstrip().split("\t")
            m=1000000000000
            for x in tssDict[line[1].rstrip()]:
                a=min([abs(enhancerDict[line[0].rstrip()][0]-x),abs(enhancerDict[line[0].rstrip()][1]-x)])
                if a < m:
                    m=a

            if m == 0:
                continue
            else:

                if float(m) > thred:
                    y.append(int(line[2]))
                    pred.append(1/float(m))

        y=np.array(y)
        pred=np.array(pred)
    #     print(y[:10])
    #     print(pred[:10])
    #     print('thres',thred)
        fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=1) 
    #     print(metrics.auc(fpr, tpr))
        aucs.append(metrics.auc(fpr, tpr))
        links.close()
    plt.figure()
    plt.plot(threds, aucs)
    plt.xlim([0.0, max])
    plt.ylim([0.0, 1.05])
    plt.xlabel('distance threds')
    plt.ylabel('auc')
    plt.title('CHIA-PET: auc after filtering the dataset based on the threds')
    plt.savefig('/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Unsupervised_result/figures/roc_chipet_threds.jpg')
    plt.show()
    
    
def filtered_mode(filename,tssDict,enhancerDict):
    tss_filterd=open('/Users/wangshili/Desktop/2021_spring/genomics/project/'+
                         'BENGI/Unsupervised_result/data/hic_tss_filterd.txt',"w+")
    ccre_filterd=open('/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/'+
                                      'Unsupervised_result/data/hic_ccre_filterd.txt',"w+")
    benchmark_filterd=open('/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/'+
                                                   'Unsupervised_result/data/hic_benchmark_filterd.txt',"w+")
    y=[]
    pred=[]  
    links=open(filename)
    for line in links:
        line=line.rstrip().split("\t")
        m=1000000000000
        for x in tssDict[line[1].rstrip()]:
            a=min([abs(enhancerDict[line[0].rstrip()][0]-x),abs(enhancerDict[line[0].rstrip()][1]-x)])
            if a < m:
                m=a
                
        if m == 0:
             continue
        else:
            if float(m) > 1000000.0:
                y.append(int(line[2]))
                pred.append(1/float(m))

                print(line[2]+"\t"+str(1/float(m)),file=output)

                for x in tssDict[line[1].rstrip()]:
                    print(str(x) +"\t"+line[1],file=tss_filterd)
                    print(str(enhancerDict[line[0].rstrip()][0])+"\t"+str(enhancerDict[line[0].rstrip()][1])+"\t"+str(line[0]),file=ccre_filterd)
                    print(line[0]+"\t"+line[1]+"\t"+line[2],file=benchmark_filterd)
    tss_filterd.close()
    ccre_filterd.close()
    benchmark_filterd.close()
    links.close()
    
    return np.array(y),np.array(pred)

def normal_mode(filename,tssDict,enhancerDict,output_file):
    y=[]
    pred=[] 
    links = open(filename)
    output=open(output_file, "w+")
    for line in links:
        line=line.rstrip().split("\t")
        m=1000000000000
        for x in tssDict[line[1].rstrip()]:
            a=min([abs(enhancerDict[line[0].rstrip()][0]-x),abs(enhancerDict[line[0].rstrip()][1]-x)])
            if a < m:
                m=a

        y.append(int(line[2]))
        if m == 0:
            pred.append(1)
            print(line[2]+"\t"+str(1),file=output)
        else:
            pred.append(1/float(m))
            print(line[2]+"\t"+str(1/float(m)),file=output)
    links.close()
    output.close()
    
    return np.array(y),np.array(pred)
def normal2_mode(filename,tssDict,enhancerDict):
    y=[]
    pred=[] 
    links = open(filename)

    for line in links:
        line=line.rstrip().split("\t")
        m=1000000000000
        for x in tssDict[line[1].rstrip()]:
            a=min([abs(enhancerDict[line[0].rstrip()][0]-x),abs(enhancerDict[line[0].rstrip()][1]-x)])
            if a < m:
                m=a

        y.append(int(line[2]))
        if m == 0:
            pred.append(1)
           
        else:
            pred.append(1/float(m))
            
    links.close()
    
    return np.array(y),np.array(pred)


def plot_auc(figure_name,y,pred):
    fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=1)


    print("compute area:",metrics.auc(fpr, tpr))
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.2f)' % metrics.auc(fpr, tpr))
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('The ROC curve for Hi-C data filter cCRE with close distance to TSS')
    plt.legend(loc="lower right")
    plt.savefig('/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Unsupervised_result/figures/'+figure_name+'.jpg')
    plt.show()
    
def plot2_auc(figure_name,y1,pred1,y2,pred2):
    fpr1, tpr1, _ = metrics.roc_curve(y1, pred1, pos_label=1)
    fpr2, tpr2, _ = metrics.roc_curve(y2, pred2, pos_label=1)

    print("compute area1:",metrics.auc(fpr1, tpr1))
    print("compute area2:",metrics.auc(fpr2, tpr2))
    group=['Hi-C','CHIA-PET']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cdict = {'Hi-C': 'yellow', 'CHIA-PET': 'green'}
    
    fig, ax = plt.subplots()
    
    ax.plot(fpr1, tpr1, c = 'yellow', label = 'Hi-C')
    ax.plot(fpr2, tpr2, c = 'green', label = 'CHIA-PET')
    ax.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.title('The ROC curve for Hi-C data and CHIA-PET data')
    ax.legend()
    plt.savefig('/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Unsupervised_result/figures/'+figure_name+'.jpg')
    
    
    
    
    
#     plt.figure()
#     lw = 2
# #     plt.plot(fpr, tpr, color='darkorange',
# #              lw=lw, label='ROC curve (area = %0.2f)' % metrics.auc(fpr, tpr))
#     plt.plot(fpr1, tpr1, color='darkorange')
#     plt.plot(fpr2, tpr2, color='green')
#     plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
#     plt.xlim([0.0, 1.0])
#     plt.ylim([0.0, 1.05])
#     plt.xlabel('False Positive Rate')
#     plt.ylabel('True Positive Rate')
#     plt.title('The ROC curve for Hi-C data filter cCRE with close distance to TSS')
#     plt.legend(loc="lower right")
#     plt.savefig('/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Unsupervised_result/figures/'+figure_name+'.jpg')
#     plt.show()    
    
    
    
tss=sys.argv[1]
tssDict=Create_TSS_Dict(tss)
print(tss)
enhancers=sys.argv[2]
enhancerDict=Create_Enhancer_Dict(enhancers)
print(enhancers)
filename = sys.argv[3]
#outputfile=open(sys.argv[4], "w+")
mode = sys.argv[5]
#print(filename)
#print(output)



if mode == 'threds':
    threds_mode(filename,tssDict,enhancerDict)
    
if mode == 'filtered':
    y,pred = filtered_mode(filename,tssDict,enhancerDict)
    plot_auc(figure_name,y,pred)

if mode == 'normal':
    y,pred = normal_mode(filename,tssDict,enhancerDict,sys.argv[4])
    plot_auc('filename',y,pred)
if mode == 'preprocess':
    #f = open("hictssDict.pkl","wb")
    f = open("chiatssDict.pkl","wb")
    pickle.dump(tssDict,f)
    f.close()
    
    #f = open("hicenhancerDict.pkl","wb")
    f = open("chiaenhancerDict.pkl","wb")
    pickle.dump(enhancerDict,f)
    f.close()

if mode == 'two':  
    f = open("hictssDict.pkl","rb")
    hictssDict = pickle.load(f)
    f.close()

    f = open("hicenhancerDict.pkl","rb")
    hicenhancerDict = pickle.load(f)
    f.close()
    
    f = open("chiatssDict.pkl","rb")
    chiatssDict = pickle.load(f)
    f.close()

    f = open("chiaenhancerDict.pkl","rb")
    chiaenhancerDict = pickle.load(f)
    f.close()
    
    y1,pred1 = normal_mode('/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Benchmark/Remove-Ambiguous-Pairs.Fixed-Ratio/GM12878.HiC-Benchmark.v2.txt',hictssDict,hicenhancerDict)
    y2,pred2 = normal_mode('/Users/wangshili/Desktop/2021_spring/genomics/project/BENGI/Benchmark/All-Pairs.Fixed-Ratio/GM12878.RNAPII-ChIAPET-Benchmark.v4.txt',chiatssDict,chiaenhancerDict)
    plot2_auc('2data',y1,pred1,y2,pred2)
    
# output.close()


#print('nonzero percent',np.count_nonzero(y)/len(y))


        

