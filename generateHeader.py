def readFile(path):
    f = open(path, 'r')
    Lines = f.readlines()
    predicts = []
    labels = []
    f1 = open("header.txt", 'w')
    for line in Lines:
        if line.startswith("MOTIF"):
            f1.write(line.strip()[6:]+"\n")
    f.close()
    f1.close()
            


readFile("/Users/xuxi/Documents/cg/BENGI/Scripts/Supervised-Methods/PEP/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme")