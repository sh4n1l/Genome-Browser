import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches 
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i1', '--pslinput_file1')
parser.add_argument('-i2', '--pslinput_file2')
parser.add_argument('-g', '--gtfinput_file')
parser.add_argument('-s', '--style_sheet', default='BME163')
parser.add_argument('-o', '--output_file')
args = parser.parse_args()
plt.style.use(args.style_sheet)

def readPSL(inFile):    
    readList=[]    
    openFile=open(inFile,'r')    
    for line in openFile:        
        a=line.strip().split('\t')        
        chromosome=a[13]        
        start=int(a[15])        
        end=int(a[16])        
        blockstarts=np.array(a[20].split(',')[:-1],dtype=int)        
        blockwidths=np.array(a[18].split(',')[:-1],dtype=int)        
        read=[chromosome,start,end,blockstarts,blockwidths,False,'NA']        
        readList.append(read)    
    return readList

def readGTF(inFile):
    gtfDic = {}
    transcriptList = []
    for line in open(inFile):
        if line.startswith('#') == False:
            myString=line.strip().split('\t')
            chromosome = myString[0]
            type1 = myString[2]
            if type1 in ['exon', 'CDS']:
                start = int(myString[3])
                end = int(myString[4])
                transcript = myString[8].split(' transcript_id "')[1].split('"')[0]
                if transcript not in gtfDic:
                    gtfDic[transcript] = []
                gtfDic[transcript].append([chromosome,start,end,type1])
        
    for transcript,parts in gtfDic.items():        
        starts=[]        
        ends=[]        
        blockstarts=[]        
        blockwidths=[]        
        types = []        
        for part in parts:            
            starts.append(part[1])            
            ends.append(part[2])            
            blockstarts.append(part[1])            
            blockwidths.append(part[2]-part[1])            
            chromosome=part[0]            
            types.append(part[3])
        transcriptList.append([chromosome,min(starts),max(ends),blockstarts,blockwidths,False,types])
        
    return transcriptList
                
def plotReads(panel,readList,target):
    genome_chromosome,genome_start,genome_end=target[0],target[1],target[2]    
    bottom=1   
    filterList = []
    numplotted = 0
    for read in readList:
        chromosome,start,end,blockstarts,blockwidths,plotted,type1=read[0],read[1],read[2],read[3],read[4],read[5],read[6]       
        if chromosome==genome_chromosome:            
            if genome_start<start<genome_end or genome_start<end<genome_end:  
                filterList.append(read)
    
    for y in range(1, len(filterList)):
        lastplot =0
        for read in filterList:
            chromosome,start,end,blockstarts,blockwidths,plotted,type1=read[0],read[1],read[2],read[3],read[4],read[5],read[6]
            if plotted == False:
                if start > lastplot:
                    rectangle1=mplpatches.Rectangle((start,y),end-start, 0.1,facecolor='black',edgecolor='black',linewidth=0)
                    panel.add_patch(rectangle1) 
                    numplotted+=1
                    for index in np.arange(0,len(blockstarts),1):
                        if type1 == 'NA':
                            blockstart=blockstarts[index]                    
                            blockwidth=blockwidths[index]                    
                            rectangle1=mplpatches.Rectangle((blockstart,y-.2),blockwidth,0.5,facecolor='black',edgecolor='black',linewidth=0)                    
                            panel.add_patch(rectangle1)
                        elif type1[index] == 'CDS':
                            blockstart=blockstarts[index]                    
                            blockwidth=blockwidths[index]                    
                            rectangle1=mplpatches.Rectangle((blockstart,y-.2),blockwidth,0.5,facecolor='black',edgecolor='black',linewidth=0)                    
                            panel.add_patch(rectangle1)
                        elif type1[index] == 'exon':
                            blockstart=blockstarts[index]                    
                            blockwidth=blockwidths[index]                    
                            rectangle1=mplpatches.Rectangle((blockstart,y-.075),blockwidth,0.25,facecolor='black',edgecolor='black',linewidth=0)                    
                            panel.add_patch(rectangle1)
                    read[5]=True
                    lastplot = end
        if len(filterList) == numplotted:
            break
        
    return y

chromosome = 'chr' + input('Chromosome #?')
start = int(input("Start Location?"))
end = int(input("End Location?"))
#target=['chr7',45232945,45240000]
target=[chromosome,start,end]
figWidth = 10
figHeight = 5

plt.figure(figsize=(figWidth, figHeight))

panWidth = 10/figWidth
panHeight = 1.25/figHeight

myPanel3 = plt.axes([(0), (.05), panWidth, panHeight])
myPanel2 = plt.axes([(0), (.35), panWidth, panHeight])
myPanel1 = plt.axes([(0), (.65), panWidth, panHeight])

readList2 = readPSL(args.pslinput_file2)
sortList2 = sorted(readList2, key = lambda x: x[1])
ylim2=plotReads(myPanel2,sortList2,target)
myPanel2.set_xlim(target[1],target[2])
myPanel2.set_ylim(0,ylim2+7)

readList3 = readPSL(args.pslinput_file1)
sortList3 = sorted(readList3, key = lambda x: x[2])
ylim3=plotReads(myPanel3,sortList3,target)
myPanel3.set_xlim(target[1],target[2])
myPanel3.set_ylim(0,ylim3+40)

transcriptList = readGTF(args.gtfinput_file)
sortList1 = sorted(transcriptList, key = lambda x: x[2])
ylim1=plotReads(myPanel1,sortList1,target)
myPanel1.set_xlim(target[1],target[2])
myPanel1.set_ylim(0,ylim1+1.9)


myPanel1.tick_params(bottom=False, labelbottom=False,
                left=False, labelleft=False,                   
                right=False, labelright=False,                   
                top=False, labeltop=False,)
myPanel2.tick_params(bottom=False, labelbottom=False,
                left=False, labelleft=False,                   
                right=False, labelright=False,                   
                top=False, labeltop=False,)
myPanel3.tick_params(bottom=False, labelbottom=False,
                left=False, labelleft=False,                   
                right=False, labelright=False,                   
                top=False, labeltop=False,)

plt.savefig(args.output_file, dpi=1200)
