import os
import subprocess	
from decimal import *

threshold = 0
nSpikes = 0
spikeBegin = False
spikeEnds = False
vmdata = []
dendMax = []
headMax = []
#opFile = []
voltage = -0.06

#morphPath = '/home/sarathy/Work/sample/modelNelec/cells'
#for i in os.listdir(morphPath):
    #if i.endswith('.swc'):
        #opFile.append(i.split('.')[0]+'.txt')

#Compute firing rate
for kk in moose.wildcardFind( "/graphs/VmTab[0]" ):
	vmdata.append(kk.vector)
	
for xx in vmdata[0]:
	if xx > threshold:
		if not spikeBegin:
			spikeBegin = True
			spikeEnds = False
		else: pass
	else:
		if spikeBegin:
			spikeEnds = True
			spikeBegin = False
			nSpikes += 1
			
spikeStartIndex = numpy.where(vmdata[0] > voltage)[0][0]			
spikeEndIndex = numpy.where(vmdata[0] > voltage)[0][-1]	
			
Time = numpy.arange( 0, len(vmdata[0]), 1 ) * 0.0001

spikeStartTime = Time[spikeStartIndex]
spikeEndTime = Time[spikeEndIndex]
simuTime = spikeEndTime - spikeStartTime
	
firingRate = str(nSpikes/simuTime)

for ii  in moose.wildcardFind( "/graphs/caDendTab[]" ):
    dendMax.append(numpy.max(ii.vector))	#collecting only the max of each dendrite
pylab.hist(dendMax)
pylab.show()    

meandendMax = numpy.mean(dendMax)
stddendMax = numpy.std(dendMax)

for jj in moose.wildcardFind( "/graphs/caHeadTab[]" ):
    headMax.append(numpy.max(jj.vector))	#collecting only the max of each spine head
pylab.hist(headMax)
pylab.show()  
print '-------------'
print len(dendMax)
print len(headMax)
meanheadMax = numpy.mean(headMax)
stdheadMax = numpy.std(headMax)

data = [nSpikes, firingRate, meandendMax, stddendMax, meanheadMax, stdheadMax]
#for fil in opFile:	
fp = open(opFile, "w")
fp.write('#Spike,FiringRate,[Ca]Dend_mean,[Ca]Dend_std,[Ca]Head_mean,[Ca]Head_std\n')
for mm, ll in enumerate(data):
	if mm==len(data)-1:
		fp.write(str(ll))
		fp.write('\n')
	else:
		fp.write(str(ll))
		fp.write(',')
fp.close()

