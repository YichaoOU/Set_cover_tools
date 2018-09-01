from __future__ import division
import pandas as pd

# python MS_greedy.py stress_responsive_genes_test_dataset.csv test.list tomtom.txt 0.02


import uuid
import os
import sys
import argparse
import shutil
import re
from Emotif.algo import greedy

from Emotif.utils import *
from Emotif.utils import Tomtom
from Emotif.motif_filtering import *
from Emotif.motif_output import *




def outputResults(outFileName, depthDict, filterMinNumSeq, idMotifDict, foreNumSeqs):

	
	
	# depthOut = open(outFileName, 'wb')
	ListOut = open(outFileName, 'wb')
	# depthOut.write('#Name,foreground_motif_seqCount,foreground_motif_seqCoverage(%),foreground_seqs_added,cumulative_foreground_seq_cov(%)\n')
	cumCov = 0
	motifIdList = depthDict[1][0]
	motifSeqSet = depthDict[1][1]
	newSeqsAddedList = depthDict[1][2]
	#add new information to the motif objects as well
	selected_motifs = []
	for i in range(len(motifIdList)):
		lineList = []
		motifId = motifIdList[i]
		motifSeqset = motifSeqSet[i]
		seqsAdded = newSeqsAddedList[i]
		#apply the filter threshold
		if seqsAdded < filterMinNumSeq:
			continue
		if motifId not in idMotifDict:
			motifName = motifId
		else:
			motifName = idMotifDict[motifId]
		print >>ListOut,motifName
		lineList.append(motifName)
		selected_motifs.append(motifName)
		lineList.append(str(len(motifSeqset)))
		cov = 100 * (len(motifSeqset)/foreNumSeqs)
		lineList.append(str(cov))
		lineList.append(str(seqsAdded))
		cumCov += (seqsAdded/foreNumSeqs)*100
		lineList.append(str(cumCov))
		line = ','.join(lineList)
		# depthOut.write(line + '\n')
	# depthOut.close()
	return selected_motifs


def csv_to_mhit(csv_file):
	mhit_file = str(uuid.uuid1())
	out_file = open(mhit_file,"wb")
	motif_name = []
	mhit_dict = {}
	flag = True
	count = 0 
	foreNumSeqs=0
	with open(csv_file) as f:
		for line in f:
			if flag:
				motif_name = line.strip().split(",")[1:-1]
				for m in motif_name:
					mhit_dict[m] = []
				flag = False
				continue
			
			label = line.strip().split(",")[-1]
			# print line,label
			seq_id = line.split(",")[0]
			if float(label) == 1:
				count += 1
				# print >>out_file,">" + line.split(",")[0]
				my_index = 0
				for m in line.split(",")[1:-1]:
					if float(m) == 1:
						# print >>out_file,motif_name[my_index]
						mhit_dict[motif_name[my_index]].append(seq_id)
					my_index+= 1
	foreNumSeqs = count
	for m in mhit_dict:
		print >>out_file,">" + m
		print >>out_file,"\n".join(mhit_dict[m])
	out_file.close()
	return [mhit_file,foreNumSeqs]
	










# input
[input,foreNumSeqs] = csv_to_mhit(sys.argv[1])
output = sys.argv[2]
tomtom_file = sys.argv[3]
threshold = float(sys.argv[4])


# import uuid


# set cover
tomtomDict = Tomtom.parse(tomtom_file)
motifIdDict, idMotifDict, seqIdDict, idSeqDict, Uset, Sdict = general_utils.processMotifHitFile_1(input)
filterMinNumSeq = int(foreNumSeqs * threshold)
depthDict = greedy.callGreedyDepthRefine(Uset, Sdict, 1, tomtomDict, motifIdDict, idMotifDict, filterMinNumSeq)




# output
outputResults(output, depthDict, filterMinNumSeq, idMotifDict, foreNumSeqs)
os.remove(input)



























