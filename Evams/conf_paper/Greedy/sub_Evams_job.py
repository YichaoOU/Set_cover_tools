#!/usr/bin/env python
import subprocess 

import glob
import os

with open("/users/PHS0293/ohu0404/project/motif_selection_paper/Motifs/TF_groups_with_known_motifs.list") as f:
	for line in f:
		print line
		line = line.strip()
		if len(line) <= 0:
			continue
		input = "/users/PHS0293/ohu0404/project/motif_selection_paper/Motif_Selection/" + line + ".Motif_Selection_Paper.csv"
		try:
			open(input)
		except:
			print line,"is not ready"
			continue
		folder = line + "_Evams_result"
		os.makedirs(folder)
		lines= open("/users/PHS0293/ohu0404/project/motif_selection_paper/2017_10_3_new_result/Greedy/TF_group.Evams.conf").readlines()
		lines[1] = lines[1].replace("TF_group1",line)
		lines[27] = lines[27].replace("TF_group1",line)
		conf = line+".Evams.conf"
		f = open(conf,"wb")
		print >>f,"".join(lines)
		f.close()
		name = line+".job"
		qsub_command = "qsub -v folder={0},TF_group={1},conf={2},job_name={3} Evams_job".format(folder,line,conf,name)
		exit_status = subprocess.call(qsub_command, shell=True)
		if exit_status is 1:  # Check to make sure the job submitted
			print file,": is failed to submit"
			print qsub_command
		

















