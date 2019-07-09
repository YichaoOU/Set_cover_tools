from __future__ import division
# python MSDC_wrapper.py fore.list back.list motif.mapping .k .j
import uuid
import os
import sys
import shutil
import re

def csv_to_input_files(csv_file):
	fore_list = "fore"+str(uuid.uuid1())
	back_list = "back"+str(uuid.uuid1())
	motif_mapping = str(uuid.uuid1())
	fore_list_file = open(fore_list,"wb")
	back_list_file = open(back_list,"wb")
	motif_mapping_file = open(motif_mapping,"wb")
	print >>motif_mapping_file,"\t".join(["motif","seq"])
	motif_name = []
	flag = True
	with open(csv_file) as f:
		for line in f:
			if flag:
				motif_name = line.strip().split(",")[1:-1]
				flag = False
				continue
			label = line.strip().split(",")[-1]
			seq_id = line.split(",")[0]
			if float(label) == 1:
				print >>fore_list_file,seq_id
				my_index = 0
				for m in line.split(",")[1:-1]:
					if float(m) == 1:
						print >>motif_mapping_file,"\t".join([motif_name[my_index],seq_id])
					my_index+= 1
			if float(label) == -1:
				print >>back_list_file,seq_id
				my_index = 0
				for m in line.split(",")[1:-1]:
					if float(m) == 1:
						print >>motif_mapping_file,"\t".join([motif_name[my_index],seq_id])
					my_index+= 1
	fore_list_file.close()
	back_list_file.close()
	motif_mapping_file.close()
	return [fore_list,back_list,motif_mapping]
	

def parse_output(temp_output,output):
	out = open(output,"wb")
	lines = open(temp_output).readlines()
	if lines[1].split()[0] == "Failed":
		out.close()
		return 1
	print >>out,"".join(lines[1:-1])
	return 1



# Failed



# input
[fore_list,back_list,motif_mapping] = csv_to_input_files(sys.argv[1])
output = sys.argv[2]
temp_output = str(uuid.uuid1())
k = sys.argv[3]
j = sys.argv[4]

command = "msdc " + fore_list + " " + back_list + " " + motif_mapping + " ." + k + " ." + j + " > " + temp_output
os.system(command)
parse_output(temp_output,output)


os.remove(fore_list)
os.remove(back_list)
os.remove(motif_mapping)
os.remove(temp_output)



























































