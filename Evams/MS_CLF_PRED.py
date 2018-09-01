import subprocess
from sklearn.ensemble import RandomForestClassifier as RFC
import os
import numpy as np
def MCP(data,train_x_inner,validation_x_inner,MS_program,params,outer_CV_count,inner_CV_count):
	MS_parameters = {}
	for key in params:
		if not key == "prob_cut":
			MS_parameters[key] = params[key]
	RF_thres = float(params['prob_cut']) 
	prefix = ".".join(map(lambda x:str(x),[outer_CV_count,inner_CV_count]+MS_parameters.values()+[RF_thres]))
	MS_input = data.iloc[list(train_x_inner)]
	MS_input_file =  prefix + ".MS_input.csv"
	MS_output_file =  prefix + ".MS_output.list"
	MS_input.to_csv(MS_input_file)
	MS_parameters['input'] = MS_input_file
	MS_parameters['output'] = MS_output_file
	selected_motifs = Motif_Selection(MS_program,MS_parameters)
	print "selected_motifs",selected_motifs
	if len(selected_motifs) == 0:
		os.remove(MS_input_file)
		os.remove(MS_output_file)
		return [0]*5
	# train = data.iloc[list(train_x_inner)]
	# Clf_X = train[selected_motifs]
	# Clf_Y = train[train.columns[-1]]
	# print Clf_Y
	# print "Clf_Y"
	# clf = Classification(Clf_X,Clf_Y)
	valid = data.iloc[list(validation_x_inner)]
	valid_X = valid[selected_motifs]
	valid_Y = valid[valid.columns[-1]]
	# print valid_Y
	# print "valid_Y"
	[se,sp,acc,precision,f1] = Prediction(valid_X,valid_Y)
	# print [se,sp,acc,precision,f1]
	os.remove(MS_input_file)
	os.remove(MS_output_file)
	return [se,sp,acc,precision,f1]
	
def multipleReplace(text, wordDict):
	for key in wordDict:
		text = text.replace(key, wordDict[key])
	return text
	
def parse_list(file):
	temp = []
	with open(file) as f:
		for line in f:
			try:
				line = line.strip().split()[0]
				# print line
			except:
				continue
			if len(line) > 0:
				temp.append(line)
	return temp
	
def Motif_Selection(program,parameters):
	# program is string
	# parameters is dict
	# return a list of selected motifs
	command = multipleReplace(program,parameters)
	print command
	# subprocess.Popen(command,stdout=subprocess.PIPE)
	os.system(command)
	
	return parse_list(parameters['output'])
	
	
def prob_to_label(pred,thres):
	temp = []
	for i in pred:
		if i[1] < thres:
			temp.append(-1)
		else:
			temp.append(1)
	return temp

def Classification(X,Y):
	clf=RFC(n_estimators=100)
	clf.fit(X, Y)
	return clf
	
	
def calculate_SE_SP_ACC(y_pred,y_true,wrt=1):
	# diabetes are 1
	#print y_pred,y_true
	tp = 0.0
	tn = 0.0
	correct = 0.0
	if len(y_pred) != len(y_true):
		print "len(y_pred) != len(y_true)!"
		exit()
	for i in range(len(y_pred)):
		if y_pred[i] == y_true[i]:
			if y_pred[i] == wrt:
				tp += 1
			else:
				tn += 1
			correct += 1
	total = float(len(y_true))
	T = y_true.count(wrt)
	N = total - T
	#print T,N,total
	if T == 0:
		
		return [1,tn/N,correct/total,0,tn,0,0]
	if N == 0:
		return [tp/T,1,correct/total,tp,0,1,1]
	#print total
	fp = N - tn
	se = tp/T
	try:
		precision = tp/(tp+fp)
	except:
		precision = 1
	f1 = (2*precision*se)/(precision+se)
	
	return [se,tn/N,correct/total,tp,tn,precision,f1]

from sklearn.metrics import auc,precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import recall_score
def boolean_mapping(x):
	if x:
		return 1
	else:
		return -1

def Prediction(X,Y):
	y_true = Y
	# y_pred_prob = clf.predict_proba(X)
	y_pred = map(lambda x:boolean_mapping(x),X.sum(axis=1)>0)
	acc=accuracy_score(y_true, y_pred)
	[se,sp,acc1,tp,tn,precision,f1] = calculate_SE_SP_ACC(map(lambda x:int(x),list(y_pred)),map(lambda x:int(x),list(y_true)))
	if acc != acc1:
		print "ALERT! SOMETHING WRONG in ACC"
		print acc,acc1
	return [se,sp,acc,precision,f1]
	
	




















