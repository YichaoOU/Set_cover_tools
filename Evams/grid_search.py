from MS_CLF_PRED import MCP
import pandas as pd
import itertools as it
import numpy as np
import joblib
# MCP(data,train_x_inner,validation_x_inner,MS_program,MS_parameters,0.4,1,1)
# MCP(data,train_x_inner,validation_x_inner,MS_program,MS_parameters,RF_thres,outer_CV_count,inner_CV_count):

def parse_to_parameter_list(confDict):
	parameter_list = []
	my_dict = {}
	my_name = []
	for key in confDict['Motif_Selection']:
		if key == "program":
			continue
		my_dict[key] = confDict['Motif_Selection'][key].split()
		my_name.append(key)
	my_name.append('prob_cut')
	my_dict['prob_cut'] = confDict['prediction']['prob_cut'].split()
	# my_dict={'A':['D','E'],'B':['F','G','H'],'C':['I','J']}
	# allNames = sorted(my_dict)
	combinations = it.product(*(my_dict[Name] for Name in my_name))
	for elem in list(combinations):
		temp = {}
		count = 0
		for ele in elem:
			temp[my_name[count]] = ele
			count += 1
		parameter_list.append(temp)
	return parameter_list

def grid_search(data,skf_inner,confDict,outer_count,file_list):
	parameter_list = parse_to_parameter_list(confDict)
	
	print "parameter_list"
	print parameter_list
	scoring_matrix = pd.DataFrame(np.zeros(shape=[int(confDict['nested_CV']['inner_cv_folds'])+1, len(parameter_list)*5]))
	column_names = []
	for value in parameter_list:
		label = ".".join(map(lambda x:str(x),value.values()))
		column_names.append("Sensitivity - " + label)
		column_names.append("Specificity - " + label)
		column_names.append("Accuracy - " + label)
		column_names.append("precision - " + label)
		column_names.append("f1 - " + label)
	scoring_matrix.columns = column_names
	scoring_matrix.index = ["inner_CV_" + str(i) for i in range(1,int(confDict['nested_CV']['inner_cv_folds'])+1)] + ['mean']
	inner_count = 0
	
	for train_index_inner, validation_index_inner in skf_inner:
		inner_count += 1
		if confDict['nested_CV']['cache'] == "true":
			MCP(data,train_index_inner,validation_index_inner,confDict['Motif_Selection']['program'],parameter_list[0],outer_count,inner_count)
		my_scores = joblib.Parallel(n_jobs=int(confDict['nested_CV']['n_jobs']),verbose=int(confDict['nested_CV']['grid_search_verbosity']))(joblib.delayed(MCP)(data,train_index_inner,validation_index_inner,confDict['Motif_Selection']['program'],params,outer_count,inner_count) for params in parameter_list)
		# the score order is the same as parameter_list
		print my_scores
		# exit()
		column_index = 0
		
		for i in range(len(my_scores)):
			[se,sp,acc,precision,f1] = my_scores[i]
			
			scoring_matrix.set_value("inner_CV_" + str(inner_count),column_names[column_index],se)
			column_index += 1
			scoring_matrix.set_value("inner_CV_" + str(inner_count),column_names[column_index],sp)
			column_index += 1
			scoring_matrix.set_value("inner_CV_" + str(inner_count),column_names[column_index],acc)
			column_index += 1
			scoring_matrix.set_value("inner_CV_" + str(inner_count),column_names[column_index],precision)
			column_index += 1
			scoring_matrix.set_value("inner_CV_" + str(inner_count),column_names[column_index],f1)
			column_index += 1
			
			
	# mean_value = list(scoring_matrix.mean())
	mean_value = list(scoring_matrix.iloc[range(int(confDict["nested_CV"]["inner_cv_folds"]))].mean(axis=0))	

	scoring_matrix.loc['mean'] = mean_value
	max_value = [0,0]
	metric_number = int(confDict['nested_CV']['metric'])
	for i in range(metric_number,len(column_names),5):
		if mean_value[i] > max_value[0]:
			max_value = [mean_value[i],i]
	best_params = parameter_list[int(max_value[1]/5)]
	my_start = int(max_value[1]/5) * 5
	best_scores = [mean_value[my_start],mean_value[my_start+1],mean_value[my_start+2],mean_value[my_start+3],mean_value[my_start+4]]
	scoring_matrix_file = "outer_cv_" + str(outer_count) + ".evaluation_scores.csv"
	scoring_matrix.to_csv(scoring_matrix_file)
	file_list.append(scoring_matrix_file)
	return [best_params,best_scores]
	
def grid_search_backup(data,skf_inner,confDict,outer_count):
	parameter_list = parse_to_parameter_list(confDict)
	
	print "parameter_list"
	print parameter_list
	scoring_matrix = pd.DataFrame(np.zeros(shape=[11, len(parameter_list)*5]))
	column_names = []
	for value in parameter_list:
		label = ".".join(map(lambda x:str(x),value.values()))
		column_names.append("Sensitivity - " + label)
		column_names.append("Specificity - " + label)
		column_names.append("Accuracy - " + label)
		column_names.append("precision - " + label)
		column_names.append("f1 - " + label)
	scoring_matrix.columns = column_names
	scoring_matrix.index = ["inner_CV_" + str(i) for i in range(1,11)] + ['mean']
	inner_count = 0
	
	for train_index_inner, validation_index_inner in skf_inner:
		inner_count += 1
		
		my_scores = joblib.Parallel(n_jobs=-1,verbose=10)(joblib.delayed(MCP)(data,train_index_inner,validation_index_inner,confDict['Motif_Selection']['program'],params,outer_count,inner_count) for params in parameter_list)
		# the score order is the same as parameter_list
		print my_scores
		# exit()
		column_index = 0
		
		for i in range(len(my_scores)):
			[se,sp,acc,precision,f1] = my_scores[i]
			
			scoring_matrix.set_value("inner_CV_" + str(inner_count),column_names[column_index],se)
			column_index += 1
			scoring_matrix.set_value("inner_CV_" + str(inner_count),column_names[column_index],sp)
			column_index += 1
			scoring_matrix.set_value("inner_CV_" + str(inner_count),column_names[column_index],acc)
			column_index += 1
			scoring_matrix.set_value("inner_CV_" + str(inner_count),column_names[column_index],precision)
			column_index += 1
			scoring_matrix.set_value("inner_CV_" + str(inner_count),column_names[column_index],f1)
			column_index += 1
			
			
	mean_value = list(scoring_matrix.mean())
	scoring_matrix.loc['mean'] = mean_value
	# selected best params for f1
	max_value = [0,0]
	for i in range(4,len(column_names),5):
		if mean_value[i] > max_value[0]:
			max_value = [mean_value[i],i]
	best_params = parameter_list[int(max_value[1]/5)]
	best_scores = [mean_value[i-4],mean_value[i-3],mean_value[i-2],mean_value[i-1],mean_value[i]]
	scoring_matrix_file = "outer_cv_" + str(outer_count) + ".evaluation_scores.csv"
	scoring_matrix.to_csv(scoring_matrix_file)
	return [best_params,best_scores]
	



