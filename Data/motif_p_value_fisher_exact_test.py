
import glob

import pandas as pd
from scipy.stats import chi2_contingency
import numpy as np

def output_p_value_list(file):
	df = pd.read_csv(file,index_col=0)
	output_dict = {}
	for m in df.columns.tolist()[:-1]:
		A = df[df['class']==1][m].sum()
		B = df[df['class']==1].shape[0] - A
		C = df[df['class']==-1][m].sum()
		D = df[df['class']==-1].shape[0] - C
		print (m,[[A,B],[C,D]])
		if A == 0:
			output_dict[m]=np.nan
		else:
			odd,p_value,t1,t2 = chi2_contingency([[A,B],[C,D]])
			output_dict[m]=p_value
	
	df = pd.DataFrame.from_dict(output_dict, orient='index')
	df.columns = ['p-value']
	df.sort_values('p-value',ascending=True,inplace=True)
	df.to_csv(file.replace(".csv.gz",'.p_value.csv'))
	return output_dict
	

for i in glob.glob("*.gz"):
	print (i)
	output_p_value_list(i)









