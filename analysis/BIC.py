import math 
import numpy as np
import matplotlib.pyplot as plt
import time, load, merge_data_types as mdt
import sys



def bic_function(ll, n, K, penality):
	return -2.*ll + K*4*math.log(n)*penality
def get_best_model(I, penality , diff_threshold):
	models 	 	= dict([ (k, MAX([(model.ll, model) for model in I.models[k] if model.diff < diff_threshold ])[1]) for k in I.models])
	BIC_best 	= min([ (bic_function(models[k].ll, I.N, k*9, penality), k ) if models[k] is not None else (np.inf, k) for k in models  ])[1]
	return models[BIC_best]


def MAX(LST):
	if not LST:
		return None, None
	return max(LST)

def run(G, figName):
	data 	= list()
	for I in G:
		if I.annotation_N==1:
			models 	= dict([ (k, MAX([(model.ll, model) for model in I.models[k] if model.diff < 10 ])[1]) for k in I.models])
			data.append((models, I))
	penality 	= np.linspace(0.1, 500, 100)
	scatter 	= list()
	for p in penality:
		error=list()
		for models,I in data:
			BIC_best 	= min([ (bic_function(models[k].ll, I.N, k*9, p), k ) if models[k] is not None else (np.inf, k) for k in models  ])[1]

			error.append(BIC_best-1)
		scatter.append(error)
	F 		= plt.figure(figsize=(15,10))
	ax1 	= F.add_subplot(211)
	ax1.scatter(penality, [np.mean(s) for s in scatter] )
	ax1.fill_between(penality, [np.mean(s)- np.std(s) for s in scatter], [np.mean(s) + np.std(s) for s in scatter] , color="grey", alpha=0.5 )
	ax1.grid()
	ax1.set_xlabel("BIC Penality")
	ax1.set_ylabel("Error (Mean)")
	ax2 	= F.add_subplot(212)
	
	ax2.scatter(penality, [np.sum(s) for s in scatter] )
	ax2.grid()
	ax2.set_xlabel("BIC Penality")
	ax2.set_ylabel("Error (sum)")
	plt.savefig(figName)


if __name__=="__main__":
	if len(sys.argv)==1:
		merged_file 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/merged_data_file_100.txt"
		figName 		= "/Users/joeyazo/Desktop/BIC"
	else:
		merged_file 	= sys.argv[1]
		figName 		= sys.argv[2]
	

	G 				= load.merge_data_out(merged_file)
	run(G, figName)