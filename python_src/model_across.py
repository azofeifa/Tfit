from model import EMGU
import multiprocessing as mp

import os
def checkFileExists(FILE, i):
	if os.path.exists(FILE + str(i)):
		return checkFileExists(FILE, i+1)
	return FILE + str(i)

def wrapper_fit_function(X, k, 
max_iterations=200, convergence_thresh=0.0001,
move_uniform=0):
	clf 	= EMGU(max_ct=convergence_thresh, max_it=max_iterations, K=k, bayes=False, noise=True, 
			noise_max=0.1, moveUniformSupport=0, cores=4)
	clf.fit(X)
	return clf.ll , clf.rvs, clf.converged, clf.resets,clf

	

def run(D, bic, rounds, max_k, 
	standardize, convergence_thresh,
			max_iterations, move_uniform, 
			write_out_dir):
	#if BIC is 0 don't perform model selection and output each model 
	#from 1  to max_k individually
	FILE 	= write_out_dir+"EMG_model_fits_"
	FILE 	= checkFileExists(FILE, 0)
	FHW 	= open(FILE, "w")
	for d in D:
		d.X[:,0]-=min(d.X[:,0])
		d.X[:,0]/=standardize
		if move_uniform == 0: #lets parrallelize the rest of this
			models 	= list()
			for k in range(max_k):
				output 	= mp.Queue()
				def wrapper_fit_function_pp(X, k, output,
					max_iterations=max_iterations, convergence_thresh=convergence_thresh,
					move_uniform=move_uniform):
	
					clf 	= EMGU(max_ct=convergence_thresh, max_it=max_iterations, K=k, bayes=False, noise=True, 
							noise_max=0.1, moveUniformSupport=0, cores=4)
					clf.fit(X)
					output.put((clf.ll , clf.rvs, clf.converged, clf.resets, clf))
				processes 		= [mp.Process(target=wrapper_fit_function_pp, args=(d.X, k,output) )for i in range(rounds)]
				for p in processes:
				    p.start()
				for p in processes:
					p.join()
				keepers 	= [output.get() for p in range(rounds)]
				models.append(min(keepers))
		else:

			models 	= [max([ wrapper_fit_function(d.X, k,
				max_iterations=max_iterations,
				convergence_thresh=convergence_thresh,
				move_uniform=move_uniform ) for r in range(rounds)]) for k in range(max_k+1)]
			

		FHW.write("#"+d.print_info())
		for ll,model, converged, resets,clf in models:
			FHW.write("~"+str(len(model)) + "," + str(ll) + "," + str(converged) + "," + str(resets)+ "\n")
			model_txt = "\n".join([ m.__str__() for m in model])
			FHW.write(model_txt+"\n")
		FHW.flush()