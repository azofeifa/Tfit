import merge_data_types as mdt
import BIC
import matplotlib.pyplot as plt
import numpy as np
def display(G):
	for I in G:
		F 	= plt.figure(figsize=(15,10))
		

		X_f 	= np.array([x for x,y in getattr(I, "gro_f") ])
		Y_f 	= np.array([y for x,y in getattr(I, "gro_f") ])
		X_r 	= np.array([x for x,y in getattr(I, "gro_r") ])
		Y_r 	= np.array([y for x,y in getattr(I, "gro_r") ])
			
		pi 		= sum(Y_f)/(sum(Y_f) + sum(Y_r))
		Y_f, X_f= np.histogram(X_f, weights=Y_f, normed=1,bins=100)
		Y_f*=pi
		Y_r, X_r= np.histogram(X_r, weights=Y_r, normed=1,bins=100)
		Y_r*=pi
		X_f 	= (X_f[:-1] + X_f[1:]) / 2.
		X_r 	= (X_r[:-1] + X_r[1:]) / 2.

		

		minX, maxX 	= min(min(X_f), min(X_r)), max(max(X_f), max(X_r))
		width 		= (maxX-minX) / len(X_f)
		xs 			= np.linspace(minX, maxX,1000)
		
		RV 			= max([(model.ll, model) for model in I.models[1] if  model.diff < 10. and model.rvs[0].si < 5] )[1]
		
		N 			= sum(Y_f) + sum(Y_r)
		plt.bar(X_f, Y_f,color="blue", width=width)
		plt.bar(X_r, -Y_r,color="red", width=width)
		
		ys 			= map(lambda x : RV.pdf(x,1), xs)
		ys_r 		= map(lambda x : -RV.pdf(x,-1), xs)
		plt.plot(xs,ys, linewidth=1.)
		plt.plot(xs,ys_r, linewidth=1.)
		
		plt.show()
def run(EMG_IN,EMG_OUT):
	with open(EMG_OUT)  as FH:
		I 	= None
		h 	= ""
		G 	= {}
		for line in FH:
			if "#" == line[0]:
				if I is not None:
					G[h] 	= I
				chrom,info 			= line[1:].strip("\n").split(":")
				start_stop, N  		= info.split(",")
				start,stop 			= start_stop.split("-")
				h 					= "#" + chrom + "," + str(start) + "," + str(stop)
				I 					= mdt.segment(chrom,int(start),int(stop),float(N))
			elif "~" == line[0]:
				I.insert_model_info(line)
			elif "N:"==line[:2] or "U:"==line[:2]:
				I.insert_component(line)
		if I is not None:
			G[h] 	= I
			
	with open(EMG_IN) as FH:
		I 	= None
		for line in FH:
			if "#" == line[0]:
				chrom,start, stop 			= line[1:].strip("\n").split(",")
				start, stop 				= int(float(start)), int(float(stop))
				h 							= "#" + chrom + "," + str(start) + "," + str(stop)
				if h in G:
					I 							= G[h]
				else:
					I 						= None
			elif I is not None:
				if "~forward" in line:
					collect_forward 	= True
				elif collect_forward and "~reverse" not in line:
					if not hasattr(I, "forward"):
						setattr(I, "forward", list())
					x,y 	= line.strip("\n").split(",")
					getattr(I, "forward").append((float(x),float(y)))

				elif collect_forward:
					collect_forward =False
				else:
					if not hasattr(I, "reverse"):
						setattr(I, "reverse", list())
					x,y 	= line.strip("\n").split(",")
					getattr(I, "reverse").append((float(x),float(y)))
	for I in G.values():
		model 		= BIC.get_best_model(I, 100, 2)
		print model.diff, len(model.rvs)
		model.rvs 	= [rv for rv in model.rvs if rv.w >0.01]
		print len(model.rvs)
		minX 		= min(min(I.forward)[0], min(I.reverse)[0])
		I.forward 	= [((x-minX)/100.,y) for x,y in I.forward]
		I.reverse 	= [((x-minX)/100.,y) for x,y in I.reverse]
		
		f_ct,f_e 	= np.histogram([x for x,y in I.forward ], weights=[y for x,y in I.forward ], bins=200, normed=1)
		r_ct,r_e 	= np.histogram([x for x,y in I.reverse ], weights=[y for x,y in I.reverse ], bins=200, normed=1)
		pi 			= sum(f_ct) / (sum(f_ct) + sum(r_ct))
		minX,maxX 	= min(min(f_e), min(r_e)),max(max(f_e), max(r_e) )
		xs 			= np.linspace(minX, maxX, 1000)
		f_ys 		= map(lambda x: model.pdf(x,1), xs)
		r_ys 		= map(lambda x: -model.pdf(x,-1), xs)
		plt.plot(xs, f_ys,linewidth=2, color="black")
		plt.plot(xs, r_ys,linewidth=2, color="black")
		plt.bar((f_e[:-1] + f_e[1:]) / 2., f_ct, width=(maxX-minX)/200.,alpha=0.3)
		plt.bar((r_e[:-1] + r_e[1:]) / 2., -r_ct, width=(maxX-minX)/200., color="red", alpha=0.3)
		plt.show()



if __name__ == "__main__":
	EMG_IN 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/EMG_out_files/test_file_2.tsv"
	EMG_OUT = "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/EMG_out_files/model_fits_out_all_11"
	run(EMG_IN, EMG_OUT)
	pass





