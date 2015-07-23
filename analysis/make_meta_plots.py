import numpy as np
import matplotlib.pyplot as plt
def load(FILE, w_thresh=0.01, s_thresh=2, l_thresh=2., pi_thresh=0.05, PEAK=False):
	G 	= {}
	N 	= 0
	with open(FILE) as FH:
		collect=False
		for line in FH:
			if "#" == line[0]:
				collect = False
			elif "N: " == line[:3]:
				mu, si, l, w, pi 	= line[3:].strip("\n").split(",")
				mu, si, l, w, pi 	= float(mu), float(si), float(l), float(w), float(pi)
				if si < s_thresh and l < l_thresh and pi_thresh<= pi<=(1-pi_thresh):
					collect=True
			elif collect:
				lineArray = line.strip("\n").split(",")
				data_type, peak, data 	= lineArray[0], lineArray[1], lineArray[2:]
				if data_type not in G:
					G[data_type] 	= np.zeros((len(data), ))
				NN 	= sum(np.array([float(d.split("-")[1]) for d in data ]))
				if NN:
					G[data_type]+=((np.array([float(d.split("-")[1]) for d in data ]))/NN)
				N+=1
	for data_type in G:
		F 	= plt.figure(figsize=(15,10))
		plt.title(data_type)
		plt.bar(np.linspace(-4,4, len(G[data_type])), G[data_type] , width=8./len(G[data_type]) )
		plt.grid()
		plt.show()







if __name__ == "__main__":
	FILE 	= "/Users/joeyazo/Desktop/Lab/EMG_files/old_counted_data_100_5_4_2_0.01_0.1"
	load(FILE)