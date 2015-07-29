import os
import merge_data_types as mdt
import BIC
import matplotlib.pyplot as plt
#=============================================
#GLOBAL VARIABLES
si_thresh 	= 4
l_thresh 	= 2
w_thresh 	= 0.01
pi_thresh 	= 0.1
#=============================================
penality 		= 50
diff_threshold  = 5


def check_bidir_component(rv, si_thresh=4, l_thresh=2., w_thresh=0.01, pi_thresh=0.1):
	if rv.type == "N":
		if rv.si < si_thresh and rv.l < l_thresh and rv.w > w_thresh and pi_thresh <= rv.pi<= (1-pi_thresh):
			return True
	return False


def output(I, G):
	model 	= BIC.get_best_model(I, penality , diff_threshold)
	bidirs 	= [rv for rv in model.rvs if check_bidir_component(rv,si_thresh=si_thresh, l_thresh=l_thresh, w_thresh=w_thresh, pi_thresh=pi_thresh)]
	for c in bidirs:
		G["mu"].append(c.mu)
		G["si"].append(c.si)
		G["l"].append(c.l)
		G["pi"].append(c.pi)
		
	pass
def parse_file(FILE, G):
	with open(FILE) as FH:

		header,I	= True,None
		for line in FH:
			if header:
				if "#" == line[0]:
					if I is not None:
						output(I, G)
					chrom,info 			= line[1:].strip("\n").split(":")
					start_stop, N  		= info.split(",")
					start,stop 			= start_stop.split("-")
					I 					= mdt.segment(chrom,int(start),int(stop),float(N) )
				elif "~" == line[0]:
					I.insert_model_info(line)
				elif "N:"==line[:2] or "U:"==line[:2]:
					I.insert_component(line)
			else:
				header=False


def load(DIR,test=3):
	G 	= {"mu":list(), "si":list(), "l":list(), "pi":list()}
	t 	= 0
	for FILE in os.listdir(DIR):
		parse_file(DIR+FILE, G)
		if test is not None and t > test:
			break
		t+=1
	return G

def run(DIR):
	G 	= load(DIR)
	F 	= plt.figure(figsize=(15,10))
	ax1 = F.add_subplot(2,2,1)
	ax1.set_title("Displacement of mu")

	ax2 = F.add_subplot(2,2,2)
	ax2.set_title("sigma")
	ax2.hist(G["si"], bins=100)
	
	ax3 = F.add_subplot(2,2,3)
	ax3.set_title("lambda")
	ax3.hist(G["l"], bins=100)
	ax4 = F.add_subplot(2,2,4)
	ax4.set_title("pi")
	ax4.hist(G["pi"], bins=100)
	plt.show()



if __name__ == "__main__":
	DIR 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/EMG_out_files/DMSO_ND_intervals_model_fits_1/"
	run(DIR)
	pass