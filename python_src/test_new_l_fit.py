import simulate
from model import component_bidir, component_elongation
import numpy as np
import matplotlib.pyplot as plt
import load
import math as m
def grab_data(forward, reverse):
	D 	= (list(),list())
	for i,FILE in enumerate((forward, reverse)):
		with open(FILE) as FH:
			for line in FH:
				chrom,start, stop, cov 	= line.strip("\n").split("\t")
				x 						= (int(stop) + int(start ))/2.
				D[i].append((x, float(cov)))
	X 	= simulate.BIN(D[0], D[1], 200)
	X[:,0] 	-= X[0,0]
	X[:,0] /=100.
	return X
def draw(components, X):
	F 			= plt.figure(figsize=(15,10))
	xs 	= np.linspace(X[0,0], X[-1,0], 1000)
	ysf = map(lambda x: sum([c.pdf(x,1) for c in components]) , xs )
	ysr = map(lambda x: -sum([c.pdf(x,-1) for c in components]) , xs )
	bins 		= 300
	fct, fe 	= np.histogram(X[:,0], weights=X[:,1,] ,normed=1, bins=bins)
	fe 			= (fe[:-1] + fe[1:])/2.
	rct, re 	= np.histogram(X[:,0], weights=X[:,2], normed=1, bins=bins)
	re 			= (re[:-1] + re[1:])/2.
	widthf 		= (fe[-1] - fe[0])/bins
	widthr 		= (re[-1] - re[0])/bins
	
	plt.plot(xs, ysf,linewidth=2., color="black")
	plt.plot(xs, ysr,linewidth=2., color="black")
	plt.bar(fe, fct, width=widthf, edgecolor="blue", color="blue", alpha=0.5)
	plt.bar(re, -rct, width=widthr, edgecolor="red", color="red", alpha=0.5)
	plt.grid()
	plt.show()
def compute_log_likelihood(X, components,  move_as=None, move_bs=None):
	forward, reverse 	= 0. , 0.
	if move_as is None:
		move_as 		= np.zeros((len(components, )))
	if move_bs is None:
		move_bs 		= np.zeros((len(components, )))
	
	for i in range(X.shape[0]):
		x,forward_ct, reverse_ct 	= X[i,:]
		forward+=m.log(sum([c.pdf(x,1, move_a=move_as[i], move_b=move_bs[i]) for i,c in enumerate(components)]))*forward_ct

	for i in range(X.shape[0]):
		x,forward_ct, reverse_ct 	= X[i,:]
		reverse+=m.log(sum([c.pdf(x,-1, move_a=move_as[i], move_b=move_bs[i]) for i,c in enumerate(components)]))*reverse_ct
	return forward + reverse




def move_l(X, parameters):
	minX, maxX 	= X[0,0], X[-1,0]
	K 			= len(parameters)
	noise_max 	= 0.01
	bidirs 		= [component_bidir(mu, si,l, w,pi,None) for mu, si,l, w, pi in parameters] 
	right_bounds= [ mu for mu, si,l, w, pi in parameters[1:] ] + [maxX]
	left_bounds = [minX] + [ mu for mu, si,l, w, pi in parameters[:-1] ] 
	uniforms    = [component_elongation(np.random.uniform(left_bounds[k], parameters[k][0]), parameters[k][0], 0.111, 0., bidirs[k], "reverse",None ,0 ) for k in range(K)]
	uniforms   += [component_elongation(parameters[k][0],np.random.uniform(parameters[k][0], right_bounds[k]) , 0.111, 1., bidirs[k], "forward",None, X.shape[0] ) for k in range(K)]
	uniforms   += [component_elongation(minX, maxX, noise_max, 0.5, bidirs[0], "noise", None, X.shape[0]) ]
	components 	= bidirs + uniforms

	ll 			= compute_log_likelihood(X, components)
	t 			= 0
	EX 			= np.zeros((len(components), ))
	weights 	= np.zeros((len(components),))
	current 	= np.zeros((len(components), 2))
	prevll 		= -np.inf
	converged 	= False
	while t < 300 and not converged:
		#update weights
		for i in range(X.shape[0]):
			x,f,r 		= X[i,:]
			normed_f 	= 0
			normed_r 	= 0
			for k,c in enumerate(components):
				prob_f 	= c.pdf(x,1)
				prob_r 	= c.pdf(x,-1)
				normed_f+=prob_f
				normed_r+=prob_r
				current[k,0] 	= prob_f
				current[k,1] 	= prob_r
			for k,c in enumerate(components):
				if (normed_f):
					EX[k]+=(current[k,0]/normed_f)*f
				if (normed_r):
					EX[k]+=(current[k,1]/normed_r)*r
		N 	= 0
		for k,c in enumerate(components):
			N+=EX[k]
		W 	= 0
		for k,c in enumerate(components):

			c.w 	= EX[k]/N	
			if c.type == "noise":
				c.w = min(c.w, noise_max)
			else:
				W+=c.w
			EX[k] 	= 0
		for c in components:
			if c.type!="noise":
				c.w/=W
		#now move lls
		ll 			= compute_log_likelihood(X, components)
		change 		= False
		for c in components:

			if c.type=="forward":
				oldb 	= c.b
				c.b 	+=np.random.normal(0,10)
				nll 	= compute_log_likelihood(X, components)
				if nll < ll:
					c.b = oldb
				else:
					change=True
					ll 	= nll
			elif c.type=="reverse":
				olda 	= c.a
				c.a 	+=np.random.normal(0,10)
				nll 	= compute_log_likelihood(X, components)
				if nll < ll:
					c.a = olda
				else:
					change=True
					ll 	= nll
		print t


		if abs(ll - prevll) < pow(10,-1):
			converged 	= True
			draw(components, X)

		prevll 		= ll



		t+=1

	draw(components, X)




if __name__=="__main__":
	# D 	= simulate.runOne(mu=0, s=1, l=5, lr=100, ll=-100, we=0.5,wl=0.25, 
	# 	wr=0.25, pie=0.5, pil=0.1, pir=0.9, N=1000, SHOW=False , bins=200, noise=True )
	forward 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/chr1_6229860_6303055_DMSO2_3.pos.BedGraph"
	reverse 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/chr1_6229860_6303055_DMSO2_3.neg.BedGraph"
	X 			= grab_data(forward, reverse)
	parameters 	= [	(298.119448,1.036515,0.228770,0.111111,0.423449), 
					(358.097039,3.538068,0.129734,0.111111,0.391122), 
					(660.517058,0.063181,0.187830,0.111111,0.379944)]
	move_l(X, parameters)