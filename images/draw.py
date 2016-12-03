import sys,numpy as np
import matplotlib.pyplot as plt
import math,time
import scipy.stats as ss
def load(FILE):
	X,Y = list(),list()
	for line in open(FILE):
		line_array = line.strip("\n").split("\t")
		vals       = map(float, line_array[3:])
		if vals[-1] and vals[-2] and np.random.uniform(0,1) < 0.5:
			x,y 		=  math.log(sum(vals[1:]),10),vals[0]
			X.append(x)
			Y.append(y+0.65)

	XY 	= zip(X,Y)
	return XY
def normal(x, mu,sigma):
	return 1.0/(math.sqrt(2.0*math.pi)*sigma)*math.exp(-pow(x-mu,2)/(2.0*pow(sigma,2)))
def exponential(x, lam, mu):
	if x >= mu:
		return lam*math.exp(-(x-mu)*lam)
	return 0.0
def pareto(x,mu, alpha):
	if (x >= mu):
		return alpha*pow(mu,alpha) / pow(x, alpha+1)
	return 0.0


def EM(X,P=True, E=False):
	mu,sigma,lam,w,alpha 			= 0.5,0.1,1.0,0.1,5.2
	t,T 							= 0,1030
	l 								= X[-1,0]
	while t < T:
		EX,EX2,EL, EY,R1, R2 	= 0.0, 0.0,0.0,0.0,0.0,0.0
		for i in range(X.shape[0]):
			x,y 				= X[i,:]
			if E:
				p1, p2 			= (1-w)*normal(x,mu,sigma),w*exponential(x, lam,mu)
			elif P:
				p1, p2 			= (1-w)*normal(x,mu,sigma),w*pareto(x, mu,alpha)
			if p1 + p2:
				r1,r2 			= p1 / (p1 + p2) , p2 / (p1 + p2)
				R1,R2 			= R1 + r1*y,R2 + r2*y
				EX 				+=(r1*x)*y
				if E:
					EY 				+=(r2*x)*y
				elif P:
					EL 				+= (math.log(x) - math.log(mu))*p2*y
				EX2 				+=(r1*pow(x-mu,2))*y
		mu 	= EX /R1
		sigma = math.sqrt(EX2 / R1)
		if E:
			lam 	= R2/EY
		elif P:
			alpha = R2 / EL
		w 		= R2 / (R1 + R2)
		t+=1
	if E:
		f 	= lambda x: (1-w)*normal(x,mu, sigma) + w*exponential(x,lam, mu)
		h,g= lambda x: (1-w)*normal(x,mu, sigma), lambda x: w*exponential(x,lam, mu)
	else:
		f 	= lambda x: (1-w)*normal(x,mu, sigma) + w*pareto(x,mu,alpha)
		h,g= lambda x: (1-w)*normal(x,mu, sigma), lambda x: w*pareto(x,mu,alpha)
	return mu,sigma, w,f,h,g

FILE = "/Users/joazofeifa/test.scores"
XY 	 	= load(FILE)
counts,edges 	= np.histogram([y for x,y in XY],bins=300,normed=1)
X 					= np.zeros((len(counts), 2))
X[:,0] 			= (edges[1:] + edges[:-1])/2.0
X[:,1] 			= counts

mu,sigma, w,f,h,g 	 = EM(X)
F 		= plt.figure(facecolor="white",figsize=(15,6))
ax1 	= plt.gca()
ax1.bar(X[:,0], X[:,1], width=(X[-1,0]-X[0,0])/float(len(counts)),
				alpha=0.25,edgecolor="white",label="random genomic\ndraws")
xs 					 = np.linspace(0.5, X[-1,0], 400)
ys 					 = map(f, xs)
ax1.set_title("$\{\mu,\sigma,w\}=$" + str(mu)[:4] + "," + str(sigma)[:4] + "," + str(w)[:4] )
#ax1.plot(xs, ys,lw=2,label="mixture")
ax1.plot(xs, map(h,xs) ,lw=2,label="normal")
ax1.plot(xs, map(g,xs) ,lw=2,label="pareto")
ax1.set_xlabel("$LLR$",fontsize=20)
ax1.set_ylabel("Normalized Frequency\nDensity",fontsize=20)
# Hide the right and top spines
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax1.legend(loc="best")
plt.show()
