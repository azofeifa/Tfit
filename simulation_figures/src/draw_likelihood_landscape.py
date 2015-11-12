import model, simulate, matplotlib.pyplot as plt
import numpy as np, math
import matplotlib as mpl
def log_likelihood(X, rv):
	ll 	= 0
	for i in range(X.shape[0]):
		vl1, vl2 	=  rv.pdf(X[i,0],1)+0.001, rv.pdf(X[i,0],-1)+0.001
		if X[i,1] and vl1:
			ll+=math.log(vl1)*X[i,1]
		if X[i,2] and vl2:
			ll+=math.log(vl2)*X[i,2]

	if not ll:
		return -100
	return ll
			


def mu_by_lambda(ax,X,res):
	A 		= np.zeros((res, res))
	for i,mu in enumerate(np.linspace(-4, 4,res)):
		for j, l in enumerate(np.linspace(1,10,res)):
			l 	= 1.0 / l
			rv 	= model.component_bidir(mu,7.3, l, 1.0, 0.5, None, foot_print=4)
			ll 		= abs(log_likelihood(X, rv))
			A[i,j] 	= math.log(ll,10)
	A 		= A[:,::-1]
	heatmap = ax.imshow(A, cmap=plt.cm.jet_r,vmin=A.min(), vmax=A.max(),aspect=0.85 )
	ax.set_xticklabels([str(x)[:4] for x in  np.linspace(1,10,len(ax.get_xticklabels())) ], rotation=45)
	ax.set_yticklabels([str(x)[:4] for x in  np.linspace(-4, 4,len(ax.get_xticklabels())) ], rotation=45)
	ax.set_xlabel(r'$1/\lambda$',fontsize=20)
	ax.set_ylabel(r'$\mu$',fontsize=20)
def sigma_by_lambda(ax,X,res):
	A 		= np.zeros((res, res))
	for i,si in enumerate(np.linspace(1, 7,res)):
		for j, l in enumerate(np.linspace(1,20,res)):
			l 	= 1.0 / l
			rv 		= model.component_bidir(0,si, l, 1.0, 0.5, None, foot_print=4)
			ll 		= abs(log_likelihood(X, rv))
			A[i,j] 	= math.log(ll,10)
	heatmap = ax.imshow(A, cmap=plt.cm.jet_r,vmin=A.min(), vmax=A.max(),aspect=0.85 )
	ax.set_xticklabels([str(x)[:4] for x in  np.linspace(1,10,len(ax.get_xticklabels())) ], rotation=45)
	ax.set_yticklabels([str(x)[:4] for x in  np.linspace(1, 7,len(ax.get_xticklabels())) ], rotation=45)
	ax.set_xlabel(r'$1/\lambda$',fontsize=20)
	ax.set_ylabel(r'$\sigma$',fontsize=20)

def fp_by_lambda(ax,X,res):
	A 		= np.zeros((res, res))
	for i,fp in enumerate(np.linspace(-2, 10,res)):
		for j, l in enumerate(np.linspace(1,20,res)):
			l 	= 1.0 / l
			rv 	= model.component_bidir(0,5, l, 1.0, 0.5, None, foot_print=fp)
			ll 		= abs(log_likelihood(X, rv))
			A[i,j] 	= math.log(ll,10)
	heatmap = ax.imshow(A, cmap=plt.cm.jet_r,vmin=A.min(), vmax=A.max(),aspect=0.85 )
	ax.set_xticklabels([str(x)[:4] for x in  np.linspace(1,50,len(ax.get_xticklabels())) ], rotation=45)
	ax.set_yticklabels([str(x)[:4] for x in  np.linspace(-2, 10,len(ax.get_xticklabels())) ], rotation=45)
	ax.set_xlabel(r'$\frac{1}{\lambda}$',fontsize=20)
	ax.set_ylabel(r'$fp$',fontsize=20)

def fp_by_sigma(ax,X,res):
	A 		= np.zeros((res, res))
	for i,fp in enumerate(np.linspace(0, 8,res)):
		for j, si in enumerate(np.linspace(1,10,res)):
			rv 	= model.component_bidir(0,si, 0.1, 1.0, 0.5, None, foot_print=fp)
			ll 		= abs(log_likelihood(X, rv))
			A[i,j] 	= math.log(ll,10)
	heatmap = ax.imshow(A, cmap=plt.cm.jet_r,vmin=A.min(), vmax=A.max(),aspect=0.85 )

	ax.set_xticklabels([str(x)[:4] for x in  np.linspace(1,7,len(ax.get_xticklabels())) ], rotation=45)
	ax.set_yticklabels([str(x)[:4] for x in  np.linspace(0, 8,len(ax.get_xticklabels())) ], rotation=45)
	ax.set_xlabel(r'$\sigma$',fontsize=20)
	ax.set_ylabel(r'$fp$',fontsize=20)
def mu_by_pi(ax,X,res):
	A 		= np.zeros((res, res))
	for i,mu in enumerate(np.linspace(-8, 8,res)):
		for j, pi in enumerate(np.linspace(0,1,res)):
			rv 	= model.component_bidir(mu,5, 1.0/5, 1.0, pi, None, foot_print=4)
			ll 		= -(log_likelihood(X, rv))

			A[i,j] 	= math.log(ll,10)
	heatmap = ax.imshow(A, cmap=plt.cm.jet_r,vmin=A.min(), vmax=A.max(),aspect=0.85 )
	ax.set_xticklabels([str(x)[:4] for x in  np.linspace(0,1,len(ax.get_xticklabels())) ], rotation=45)
	ax.set_yticklabels([str(x)[:4] for x in  np.linspace(-8, 8,len(ax.get_xticklabels())) ], rotation=45)
	ax.set_ylabel(r'$\mu$',fontsize=20)
	ax.set_xlabel(r'$\pi$',fontsize=20)

def draw():
	X 	= simulate.runOne(SHOW=False,we=1.0,wl=0., wr=0., foot_print=4,s=5, l=5)
	F 	= plt.figure( )
	ax1 = F.add_subplot(2,2,1)
	ax2 = F.add_subplot(2,2,2)
	ax3 = F.add_subplot(2,2,3)
	ax4 = F.add_subplot(2,2,4)
	ax_c= F.add_axes([0.35,0.52,0.35,0.02])
	ax_c.set_xticks([])
	ax_c.set_yticks([])
	cmap = mpl.cm.jet
	norm = mpl.colors.Normalize(vmin=-1000, vmax=0)

	# ColorbarBase derives from ScalarMappable and puts a colorbar
	# in a specified axes, so it has everything needed for a
	# standalone colorbar.  There are many more kwargs, but the
	# following gives a basic continuous colorbar with ticks
	# and labels.
	cb1 = mpl.colorbar.ColorbarBase(ax_c, cmap=cmap,
	                            norm=norm,
	                            orientation='horizontal')
	ax_c.set_xticklabels([r"$-10^{30}$"]+ [" "," "," "," ",r"$-10^{15}$"," "," "," "," " ]+ [r'$-10^3$'] , rotation=15 )
	
	ax_c.set_title("log-likelihood")
	res = 3
	mu_by_lambda(ax1, X,res)
	mu_by_pi(ax2,X,res)
	sigma_by_lambda(ax3,X,res)
	fp_by_sigma(ax4,X,res)


	plt.tight_layout()
	plt.show()

if __name__ == "__main__":
	draw()
	pass


