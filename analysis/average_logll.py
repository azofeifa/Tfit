import matplotlib.pyplot as plt
import load,math
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.ticker
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
import random as r
rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from scipy import stats

def show(L):
	LLS 	= [(-math.log((S.ll / S.N)*-1), math.log(S.N   )) for chrom in L for st,sp,S in L[chrom] if S.N and S.ll < 0]
	F 		= plt.figure(figsize=(15,10))
	ax1 	= F.add_subplot(1,2,1)
	ax2 	= F.add_subplot(1,2,2)
	x,y 	= [d for ll, d in LLS], [ll for ll, d in LLS]
	xy = np.vstack([x,y])
	
	z = gaussian_kde(xy)(xy)
	bins 			= 50
	counts,edges 	= np.histogram([ll for ll, d in LLS], bins=bins)
	counts 			= np.array([float(ct) for ct in counts])
	edges 			= (edges[:-1] + edges[1:])/2.
	CDF 			= [sum(counts[:i+1]) / sum(counts) for i in range(len(edges))]
	counts/=float(max(counts))
	ax1.bar(edges,counts, alpha=0.5, width=(edges[-1]-edges[0])/bins )
	bins 			= 300
	counts,edges 	= np.histogram([ll for ll, d in LLS], bins=bins)
	counts 			= np.array([float(ct) for ct in counts])
	edges 			= (edges[:-1] + edges[1:])/2.
	CDF 			= [sum(counts[:i+1]) / sum(counts) for i in range(len(edges))]
	ax1.set_title("Log-Likelihood / N")
	ax1.plot(edges, CDF, linewidth=2.5, color="black")
	ax1.grid()
	ro 				= np.corrcoef(x,y=y)
	ax2.set_title("Log-Likelihood / N vs N")
	ax2.scatter(x,y, c=z,s=14, edgecolor='', label=str(ro[0,1]) )
	ax2.legend()
	ax2.grid()
	plt.show()
def bootstraps(L):
	F 		= plt.figure()
	ax1 	= F.add_subplot(3,2,1)
	ax1.set_title("mu")
	ax1.hist([S.variances[0]*100 for chrom in L for st, sp, S in L[chrom] if S.variances[0] < 5  ], bins=50)
	ax2 	= F.add_subplot(3,2,2)
	ax2.set_title("sigma")
	ax2.hist([S.variances[1]*100 for chrom in L for st, sp, S in L[chrom] if S.variances[1] < 5  ], bins=50)
	ax3 	= F.add_subplot(3,2,3)
	ax3.set_title("lambda")
	ax3.hist([S.variances[2]*100 for chrom in L for st, sp, S in L[chrom] if S.variances[2] < 5 ], bins=50)
	ax4 	= F.add_subplot(3,2,4)
	ax4.set_title("weight")
	ax4.hist([S.variances[3] for chrom in L for st, sp, S in L[chrom]  ], bins=50)
	ax5 	= F.add_subplot(3,2,5)
	ax5.set_title("pi")
	ax5.hist([S.variances[4] for chrom in L for st, sp, S in L[chrom]  ], bins=50)
	ax6 	= F.add_subplot(3,2,6)
	ax6.set_title("foot print")
	ax6.hist([S.variances[5] for chrom in L for st, sp, S in L[chrom]  ], bins=50)
	plt.tight_layout()
	plt.show()
def correlate_density(L):
	F 		= plt.figure()
	ax1 	= F.add_subplot(2,2,1)
	ax1.set_title(r"$\mu$" )
	xy 		= [(math.log(S.variances[0]*10, 10), math.log(S.N+1,10)) for chrom in L for st, sp, S in L[chrom] if r.uniform(0,1) < 0.1  ]
	x,y 	= [y for x,y in xy],[x for x,y in xy]
	xy 		= np.vstack([x,y])
	ro 		= np.corrcoef(x,y=y)
	z 		= gaussian_kde(xy)(xy)
	ax1.scatter(x,y, s=10, c=z, edgecolor='', label=str(ro[0,1])[:6] )
	ax1.legend(loc=(0.1,0.1))
	ax1.set_ylabel(r'$\log_{10}(bp)$')
	ax2 	= F.add_subplot(2,2,2)
	ax2.set_title(r"$\sigma^2$"   )
	xy 		= [(math.log(S.variances[1]*10, 10), math.log(S.N+1,10)) for chrom in L for st, sp, S in L[chrom] if r.uniform(0,1) < 0.1  ]
	x,y 	= [y for x,y in xy],[x for x,y in xy]
	xy 		= np.vstack([x,y])
	ro 		= np.corrcoef(x,y=y)
	z 		= gaussian_kde(xy)(xy)
	ax2.scatter(x,y,c=z,s=10, edgecolor='', label=str(ro[0,1])[:6] )
	ax2.legend(loc=(0.1,0.1))
	ax2.set_ylabel(r'$\log_{10}(bp)$')
	
	ax3 	= F.add_subplot(2,2,3)
	ax3.set_title(r"$\lambda$"  )
	xy 		= [(math.log(S.variances[2]*10, 10), math.log(S.N+1,10)) for chrom in L for st, sp, S in L[chrom] if r.uniform(0,1) < 0.1 and S.variances[2]>0  ]
	x,y 	= [y for x,y in xy],[x for x,y in xy]
	xy 		= np.vstack([x,y])
	ro 		= np.corrcoef(x,y=y)
	z 		= gaussian_kde(xy)(xy)
	ax3.set_ylabel(r'$\log_{10}(bp)$')

	ax3.scatter(x,y,c=z,s=10, edgecolor='', label=str(ro[0,1])[:6], )
	ax3.legend(loc=(0.1,0.1))
	ax4 	= F.add_subplot(2,2,4)
	ax4.set_title(r"$\pi$")
	xy 		= [(S.variances[4] , math.log(S.N+1, 10)) for chrom in L for st, sp, S in L[chrom] if r.uniform(0,1) < 0.1   ]
	x,y 	= [y for x,y in xy],[x for x,y in xy]
	xy 		= np.vstack([x,y])
	ro 		= np.corrcoef(x,y=y)
	z 		= gaussian_kde(xy)(xy)
	ax4.scatter(x,y,c=z,s=10, edgecolor='', label=str(ro[0,1])[:6], )
	ax4.legend(loc=(0.1,0.1))
	for ax in (ax1, ax2,ax3,ax4):
		ax.grid()

	F.text(0.5, 0.04, "Data Size\n\n"+ r'$\log_{10}$'+'(coverage)', ha='center', va='center')
	F.text(0.06, 0.5, 'Bootstrap\nSample Variance', ha='center', va='center', rotation='vertical')
	plt.show()


if __name__ == "__main__":
	DIR 			="/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/"
		
	REF 						= "/Users/joazofeifa/Lab/genome_files/RefSeqHG19.txt"
	DMSO2_3 					="Allen2014_boot_DMSO2_3-2_bootstrapped_bidirectional_hits_intervals.bed"
	DMSO2_3_L,DMSO2_3_G 			= load.load_model_fits_bed_file(DIR+DMSO2_3)
	#bootstraps(DMSO2_3_L)
	correlate_density(DMSO2_3_L)
	