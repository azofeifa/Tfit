import matplotlib.pyplot as plt
import time
import numpy as np
import math
from scipy.stats import gaussian_kde
import matplotlib.ticker
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

def match_UP(A_L, B_L):
	A,B ={}, {}

	for chrom in A_L:
		if chrom not in A:
			A[chrom] 	= list()
		for st, sp, S in A_L[chrom]:
			for m in S.models:
				st, sp 	= m.mu - m.std, m.mu+m.std
				m.annotated 	= S.annotated
				m.p53_site 		= S.p53_site
				m.dist 			= S.dist
			
				A[chrom].append((st,sp, m))
	for chrom in B_L:
		if chrom not in B:
			B[chrom] 	= list()
		for st, sp, S in B_L[chrom]:
			for m in S.models:
				m.annotated 	= S.annotated
				m.p53_site 		= S.p53_site
				m.dist 			= S.dist
				st, sp 	= m.mu - m.std, m.mu+m.std
				B[chrom].append((st,sp, m))
	overlaps 	= {}
	for chrom in A:
		overlaps[chrom] 	= list()
		if chrom in B:
			a,b 	= A[chrom],B[chrom]
			j,N 	= 0,len(b)
			for i,(a_st, a_sp, a_m) in enumerate(a):
				while j < N and b[j][1] < a_st:
					j+=1
				u 	= j
				o 	= list()
				while u < N and b[u][0] < a_sp:
					o_st,o_sp 	= max(b[u][0], a_st), min(b[u][1], a_sp)
					o.append(( o_sp-o_st, a[i],  b[u] ))
					u+=1
				if o:
					o 	= max(o)
					if abs(o[1][2].mu - o[2][2].mu) < 150:
						overlaps[chrom].append((o[1], o[2] ))
	return [(a[2], b[2]) for chrom in overlaps for a,b in overlaps[chrom] ]
				

def overlaps_and_not(A_L, B_L):
	A,B ={}, {}

	for chrom in A_L:
		if chrom not in A:
			A[chrom] 	= list()
		for st, sp, S in A_L[chrom]:
			for m in S.models:
				st, sp 	= m.mu - m.std, m.mu+m.std
				m.annotated 	= S.annotated
				m.p53_site 		= S.p53_site
				m.dist 			= S.dist
				A[chrom].append((st,sp, m))
	for chrom in B_L:
		if chrom not in B:
			B[chrom] 	= list()
		for st, sp, S in B_L[chrom]:
			for m in S.models:
				m.annotated 	= S.annotated
				m.p53_site 		= S.p53_site
				m.dist 			= S.dist
				st, sp 	= m.mu - m.std, m.mu+m.std
				B[chrom].append((st,sp, m))
	overlaps 	= {}
	for chrom in A:
		overlaps[chrom] 	= list()
		if chrom in B:
			a,b 	= A[chrom],B[chrom]
			j,N 	= 0,len(b)
			for i,(a_st, a_sp, a_m) in enumerate(a):
				while j < N and b[j][1] < a_st:
					j+=1
				u 	= j
				o 	= list()
				while u < N and b[u][0] < a_sp:
					
					a_m.overlap  	= True
					b[u][2].overlap = True

					u+=1
				if o:
					o 	= max(o)
	return A,B

def run_all(overlaps):

	F 	= plt.figure(figsize=(15,10))
	ax1 = F.add_subplot(2,2,1)
	ax1.set_title("Loading Variance")
	attr 			= "si"
	overlaps 		= [(x,y) for x,y in overlaps if getattr(x, attr) is not None and  getattr(y, attr) is not None ]
	overlaps 		= [(x,y) for x,y in overlaps if getattr(x, "wEM") >0.1 and  getattr(y, "wEM") > 0.1]
	x,y 			= [math.log(getattr(x, attr),10) for x,y in overlaps],[math.log(getattr(y, attr),10) for x,y in overlaps]
	xy = np.vstack([x,y])
	
	z = gaussian_kde(xy)(xy)
	
	ax1.scatter(x, y, c=z, s=14, edgecolor='' )
	ax1.grid()
	ax1.set_xlabel("HCT116 GRO-seq Rep1, log10 bp")
	ax1.set_ylabel("HCT116 GRO-seq Rep2, log10 bp")

	ax2 = F.add_subplot(2,2,2)
	ax2.set_title("Initiating Length")

	attr 			= "lam"
	overlaps 		= [(x,y) for x,y in overlaps if getattr(x, attr) is not None and  getattr(y, attr) is not None ]
	overlaps 		= [(x,y) for x,y in overlaps if getattr(x, "wEM") >0.1 and  getattr(y, "wEM") > 0.1]
	x,y 			= [math.log(getattr(x, attr),10) for x,y in overlaps],[math.log(getattr(y, attr),10) for x,y in overlaps]
	xy = np.vstack([x,y])
	
	z = gaussian_kde(xy)(xy)
	
	ax2.scatter(x, y, c=z, s=14, edgecolor='' )
	ax2.grid()
	ax2.set_xlabel("HCT116 GRO-seq Rep1, log10 bp")
	ax2.set_ylabel("HCT116 GRO-seq Rep2, log10 bp")

	ax3 = F.add_subplot(2,2,3)
	ax3.set_title("Strand Bias")

	attr 			= "pi"
	overlaps 		= [(x,y) for x,y in overlaps if getattr(x, attr) is not None and  getattr(y, attr) is not None ]
	overlaps 		= [(x,y) for x,y in overlaps if getattr(x, "wEM") >0.1 and  getattr(y, "wEM") > 0.1]
	x,y 			= [getattr(x, attr)  for x,y in overlaps],[ getattr(y, attr)  for x,y in overlaps]
	xy = np.vstack([x,y])
	
	z = gaussian_kde(xy)(xy)
	
	ax3.scatter(x, y, c=z, s=14, edgecolor='' )
	ax3.grid()
	ax3.set_xlabel("HCT116 GRO-seq Rep1")
	ax3.set_ylabel("HCT116 GRO-seq Rep2")
	
	ax4 = F.add_subplot(2,2,4)
	ax4.set_title("Paused Probability")

	attr 			= "wEM"
	overlaps 		= [(x,y) for x,y in overlaps if getattr(x, attr) is not None and  getattr(y, attr) is not None ]
	overlaps 		= [(x,y) for x,y in overlaps if getattr(x, "wEM") >0.1 and  getattr(y, "wEM") > 0.1]
	x,y 			= [getattr(x, attr)  for x,y in overlaps],[ getattr(y, attr)  for x,y in overlaps]
	xy = np.vstack([x,y])
	
	z = gaussian_kde(xy)(xy)
	
	ax4.scatter(x, y, c=z, s=14, edgecolor='' )
	ax4.grid()
	ax4.set_xlabel("HCT116 GRO-seq Rep1")
	ax4.set_ylabel("HCT116 GRO-seq Rep2")
	plt.tight_layout()
	plt.savefig("/Users/joazofeifa/Lab/Talks/2015/CSHL/GeneFig.svg")

def run(overlaps, attr ="si", LOG=False):


	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left+width+0.02
	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]
	F 			= plt.figure(figsize=(15,10))
	
	axScatter = plt.axes(rect_scatter)
	axHistx = plt.axes(rect_histx)
	axHisty = plt.axes(rect_histy)

	overlaps 		= [(x,y) for x,y in overlaps if getattr(x, attr) is not None and  getattr(y, attr) is not None ]

	overlaps 		= [(x,y) for x,y in overlaps if getattr(x, "wEM") >0.5 and  getattr(y, "wEM") > 0.5]
	
	A,B 			= [(x, math.log(getattr(x, attr) ,10)) if LOG else (x,getattr(x, attr)) for x,y in overlaps],[(y,math.log(getattr(y, attr) ,10)) if LOG else (y,getattr(y, attr)) for x,y in overlaps]
	
	x,y 			= [math.log(getattr(x, attr) ,10) if LOG else getattr(x, attr) for x,y in overlaps],[math.log(getattr(y, attr) ,10) if LOG else getattr(y, attr) for x,y in overlaps]
	
	var 			= np.mean([abs(a-b) for a,b in zip(x,y) ])

	AB 				= [(A[i][0], B[i][0], A[i][1], B[i][1] )  for i in range(len(A)) if abs(A[i][1] - B[i][1]) > 0*var  ]


	xy = np.vstack([x,y])
	
	z = gaussian_kde(xy)(xy)
	
	axScatter.scatter(x, y, c=z, s=14, edgecolor='' )
	axScatter.grid()
	
	axHistx.hist(x, bins=200, label=str(np.mean(x)))
	axHisty.hist(y, bins=200, orientation='horizontal', label=str(np.mean(y)))
	axHistx.legend()
	axHisty.legend()
	axHistx.set_xlim( axScatter.get_xlim() )
	axHisty.set_ylim( axScatter.get_ylim() )
	axHistx.grid()
	axHisty.grid()
	
	axHistx.set_xticks([])
	axHisty.set_yticks([])
	print "----------------------------------------------------------"
	print "All DMSO to ALL Nutlin" ,stats.ks_2samp(x,y)
	print "All DMSO to p53 DMSO", stats.ks_2samp(x,[ ab[2]  for ab in AB if ab[0].p53_site ] )
	print "All DMSO to p53 Nutlin", stats.ks_2samp(x,[ ab[3]  for ab in AB if ab[1].p53_site ] )
	print "All Nutlin to p53 DMSO", stats.ks_2samp(y, [ ab[2]  for ab in AB if ab[0].p53_site ]  )
	print "All Nutlin to p53 Nutlin", stats.ks_2samp(y,[ ab[3]  for ab in AB if ab[1].p53_site ] )
	print "----------------------------------------------------------"
		
	
	
	for i,ab in enumerate(AB):     
		label 	= AB[i][0].chrom + ":" + str(int(AB[i][0].start)) + "-" + str(int(AB[i][0].stop))                                           # <--
		xy 		= (AB[i][2], AB[i][3])
		
		#axScatter.annotate(label, xy=xy)
		if AB[i][0].p53_site:
			axScatter.scatter([xy[0]], [xy[1]], color="black", s=35)

	plt.show()

	pass

def si_lam(overlaps):

	X_si_wEM 	= [(a.lam, a.fp) for a,b in overlaps] + [(b.lam, b.fp) for a,b in overlaps]
	X_si_lam 	= [(a.si, a.lam) for a,b in overlaps] + [(b.si, b.lam) for a,b in overlaps]
	X_lam_wEM 	= [(a.lam, a.wEM) for a,b in overlaps] + [(b.lam, b.wEM) for a,b in overlaps]
	X_pi_wEM 	= [(a.si, a.fp) for a,b in overlaps] + [(b.si, b.fp) for a,b in overlaps]
	

	F 			= plt.figure(figsize=(15,10))
	
	ax1 		= F.add_subplot(2,2,1)
	x,y 		= [math.log(x,10) for x,y in X_si_wEM],[y for x,y in X_si_wEM]
	x,y 		= [x for x,y in X_si_wEM],[y for x,y in X_si_wEM]
	xy 			= np.vstack([x,y])
	z 			= gaussian_kde(xy)(xy)
	ax1.scatter(x, y, c=z, s=14, edgecolor='')
	ax1.set_xlabel("Variance in Loading")
	ax1.set_ylabel("Probability of Paused")
	
	ax1.grid()

	
	ax2 		= F.add_subplot(2,2,2)
	x,y 		= [math.log(x,10) for x,y in X_si_lam],[math.log(y,10) for x,y in X_si_lam]
	x,y 		= [x for x,y in X_si_lam],[y for x,y in X_si_lam]
	xy 			= np.vstack([x,y])
	z 			= gaussian_kde(xy)(xy)
	ax2.set_xlabel("Variance in Loading")
	ax2.set_ylabel("Length of Initiation")
	ax2.scatter(x, y, c=z, s=14, edgecolor='')
	ax2.grid()

	ax3 		= F.add_subplot(2,2,3)
	x,y 		= [math.log(x,10) for x,y in X_lam_wEM],[y for x,y in X_lam_wEM]
	xy 			= np.vstack([x,y])
	z 			= gaussian_kde(xy)(xy)
	ax3.set_xlabel("Length of Initiation")
	ax3.set_ylabel("Probability of Paused")
	ax3.scatter(x, y, c=z, s=14, edgecolor='')
	ax3.grid()

	ax4 		= F.add_subplot(2,2,4)
	x,y 		= [x for x,y in X_pi_wEM],[y for x,y in X_pi_wEM]
	xy 			= np.vstack([x,y])
	z 			= gaussian_kde(xy)(xy)
	
	ax4.set_xlabel("Strand Probability")
	ax4.set_ylabel("Probability of Paused")
	ax4.scatter(x, y, c=z, s=14, edgecolor='')
	ax4.grid()

	plt.show()


def promoter_differences_test(LS):
	sigmasA 	= [math.log(M.models[0].si, 10) for L in LS for chrom in L for st,sp, M in L[chrom] if M.annotated]
	sigmasU 	= [math.log(M.models[0].si, 10) for L in LS for chrom in L for st,sp, M in L[chrom] if not M.annotated]
	lamsA 		= [math.log(M.models[0].lam, 10) for L in LS for chrom in L for st,sp, M in L[chrom] if M.annotated]
	lamsU 		= [math.log(M.models[0].lam, 10) for L in LS for chrom in L for st,sp, M in L[chrom] if not M.annotated]
	weightsA 	= [ M.models[0].wEM  for L in LS for chrom in L for st,sp, M in L[chrom] if M.annotated]
	weightsU 	= [ M.models[0].wEM  for L in LS for chrom in L for st,sp, M in L[chrom] if not M.annotated]
	pisA 		= [ M.models[0].pi  for L in LS for chrom in L for st,sp, M in L[chrom] if M.annotated]
	pisU 		= [ M.models[0].pi  for L in LS for chrom in L for st,sp, M in L[chrom] if not M.annotated]
	
	F 			= plt.figure(figsize=(15,10))
	ax1 		= F.add_subplot(2,2,1)
	ax1.set_title("Variance in Loading")
	ax1.hist(sigmasA, alpha=0.5, bins=100, normed=1, label="promoter")
	ax1.hist(sigmasU, alpha=0.5, bins=100, normed=1, label="unannotated")
	ax1.legend(loc=(0.1,0.8))
	ax1.grid()

	ax2 		= F.add_subplot(2,2,2)
	ax2.set_title("Initiating Length")
	ax2.hist(lamsA, alpha=0.5, bins=100, normed=1, label="promoter")
	ax2.hist(lamsU, alpha=0.5, bins=100, normed=1, label="unannotated")
	ax2.legend(loc=(0.1,0.8))
	ax2.grid()

	ax3 		= F.add_subplot(2,2,3)
	ax3.set_title("Paused Probability")
	ax3.hist(weightsA, alpha=0.5, bins=100, normed=1, label="promoter")
	ax3.hist(weightsU, alpha=0.5, bins=100, normed=1, label="unannotated")
	ax3.legend(loc=(0.1,0.8))
	ax3.grid()


	ax4 		= F.add_subplot(2,2,4)
	ax4.set_title("Strand Probability")
	ax4.hist(pisA, alpha=0.5, bins=100, normed=1, label="promoter")
	ax4.hist(pisU, alpha=0.5, bins=100, normed=1, label="unannotated")
	ax4.legend(loc=(0.1,0.8))
	ax4.grid()
	plt.show()

def p53_differences_test(LS):
	sigmasA 	= [math.log(M.models[0].si, 10) for L in LS for chrom in L for st,sp, M in L[chrom] if M.p53_site and   not  M.annotated]
	sigmasU 	= [math.log(M.models[0].si, 10) for L in LS for chrom in L for st,sp, M in L[chrom] if not M.p53_site and  not   M.annotated]
	lamsA 		= [math.log(M.models[0].lam, 10) for L in LS for chrom in L for st,sp, M in L[chrom] if M.p53_site and   not M.annotated]
	lamsU 		= [math.log(M.models[0].lam, 10) for L in LS for chrom in L for st,sp, M in L[chrom] if not M.p53_site and  not   M.annotated]
	weightsA 	= [ M.models[0].wEM  for L in LS for chrom in L for st,sp, M in L[chrom] if M.p53_site and  not  M.annotated]
	weightsU 	= [ M.models[0].wEM  for L in LS for chrom in L for st,sp, M in L[chrom] if not M.p53_site and  not  M.annotated]
	pisA 		= [ M.models[0].pi  for L in LS for chrom in L for st,sp, M in L[chrom] if M.p53_site and  not M.annotated]
	pisU 		= [ M.models[0].pi  for L in LS for chrom in L for st,sp, M in L[chrom] if not M.p53_site and  not  M.annotated]
	
	print len(sigmasA)

	F 			= plt.figure(figsize=(15,10))
	ax1 		= F.add_subplot(2,2,1)
	ax1.set_title("Variance in Loading")
	ax1.hist(sigmasA, alpha=0.5, bins=100, normed=1, label="p53_site")
	ax1.hist(sigmasU, alpha=0.5, bins=100, normed=1, label="no p53_site")
	ax1.legend(loc=(0.1,0.8))
	ax1.grid()
	
	ax2 		= F.add_subplot(2,2,2)
	ax2.set_title("Initiating Length")
	ax2.hist(lamsA, alpha=0.5, bins=100, normed=1, label="p53_site")
	ax2.hist(lamsU, alpha=0.5, bins=100, normed=1, label="no p53_site")
	ax2.legend(loc=(0.1,0.8))
	ax2.grid()

	ax3 		= F.add_subplot(2,2,3)
	ax3.set_title("Paused Probability")
	ax3.hist(weightsA, alpha=0.5, bins=100, normed=1, label="p53_site")
	ax3.hist(weightsU, alpha=0.5, bins=100, normed=1, label="no p53_site")
	ax3.legend(loc=(0.1,0.8))
	ax3.grid()


	ax4 		= F.add_subplot(2,2,4)
	ax4.set_title("Strand Probability")
	ax4.hist(pisA, alpha=0.5, bins=100, normed=1, label="p53_site")
	ax4.hist(pisU, alpha=0.5, bins=100, normed=1, label="no p53_site")
	ax4.legend(loc=(0.1,0.8))
	ax4.grid()
	plt.show()


def label_p53(overlaps,attr="si",LOG=False):
	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left+width+0.02
	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]
	F 			= plt.figure(figsize=(15,10))
	
	axScatter = plt.axes(rect_scatter)
	axHistx = plt.axes(rect_histx)
	axHisty = plt.axes(rect_histy)


	x,y 		= [math.log(getattr(x, attr) ,10) if LOG else getattr(x, attr) for x,y in overlaps],[math.log(getattr(y, attr) ,10) if LOG else getattr(y, attr) for x,y in overlaps]
	
	xy = np.vstack([x,y])
	
	z = gaussian_kde(xy)(xy)
	axScatter.scatter(x, y, c=z, s=14, edgecolor='')
	axScatter.grid()
	
	axHistx.hist(x, bins=200)
	axHisty.hist(y, bins=200, orientation='horizontal')

	axHistx.set_xlim( axScatter.get_xlim() )
	axHisty.set_ylim( axScatter.get_ylim() )
	axHistx.grid()
	axHisty.grid()
	
	axHistx.set_xticks([])
	axHisty.set_yticks([])
	
	x,y 		=  [x for x,y in overlaps if x.p53_site and y.p53_site  and not x.annotated and not y.annotated],[y  for x,y in overlaps if x.p53_site and y.p53_site and not x.annotated and not y.annotated]
	x,y 		= [math.log(getattr(X, attr) ,10) if LOG else getattr(X, attr) for X in x],[ math.log(getattr(X, attr) ,10) if LOG else getattr(X, attr) for X in y]
	axScatter.scatter(x,y,s=20,c="black")

	plt.show()

def p53_binding(Nutlin2_3, DMSO2_3, overlaps):

	#get sets overlaping and not overlaping
	Nutlin2_3, DMSO2_3 	= overlaps_and_not(Nutlin2_3, DMSO2_3)
	#those that are in Nultin2_3 and not in DMS02_3
	Nutlin_N_no_overlap 				= len([1 for chrom in Nutlin2_3  for st,sp, m in Nutlin2_3[chrom] if not m.overlap and not m.annotated ])
	Nutlin_N_overlap 					= len([1 for chrom in Nutlin2_3  for st,sp, m in Nutlin2_3[chrom] if  m.overlap and not m.annotated])
	
	DMSO_N_no_overlap 				= len([1 for chrom in DMSO2_3  for st,sp, m in DMSO2_3[chrom] if not m.overlap and not m.annotated ])
	DMSO_N_overlap 					= len([1 for chrom in DMSO2_3  for st,sp, m in DMSO2_3[chrom] if  m.overlap and not m.annotated ])
	
	Nutlin_p53_no_overlap 		= [m for chrom in Nutlin2_3 for st,sp, m in Nutlin2_3[chrom] if m.p53_site and not m.overlap and not m.annotated]
	Nutlin_no_p53_no_overlap 	= [m for chrom in Nutlin2_3 for st,sp, m in Nutlin2_3[chrom] if not m.p53_site and not m.overlap and not m.annotated ]
	
	Nutlin_p53_overlap 			= [m for chrom in Nutlin2_3 for st,sp, m in Nutlin2_3[chrom] if m.p53_site and m.overlap and not m.annotated ]
	Nutlin_not_p53_overlap 		= [m for chrom in Nutlin2_3 for st,sp, m in Nutlin2_3[chrom] if not m.p53_site and m.overlap and not m.annotated ]

	DMSO2_3_p53_no_overlap 		= [m for chrom in DMSO2_3 for st,sp, m in DMSO2_3[chrom] if m.p53_site and not m.overlap and not m.annotated ]
	DMSO2_3_no_p53_no_overlap 	= [m for chrom in DMSO2_3 for st,sp, m in DMSO2_3[chrom] if not m.p53_site and not m.overlap and not m.annotated ]
	
	DMSO2_3_p53_overlap 		= [m for chrom in DMSO2_3 for st,sp, m in DMSO2_3[chrom] if m.p53_site and m.overlap and not m.annotated ]
	DMSO2_3_no_p53_overlap 		= [m for chrom in DMSO2_3 for st,sp, m in DMSO2_3[chrom] if not m.p53_site and m.overlap and not m.annotated ]


	print "-------------------------"
	print "Of the", Nutlin_N_overlap, "bidirectional calls that are in both DMSO2_3 and Nutlin2_3", float(len(Nutlin_p53_overlap)) /Nutlin_N_overlap , "percent of them overlap a p53 call"
	print "..."
	print "Of the", Nutlin_N_no_overlap, "bidirectional calls that are only in  Nutlin2_3", float(len(Nutlin_p53_no_overlap)) / Nutlin_N_no_overlap  ,"percent of them overlap a p53 call"
	print "..."
	print "Of the", DMSO_N_no_overlap, "bidirectional calls that are only in  DMSO2_3", float(len(DMSO2_3_p53_no_overlap)) /DMSO_N_no_overlap  , "percent of them overlap a p53 call"



	pass

def parameters_dist(L):
	si 	= [m.si for chrom in L for st, sp, S in L[chrom] for m in S.models if m.wEM > 0.5]
	lam = [m.lam for chrom in L for st, sp, S in L[chrom] for m in S.models if m.wEM > 0.5]
	F 	= plt.figure(figsize=(8,6))
	ax1 = F.add_subplot(2,2,1)
	ax1.hist(si,bins=150)

	ax2 = F.add_subplot(2,2,2)
	ax2.hist(lam,bins=150)
	
	plt.show()


