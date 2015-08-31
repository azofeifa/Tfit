import simulate, model
import numpy as np
import matplotlib.pyplot as plt
import load

if __name__=="__main__":
	# D 	= simulate.runOne(mu=0, s=1, l=5, lr=100, ll=-100, we=0.5,wl=0.25, 
	# 	wr=0.25, pie=0.5, pil=0.1, pir=0.9, N=1000, SHOW=False , bins=200, noise=True )
	X 	= load.grab_specific_region("chr3", 5709069,6479278, SHOW=True, bins=300 )
	