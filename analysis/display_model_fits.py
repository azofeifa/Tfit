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


