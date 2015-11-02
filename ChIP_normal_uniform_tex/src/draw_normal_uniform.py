import matplotlib.pyplot as plt
import numpy as np 

N=[x for x in np.random.normal(0,1,1000) ]
U=[x for x in np.random.uniform(-10,10, 1000)]
D=N+U
plt.title(r'$\mu=0, \sigma^2=1, \pi=0.5$')
plt.hist(D, bins=100,alpha=0.5)
plt.xlabel("relative genomic coordinates")
plt.ylabel("read coverage/frequency")
plt.grid()
plt.show()

