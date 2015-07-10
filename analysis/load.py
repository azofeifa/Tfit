class EMG:
	def __init__(self, mu, si, l, w, pi):
		self.mu, self.si, self.l 	= mu, si, l
		self.w, self.pi 			= w,pi
		self.type 					= "N"
class UNI:
	def __init__(self, a,b, w, pi):
		self.a, self.b 	= a,b
		self.w, self.pi = w,pi
		self.type 		= "U"

class model:
	def __init__(self, ll , converged, retries):
		self.ll 		= ll
		self.converged 	= converged
		self.retries 	= retries
		self.rvs 		= list()

class info:
	def __init__(self, chrom, start, stop, N):
		self.chrom 	= chrom
		self.start 	= start
		self.stop 	= stop
		self.models = {}
		self.current= None
		self.N 		= N
	def add_model(self, line):
		k,ll, converged, retries 	= line[1:].strip("\n").split(",")
		k,ll, retries 				= int(k), float(ll), float(retries)
		converged 					= bool(converged=="True")
		assert k not in self.models, "there should only be one model entry"
		self.models[k] 	= model(ll , converged, retries)
		self.current 	= k
	def add_components(self, line):
		TYPE, Info 	= line.strip("\n").split(": ")
		if TYPE=="N":
			mu, si, l, w, pi 	= Info.strip("\n").split(",")
			mu, si, l, w, pi 	= float(mu), float(si), float(l),float(w), float(pi)
			self.models[self.current].rvs.append(EMG(mu, si, l, w, pi))
		elif TYPE=="U":
			a,b,w, pi 	= Info.strip("\n").split(",")
			a,b,w, pi 	= float(a), float(b), float(w), float(pi)
			self.models[self.current].rvs.append(UNI(a, b, w, pi))

def EMG_out(FILE):
	I 		= None
	fits 	= list() 	
	with open(FILE) as FH:
		for line in FH:
			if "#"==line[0]:
				chrom, stsp,N 	= line[1:].strip("\n").split(":")
				start, stop 	= stsp.split("-")
				start, stop 	= int(start), int(stop)
				if I is not None:
					fits.append(I)
				I 				= info(chrom,start, stop, float(N))
			elif "~" == line[0]:
				I.add_model(line)
			elif I is not None:
				I.add_components(line)
	return fits










