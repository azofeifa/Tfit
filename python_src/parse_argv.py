import os
def checkFile(f):
	vl, text 	= os.path.isfile(f), ""

	if not vl:
		text 	= " does not exist..."
	return vl, text

def checkNumber(vl):
	text 	= ""
	
	try:
		a 	= float(vl)
		b 	= int(vl)
	except:
		text = " is not a number..."
		return False, text
	return True,text
	
def checkDir(path):
	vl, text 	= os.path.isdir(path), ""
	if not vl:
		text 	= " is not a valid path..."
	return vl, text


class argv_wrapper:
	def __init__(self):
		self.G_rm 	= {	"-i":("", True),
			"-c":("all", False),
			"-d":("0.", False),
			"-b":("3", False),
			"-r":("300", False),
			"-t":("5", False)
			}
		self.G_FSI 	= {	"-ref": ("", True, checkFile),
						"-ffs": ("", True, checkFile),
						"-rfs": ("", True, checkFile),
						"-fbg": ("", True, checkFile), 
						"-rbg": ("", True, checkFile), 
						"-wo" : ("", True, checkDir),
						"-pad": ("0", False, checkNumber) }
		self.G_RSO 	= { "-ref": ("", True, checkFile),
						"-fbg": ("", True, checkFile), 
						"-rbg": ("", True, checkFile), 
						"-wo" : ("", True, checkDir),
						"-pad": ("0", False, checkNumber)
						}
		self.G_FSO 	= {	
						"-ffs": ("", True, checkFile),
						"-rfs": ("", True, checkFile),
						"-fbg": ("", True, checkFile), 
						"-rbg": ("", True, checkFile), 
						"-wo" : ("", True, checkDir),
						"-pad": ("0", False, checkNumber) }
		

		self.modules = ["formatData", "runModel"]
		self.modules2= ["FStitchSingleIsoform", "RefSeqOnly", "FStitchMerged"]
		self.error_1 = "please specify the correct module, either: formatData or runModel"
		self.error_2 = "please specify the corrent formating type: " + ",".join(self.modules2)
		self.error_3 = "please specify correct path to RefSeqFile"
		self.error_4 = "please specify correct path to FStitch Forward Strand File"
		self.error_5 = "please specify correct path to FStitch Reverse Strand File"
		self.error_6 = "please specify correct path to Bed/BedGraph File Forward Strand"
		self.error_7 = "please specify correct path to Bed/BedGraph File Reverse Strand"
		self.error_8 = "please specify a path to write_out File"

		self.mod, self.mod2 		= "", ""
		self.G 						= None
		self.errors 				= list()
		self.exit 					= False
		self.i 						= 2
		self.help 					= ["-h", "--help"]
	def checkForModuleType(self,argv):

		if len(argv)<2:
			self.exit=True
			self.errors.append(self.error_1)
		else:
			mod 	= argv[1]
			if mod in self.modules:
				self.mod 	= mod 
				self.i 		= 2
				self.G 		= self.G_rm
			else:
				self.exit=True
				self.errors.append(self.error_1)
			if self.mod == "formatData":
				if len(argv)<3:
					self.exit=True
					self.errors.append(self.error_2)
				else:
					mod2 = argv[2]
					if mod2 in self.modules2:
						self.mod2 	= mod2

						self.i 		= 3
						if mod2 == self.modules2[0]:
							self.G 		= self.G_FSI
						elif mod2 == self.modules2[1]:
							self.G 		= self.G_RSO
						elif mod2 == self.modules2[2]:
							self.G 		= self.G_FSO

					else:
						self.exit=True
						self.errors.append(self.error_2)
	def checkForHelp(self, argv):
		pass

	def printErrors(self):
		for error in self.errors:
			print error



	def add(self, argv):
		current, add 	= None, None
		errors 			= list()
		for arg in argv[self.i:]:
			if "-"==arg[0] and arg not in self.G.keys() :
				errors.append(arg)
			elif "-"==arg[0]:
				current = arg
			elif current is not None:
				self.G[current] 	= (arg,self.G[current][1], self.G[current][2])
				current 			= None
	def checkImportantArgs(self):
		for i,(vl, important, f) in zip(self.G.keys(), self.G.values()):
			if important and not vl:
				self.errors.append("Did not specify value for " + str(i))
				self.exit 	= True
			if not f(vl)[0]:
				error, text 	= f(vl)
				self.errors.append(str(vl) + text)
				self.exit 		= True
		pass



def run(argv):
	aw 	= argv_wrapper()
	aw.checkForModuleType(argv)
	if aw.exit:
		return aw
	aw.add( argv)
	aw.checkImportantArgs()
	return aw
