def run(argv):
	G 	= {	"-i":"",
			"-c":"all",
			"-d":"0.",
			"-b":"3",
			"-r":"300",
			"-t":"5"
			}
	add 	= False
	current = None
	for i in argv:
		if i in argv and not add:
			add=True
			current=i
		elif current is not None:
			G[current] 	= i
			add 		= False
	return G