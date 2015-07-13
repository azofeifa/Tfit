import os

def get_np_cpu_walltime(FILE):
	np, CPU, wall 	= None,None,None
	with open(FILE) as FH:
		for line in FH:
			if "-np" in line:

				np=line.split(":")[1]
				np=int(np[1:].strip("\n"))

			elif "CPU" in line:
				info 	= line.split(":")[1]
				ticks, CPU 	= info.split(",")
				CPU 		= float(CPU[1:].strip("\n"))
			elif "Wall" in line:
				wall 	= line.split(":")[1]
				wall 		= float(wall[1:].strip("\n"))
	return np, CPU, wall


def all_files(directory):
	G 	= {}
	for FILE in os.listdir(directory):
		np, CPU, wall 	= get_np_cpu_walltime(directory+FILE)
		if np is not None:
			if np not in G:
				G[np] 	= (list(), list())
			G[np][0].append(CPU)
			G[np][1].append(wall)
	return G


		