def load(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			u,v,w 	= line.strip("\n").split(",")
			w 		= float(w)
			if u not in G:
				G[u] = {}
			if v not in G[u]:
				G[u][v]=w
			if v not in G:
				G[v] ={}
			if u not in G[v]:
				G[v][u]=w
	return G
def convert(G, OUT):
	FHW 	= open(OUT, "w")
	nodes 		= dict([(i,n)for i,n in enumerate(G.keys())])
	other_way 	= dict([(n,i)for i,n in enumerate(G.keys())])
	FHW.write("nodedef>name VARCHAR,label VARCHAR\n")
	for n in nodes:
		FHW.write("s" + str(n) + "," + nodes[n]+ "\n")
	FHW.write("edgedef>node1 VARCHAR,node2 VARCHAR\n")
	for u in G:
		for v in G[u]:
			if G[u][v] <= 2:
				
				FHW.write("s" + str(other_way[u])+ "," +"s" +  str(other_way[v])+"\n")









if __name__ =="__main__":
	FILE 	="/Users/joazofeifa/Lab/TF_predictions/adjacency_matrices/TF_network.csv"
	OUT 	= "/Users/joazofeifa/Lab/TF_predictions/adjacency_matrices/TF_network.GDF"
	G 		= load(FILE)
	convert(G, OUT)
