import time
class info:
	def __init__(self, start, stop, info=None):
		self.start 	= start
		self.stop 	= stop
		self.info 	= info

class node:
	def __init__(self,start, stop):
		self.start 				= start
		self.stop 				= stop
		self.intervals 			= list()
		self.right, self.left 	= None, None
	def search(self, interval):
		finds 					= list()
		for info in self.intervals:
			if info[0] <= interval[0] <= info[1] or info[0] <= interval[1] <=info[1] or (interval[0] <= info[0] and interval[1] >= info[1]) :
				finds.append(info)
		return finds
	def __str__(self):
		return str(self.start)+"-"+str(self.stop) + ", " + str(len(self.intervals))
	def erase_repeats(self):
		merged=True
		i 		= 0
		while  merged:
			N 		= len(self.intervals)
			merged 	= False
			i 		= 0
			while (i+1 < N and 
				self.intervals[i][0]==self.intervals[i+1][0] and 
				self.intervals[i][1]==self.intervals[i+1][1] ):
				self.intervals 	= self.intervals[:i]+self.intervals[i+1:]
				N 				= len(self.intervals)
				merged=True
				break
		pass
class treeNode:
	def __init__(self):
		self.node 	= None
		self.right 	= None
		self.left 	= None
	def build(self, nodes):
		i 			= len(nodes)/ 2
		self.node 	= nodes[i]
		if i > 0:
			self.left 	= treeNode()
			self.left.build(nodes[:i])
		if i+1 < len(nodes):
			self.right 	= treeNode()
			self.right.build(nodes[i+1:])
	def searchInterval(self, interval):
		start, stop 	= interval
		if self.node.start <= start <= self.node.stop or self.node.start <= stop <= self.node.stop or (start <= self.node.start and stop >= self.node.stop):
			return self.node.search(interval)
		if stop < self.node.start and self.left:
			return self.left.searchInterval(interval)
		if start > self.node.stop and self .right:
			return self.right.searchInterval(interval)
		return None

	def searchPoint(self, point):
		if self.node.start <= point <= self.node.stop:
			return self.node
		if point < self.node.start and self.left:
			return self.left.searchPoint(point)
		if point > self.node.stop and self .right:
			return self.right.searchPoint(point)
		return None


class tree:
	def __init__(self, *args):
		assert len(args) < 2, "either no or 1 argumnet"
		self.root 	 	= None
		if len(args)==1:
			self.build(args[0])
	def assemble(self, LST):
		nodes 	= list()
		while LST:
			i 	= 0
			N 	= len(LST)
			o_st,o_sp 	= LST[i][0],LST[i][1]
		
			while i < N and (LST[i][0] < o_sp and LST[i][1]> o_st) : #want to find where there are no overlaps split on that
				o_st, o_sp 	= min((o_st, LST[i][0])),max((o_sp, LST[i][1]))
				i+=1
			left, right 	=  LST[:i], LST[i:]
			x,y 			= [info[0] for info in left],[info[1] for info in left]
			NODE 			= node(min(x), max(y))
			NODE.intervals 	= left		
			NODE.erase_repeats()
			nodes.append(NODE)
			LST 			= right
		return nodes
			
	def build(self, LST):
		LST.sort()
		nodes 	= self.assemble(LST)
		#root is the middle of nodes
		i 			= len(nodes)/ 2
		self.root 	= treeNode()
		self.root.build(nodes)
	def searchInterval(self, interval):
		assert self.root, "interval tree has not been built yet"
		return self.root.searchInterval(interval)
	def searchPoint(self,point):
		assert self.root, "interval tree has not been built yet"
		return self.root.searchPoint(point)



if __name__ == "__main__":
	LST 	= [(1,4, True), (3,5, True), (0,7,True), (9,12, True), (11,14, True), (10,15, True), (20, 25, True), (24, 26, True)]
	T 		= tree(LST)
	finds 	= T.searchInterval( (2,4) )
	print finds
