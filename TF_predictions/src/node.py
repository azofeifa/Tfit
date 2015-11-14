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
	def searchInterval(self, interval):
		finds 					= list()
		for info in self.intervals:
			if info[0] < interval[1] and info[1] > interval[0]:# <=info[1] or info[0] <= interval[1] <=info[1] or (interval[0] <= info[0] and interval[1] >= info[1]) :
				finds.append(info)
		return finds
	def searchPoint(self, pt):
		finds 					= list()
		for info in self.intervals:
			if info[0]< pt < info[1]:
				finds.append((info))
		return finds
	def __str__(self):
		return str(self.start)+"-"+str(self.stop) + ", " + str(len(self.intervals))
class treeNode:
	def __init__(self):
		self.node 	= None
		self.right 	= None
		self.left 	= None
		self.parent = None
		self.max 	= None
	def __str__(self):
		return str(self.node.start) + "-" + str(self.node.stop) + ", max: " + str(self.max) + ", min: " + str(self.min)
	def build(self, nodes):
		i 				= min(len(nodes)/ 2, len(nodes)-1)
		self.node 		= nodes[i]
		if i > 0:

			if len(nodes[:i]) > 0:
				self.left 			= treeNode()
				self.left.parent 	= self
				self.left.build(nodes[:i])
		if i+1 < len(nodes):
			if len(nodes[i+1:])>0:
				self.right 	= treeNode()
				self.right.parent 	= self

				self.right.build(nodes[i+1:])
	def get_max(self):
		if self.right:
			return max(self.node.stop, self.right.get_max() )
		return self.node.stop
	def get_min(self):
		if self.left:
			return min(self.node.start, self.left.get_min() )
		return self.node.start
		
		
	def searchInterval(self, interval):
		if interval[0] > self.max:
			return []
		if interval[1] < self.min:
			return []
		if self.left and self.left.max > interval[0]:
			return self.node.searchInterval(interval) + self.left.searchInterval(interval)
		if self.right:
			return self.node.searchInterval(interval) + self.right.searchInterval(interval)
		return self.node.searchInterval(interval)
	def searchPoint(self, point):
		if self.node.start <= point <= self.node.stop:
			return self.node.searchPoint(point)
		if point < self.node.start and self.left:
			return self.left.searchPoint(point)
		if point > self.node.stop and self .right:
			return self.right.searchPoint(point)
		return []
	def get_all(self):
		if self.left and self.right:
			return self.node.intervals+self.left.get_all() + self.right.get_all()
		if self.left:
			return self.node.intervals+ self.left.get_all()
		if self.node.right:
			return self.node.intervals+ self.right.get_all()
		return self.node.intervals
	def set_max(self):
		
		self.max 	= self.get_max()
		if self.right:
			self.right.set_max()
		if self.left:
			self.left.set_max()
	def set_min(self):
		self.min 	= self.get_min()
		if self.right:
			self.right.set_min()
		if self.left:
			self.left.set_min()
	

class tree:
	def __init__(self, *args):
		assert len(args) < 2, "either no or 1 argumnet"
		self.root 	 	= None
		self.SEARCH 	= False
		if len(args)==1:
			self.build(args[0])
		if self.SEARCH:
			self.root.set_max()
			self.root.set_min()
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
			nodes.append(NODE)
			LST 			= right
		return nodes
	def get_all(self):
		if self.SEARCH:
			return self.root.get_all()
		return []
	def build(self, LST):
		if len(LST):
			LST.sort()
			nodes 	= self.assemble(LST)
			#root is the middle of nodes
			i 			= len(nodes)/ 2
			self.root 	= treeNode()
			self.root.build(nodes)
			self.SEARCH = True
	def searchInterval(self, interval):

		if self.SEARCH:
			return self.root.searchInterval(interval)
		else:
			return []
	def searchPoint(self,point):
		if self.SEARCH:
			return self.root.searchPoint(point)
		else:
			return []



if __name__ == "__main__":
	LST 	= [(1,4, True), (7,10, True)  , (15,25), (17,26) , (27,35)  ]
	T 		= tree(LST)
	#====================================
	#points
	finds 	= T.searchPoint( 2 )
	print finds
	finds 	= T.searchPoint( 8 )
	print finds
	finds 	= T.searchPoint( 20 )
	print finds
	finds 	= T.searchPoint( 5 )
	print finds
	finds 	= T.searchPoint( 38 )
	print finds


	#====================================
	finds 	= T.searchInterval( (1,3) )
	print finds
	finds 	= T.searchInterval( (8,9) )
	print finds
	finds 	= T.searchInterval( (2,8) )
	print finds
	finds 	= T.searchInterval( (22,38) )
	print finds
	finds 	= T.searchInterval( (38,45) )
	print finds




