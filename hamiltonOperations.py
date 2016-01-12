import graphs
import fileOperations
CELL_INSIDE = 0
CELL_OUTSIDE = 1 
class HamiltonFun:
	def __init__(self, lattice_rows, lattice_cols):
		""" initializes a graph object """
		self.__rows = lattice_rows
		self.__cols = lattice_cols
		self.__latticeGraphDict = self.latticeGraphDictGenerator()
		self.__latticeGraph = graphs.Graph(self.__latticeGraphDict)
		self.__initCycle = []
		
		self.__activePath = []
		self.__activeCells = []
		self.__splitPoints = []
		self.__neighbourCycles = []
		self.__splitCells = []
		
	def find_cycle_split_positions(cycle):
		pass
		
	def is_hamilton_cycle_existing(self):
		#only cycle possible if even number of nodes in lattice.
		return self.__rows * self.__cols % 2 == 0
		
	def set_up_first_cycle(self, cycle = None):
		#finding a cycle can take an enormous amount of time for bigger lattices. So, provide one. Or, make a simple algorithm to find a snake...
		if cycle is None:
			self.__initCycle = self.find_hamilton_cycle()
		else:
			self.__initCycle = cycle
		
		if not self.is_hamilton_cycle(self.__initCycle):
			print "ASSERT ERROR: no hamilton cycle found"
			self.print_path_ASCII(self.__initCycle) 
		return self.__initCycle
	
	def create_cell_pattern_from_hamilton_cycle(self, path):
		#if a cycle is drawn, there is an inside and an outside. 
		#each four points of the lattice define a cell. This cell is in or outside the cycle.
		
		#ASSERT hamilton cycle
		
		#default all cells are outside path
		cells = {(r,c):CELL_OUTSIDE for r in range(self.__rows-1) for c in range(self.__cols-1)}
		return self.__explore_inside_of__hamilton_cycle(path, cells, (0,0)) 
	
	def __explore_inside_of__hamilton_cycle (self, path, cells, nodeCell):
		#every cell has four directions, N, E, S, W . inside cells are all connected. Go over inside recursively
		#add nodeCell as inside in cells
		cells[nodeCell] = CELL_INSIDE
		# print path
		for drA, dcA,drB,dcB,drC,dcC in [(0,0,0,1,-1,0),(0,1,1,1,0,1),(1,1,1,0,1,0),(1,0,0,0,0,-1)]:
			nodeA = (nodeCell[0]+drA,nodeCell[1]+dcA)
			nodeB = (nodeCell[0]+drB,nodeCell[1]+dcB)
			nextNodeCell = (nodeCell[0]+drC, nodeCell[1]+dcC) #coordinate from adjecent cell. in direction indicated by node a and b
			try:
				if cells[nextNodeCell]==CELL_INSIDE:
					#if next nodecell already scanned, don't bother checking it.
					continue
			except:
				#if nextNodeCell is an invalid position, dont bother checking that direction...
				continue
				
			#check directions  nodeCell
			# try:
			pathNeighboursNodeA = self.neighbours_on_hamilton_path(path, nodeA)
			# pathNeighboursNodeB = self.neighbours_on_hamilton_path(path, nodeB)
			
			# print pathNeighboursNodeA
			# print nodeB
			if nodeB not in pathNeighboursNodeA: #	or nodeA in pathNeighB
				self.__explore_inside_of__hamilton_cycle(path,cells, nextNodeCell )
			# except:
				# #assert graph boundaries passed...
				# # print "boundary"
				# # print nodeCell
				# pass
		
		return cells
			#for every open direction, recursive.
		
	# def split(self, path, splitCell):
		

		
	def find_split_cells(self, path):
		#find all potential split cells in path or cycle.  
		#  cell (0,0) with its corner nodes:
		#    (0,0)     (0,1)  ... 
		#     
		#    (1,0)     (1,1)  ...  
		
		# at node 0,0  we investigate the drawn situation.  if there is one edge, not good, if there are three edges, not good. if there are 2 edges that don't touch each other (two vertical OR two horizontal) OK! (add to list)
		#this is a split point, because we can transform the two horizontal edges to vertical ones (or the other way around)  and have a new path.
		
		#define inside and outside of path
		cells = self.create_cell_pattern_from_hamilton_cycle(path)
		
		#we run through each node in the lattice and check it with the three nodes east, south-east and south if the situation is occuring. if so, add node to list as "found"
		#we could say we run through each cell in the lattice. Only bother if the cell is inside the cycle, otherwise, no swapping possible. (the cycle would be split!)
		splitCells = []
		for row in range(self.__rows-1):
			for col in range(self.__cols-1):
				if cells[(row,col)] == CELL_OUTSIDE:
					continue
				#nodes to examine in Lattice, East, SE, South  (referenced from the top left node, which is the same coord as the cell name)
				nodesToCheckInLattice = {(row,col+1),(row+1,col),(row+1,col+1)} #create as set
				
				#orthogonal neighbours of first node E and S
				orthoNodesInLattice =  {(row,col+1),(row+1,col)} #create as set
				
				#get neighbour nodes on path on top left node.
				neighboursOnPathTopLeftNode =  self.neighbours_on_hamilton_path(path,(row,col))
				
				#------------------------
				#there should be exactly one neighbour in the nodes to examine (E or S)
				nodesPartOfTopLeftNodePath = neighboursOnPathTopLeftNode.intersection(orthoNodesInLattice)
				
				if len(nodesPartOfTopLeftNodePath) != 1:
					#isolated if zero, corner if 2
					continue #next iteration of for loop (is not the same as break!!!)
				else:
					topLeftNodeConnectedNode = nodesPartOfTopLeftNodePath.pop()
				#------------------------
				#get neighbour nodes on path on bottom right node.
				neighboursOnPathBottomRightNode =  self.neighbours_on_hamilton_path(path,(row+1,col+1))
				
				#there should be exactly one neighbour in the nodes to examine (E or S)
				nodesPartOfBottomRightNodePath = neighboursOnPathBottomRightNode.intersection(orthoNodesInLattice)
				
				if len(nodesPartOfBottomRightNodePath) != 1:
					#isolated if zero, corner if 2
					continue #next iteration of for loop (is not the same as break!!!)
				else:
					bottomRighNodeConnectedNode = nodesPartOfBottomRightNodePath.pop()
				
				#------------------------
				if bottomRighNodeConnectedNode == topLeftNodeConnectedNode:
					#u shape (all nodes connected)
					continue #next iteration of for loop (is not the same as break!!!)
				else:
					splitCells.append((row,col))
		return splitCells
		
	def __split_search_return_one_neighbour_or_none(self,path, node):
		pass
	
	def neighbours_on_hamilton_path(self,path, node):
		#assert path is hamilton path or cycle
		
		#get index of node in path
		pos = path.index(node)
		
		neighboursOnPath = set()
		if (pos>0):
			neighboursOnPath.add(path[pos-1])
		else:
			#check if cycle if, so, add element from end of path
			#ASSERT pos == 0
			if path[0]==path[-1]:
				neighboursOnPath.add(path[-2]) #penultimate element of list
		try:
			neighboursOnPath.add(path[pos+1])
		except IndexError:
			#ASSERT pos == len(path)
			#check if cycle if, so, add element from end of path
			if path[0]==path[-1]:
				neighboursOnPath.add(path[1]) #second element of list
		return neighboursOnPath
	
		
	
	def latticeGraphDictGenerator(self):
		#of this type, where the node name is a coordinate: (row,col)
		# g = { "a" : ["d","f"],
		   # "b" : ["c","b"],
		   # "c" : ["b", "c", "d", "e"],
		   # "d" : ["a", "c"],
		   # "e" : ["c"],
		   # "f" : ["a"]
		# }

		#lattice: every node is connected to a neighbour
		latticeGraph = {}
		for row in range(self.__rows):
			for col in range(self.__cols):
				neighbours = []
				north = row -1
				if north >= 0:
					neighbours.append((north,col))
				south = row +1
				if south < self.__rows:
					neighbours.append((south,col))
				east = col + 1
				if east < self.__cols:
					neighbours.append((row,east))
				west = col - 1
				if west >= 0:
					neighbours.append((row,west))
				latticeGraph[(row,col)]=neighbours
		return latticeGraph	

	def find_path():
		print self.__latticeGraph.find_path((0,0),(0,1))
	
	def find_hamilton_cycle(self):
		cycle = self.find_hamilton_path((0,0),(0,1))
		cycle.append(cycle[0]) #close the loop
		return cycle
	
	def find_all_hamilton_cycles(self):
		'''closed hamilton path
		'''
		paths = self.__latticeGraph.find_all_paths((0,0),(0,1)) 
		self.__hamiltonCycles= []
		for path in paths:
			if len(path) == len(self.__latticeGraph.vertices()):
				path.append(path[0]) #close the path
				self.__hamiltonCycles.append( path)
		return self.__hamiltonCycles
	
	
	def find_all_hamilton_paths_topLeftToBottomRight(self):
		return self.find_all_hamilton_paths((0,0),(self.__rows-1,self.__cols-1))
		
	def find_hamilton_path_topLeftToBottomRight(self):
		return self.find_hamilton_path((0,0),(self.__rows-1,self.__cols-1))
		
	def find_hamilton_path(self,startVertex, endVertex):
		path = []
		while len(path) != len(self.__latticeGraph.vertices()):
			path = self.__latticeGraph.find_hamilton_path(startVertex, endVertex) 
		return path
		
	def find_all_hamilton_paths(self,startVertex, endVertex):
		'''paths going through each vertex only once, begin and end position set fixed.
		'''
		paths = self.__latticeGraph.find_all_paths(startVertex, endVertex) 
		self.__hamiltonPaths= []
		
		for path in paths:
			if len(path) == len(self.__latticeGraph.vertices()):
				self.__hamiltonPaths.append( path)
		return self.__hamiltonPaths
	
	def is_hamilton_path(self,path):
		print "not yet implemented"
		pass
		
	def is_hamilton_cycle(self,cycle):
		print "not yet implemented"
		pass
	
	def print_cells_ASCII(self, cells):
		printcells = ""
		for row in range(self.__rows - 1):
			printrow = ""
			for col in range(self.__cols - 1):
				if cells[(row,col)] ==1:
					printrow += "X"
				else:
					printrow += " "
			printcells += printrow + "\n"
		print printcells
		
	def print_path_ASCII(self,path):
		#from a list of vertices, and rows and cols, print path on screen
		
		#create row with points on each node
		print_empty_row = self.__cols*"X"
		print_empty_row = list(" ".join(print_empty_row))
		
		#create completely blank row
		blankRowPrint =  list(((2*self.__cols - 1)*" "))
		
		#add all rows to lattice
		latticeMinimalCoords = []
		for row in range(self.__rows):
			latticeMinimalCoords.append(print_empty_row[:])
			latticeMinimalCoords.append(blankRowPrint[:])
		latticeMinimalCoords.pop()
		
		#double all node(row,col)names to coords of path
		path = [(2*r,2*c) for r,c in path]
			
		#fill lattice with path
		previousNode = path[0]
		for node in path[1:]:
		
			rowCoord = (previousNode[0] + node[0])/2
			colCoord = (previousNode[1] + node[1])/2
			latticeMinimalCoords[rowCoord][colCoord] = "X"
			previousNode = node
			
		#print empty lattice
		for printrow in latticeMinimalCoords:
			print "".join(printrow)

			
def print_all_hamilton_paths(rows, cols, startNode=None, endNode=None):
	#print all hamilton paths for a given start and ending
	lattice = HamiltonFun(rows,cols)
	if startNode is None or endNode is None:
		paths = lattice.find_all_hamilton_paths_topLeftToBottomRight()
	else:
		paths = lattice.find_all_hamilton_paths(startNode, endNode)
	print_all_paths(lattice, paths)

def print_all_hamilton_cycles_inefficient(rows, cols):
	#print all hamilton cycles
	lattice = HamiltonFun(rows,cols)
	cycles = lattice.find_all_hamilton_cycles()
	print_all_paths(lattice, cycles)

def print_all_paths(lattice, paths):
	for path in paths:
		lattice.print_path_ASCII(path)
		print "\n"
		
if __name__ == "__main__":
	ROWS = 5
	COLS = 6
	
	lattice = HamiltonFun(ROWS, COLS)
	
	# dictje = lattice.latticeGraphDictGenerator()
	# hoi = graphs.Graph(dictje)
	# print hoi
	# print hoi.vertices()
	
	
	print lattice.is_hamilton_cycle_existing()
	path = lattice.set_up_first_cycle()
	print lattice.find_split_cells(path)
	# cells = lattice.create_cell_pattern_from_hamilton_cycle(path)
	#print_all_hamilton_paths(ROWS, COLS,(0,0),(2,3))
	# print_all_hamilton_cycles_inefficient(ROWS, COLS)
	# lattice.print_cells_ASCII(cells)