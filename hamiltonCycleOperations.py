import latticeOperations

CELL_INSIDE = 0
CELL_OUTSIDE = 1 
CELL_SPLITPOINT_POTENTIAL = 2
CELL_SPLITPOINT = 3
CELL_RECOMBINATION_CANDIDATE = 4

class HamiltonCycle():
	def __init__(self, lattice_rows, lattice_cols, cycle):
		self.__lattice = latticeOperations.Lattice(lattice_rows, lattice_cols)
		
		#let path start from (0,0)  here we don't yet really check for hamilton, but if error, we know it is wrong already
		try:
			cycle = standardized_hamilton_cycle(cycle)
		except:
			print "make sure to provide a single hamilton cycle"
			raise 
			
		self.addCycle(cycle) 
		#from here cycle is asserted to be hamilton
		self.__cells = {}
		self.__cycle = cycle 
		
		self.__splitCells = [] #all the potential splitting points.
		self.__splitPaths = [] #all the potential splitting points.
		
	def __str__(self):
		self.__lattice.string_from_paths()
		return str(self.__lattice)
		
	def addCycle(self,path):
		#add hamilton cycle, and check if it is one
		self.__lattice.new_paths_to_lattice(path)
		stats = self.__lattice.stats()
		if not stats["hamilton"] or stats["cycles"]!=1:
			print self
			raise Exception("no hamilton path provided")
	def get_split_paths(self):
		return self.__splitPaths
		
	def create_cell_pattern_from_hamilton_cycle(self):
		#if a cycle is drawn, there is an inside and an outside. 
		#each four points of the lattice define a cell. This cell is in or outside the cycle.
		
		#default all cells are outside path
		self.__cells = {(r,c):CELL_OUTSIDE for r in range(self.__lattice.rows()-1) for c in range(self.__lattice.cols()-1)}
		try:
			self.__explore_inside_of__hamilton_cycle((0,0)) 
		except:
			print "creation of cells failed"
			raise
	
	def __explore_inside_of__hamilton_cycle (self, nodeCell):
		#every cell has four directions, N, E, S, W . inside cells are all connected. Go over inside recursively
		#add nodeCell as inside in cells
		self.__cells[nodeCell] = CELL_INSIDE
		
		#four directions, (common nodes on cycle with adjecent cell and adjecent cell )
		for drA, dcA,drB,dcB,drC,dcC in [(0,0,0,1,-1,0),(0,1,1,1,0,1),(1,1,1,0,1,0),(1,0,0,0,0,-1)]:
			nodeA = (nodeCell[0]+drA,nodeCell[1]+dcA)
			nodeB = (nodeCell[0]+drB,nodeCell[1]+dcB)
			nextNodeCell = (nodeCell[0]+drC, nodeCell[1]+dcC) #coordinate from adjecent cell. in direction indicated by node a and b
			
			#check if already INSIDE or if position is valid
			try:
				if self.__cells[nextNodeCell]==CELL_INSIDE:
					#if next nodecell already scanned, don't bother checking it.
					continue
			except:
				#if nextNodeCell is an invalid position, dont bother checking that direction...
				continue
				
			#check directions  nodeCell
			pathNeighboursNodeA = self.neighbours_of_node( nodeA)
			
			if nodeB not in pathNeighboursNodeA: #	or nodeA in pathNeighB
				#if there is not path from A to B, then, the cells are connected
				#go explore the neighbour
				self.__explore_inside_of__hamilton_cycle( nextNodeCell )

	def neighbours_of_node(self, node):
		#from old neighbours_on_cycle in first program.
		#get the cycle
		path = self.__cycle 
		
		#get index of node in path
		pos = path.index(node)
		
		neighboursOnPath = set()
		
		#previous node
		if (pos>0):
			neighboursOnPath.add(path[pos-1])
		else:
			#check if cycle if, so, add element from end of path
			neighboursOnPath.add(path[-2]) #penultimate element of list
		
		#next node
		try:
			neighboursOnPath.add(path[pos+1])
		except IndexError:
			#ASSERT pos == len(path)
			neighboursOnPath.add(path[1]) #second element of list
		return neighboursOnPath
		
	def print_cells_ASCII(self):
		printcells = ""
		for row in range(self.__lattice.rows() - 1):
			printrow = ""
			for col in range(self.__lattice.cols() - 1):
				if self.__cells[(row,col)] ==CELL_INSIDE:
					printrow += "X"
				elif self.__cells [(row,col)] == CELL_SPLITPOINT_POTENTIAL:
					printrow += "x"
				elif self.__cells [(row,col)] == CELL_OUTSIDE:
					printrow += " "
				elif self.__cells [(row,col)] == CELL_SPLITPOINT:
					printrow += "+"
				elif self.__cells [(row,col)] == CELL_RECOMBINATION_CANDIDATE:
					printrow += "o"
				else:
					printrow += "?"
			printcells += printrow + "\n"
		print printcells
	
	######################################
	#neighbour finder
	######################################
	
	def find_all_neighbour_hamilton_cycles(self):
		self.find_split_cells()
		self.print_cells_ASCII()
		for cell in self.__splitCells:
			self.split(cell)
	
	def find_split_cells(self):
		#find all potential split cells in path or cycle.  
		#  cell (0,0) with its corner nodes:
		#    (0,0)     (0,1)  ... 
		#     
		#    (1,0)     (1,1)  ...  
		
		# at node 0,0  we investigate the drawn situation.  if there is one edge, not good, if there are three edges, not good. if there are 2 edges that don't touch each other (two vertical OR two horizontal) OK! (add to list)
		#this is a split point, because we can transform the two horizontal edges to vertical ones (or the other way around)  and have a new path.
		
		#we run through each node in the lattice and check it with the three nodes east, south-east and south if the situation is occuring. if so, add node to list as "found"
		#we could say we run through each cell in the lattice. Only bother if the cell is inside the cycle, otherwise, no swapping possible. (the cycle would be split!)
		splitCells = []
		for row in range(self.__lattice.rows()-1 ):
			for col in range(self.__lattice.cols()-1 ):
				if self.__cells[(row,col)] == CELL_OUTSIDE:
					continue
				#nodes to examine in Lattice, East, SE, South  (referenced from the top left node, which is the same coord as the cell name)
				nodesToCheckInLattice = {(row,col+1),(row+1,col),(row+1,col+1)} #create as set
				
				#orthogonal neighbours of first node E and S
				orthoNodesInLattice =  {(row,col+1),(row+1,col)} #create as set
				
				#get neighbour nodes on path on top left node.
				neighboursOnPathTopLeftNode =  self.neighbours_of_node((row,col))
				
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
				neighboursOnPathBottomRightNode =  self.neighbours_of_node((row+1,col+1))
				
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
		#indicate splitcells in cells
		for cell in splitCells:
			self.__cells[cell] = CELL_SPLITPOINT_POTENTIAL 
		
		self.__splitCells = splitCells
	
	
	def split(self, splitCell):
		#splitcell is also a coordinate on the cycle (left up coordinate of cell equals path coordinate)
		#split the path in two at splitcell
		
		# nodeTopLeftPositionInPath = path.index(splitCell)
		# nodeBottomRightPositionInPath = path.index((splitCell[0]+1, splitCell[1]+1))
		
		
		neighbours = self.neighbours_of_node( splitCell)
		
		
		#take the opposite node of splitcell 
		diagonalSplitNode = (splitCell[0]+1,splitCell[1]+1)
		#splitNode is same as splitcell
		
		
		#look if the splitcell has horizontal or vertical neighbours , define the path neighbours of the split node and of the diagonal split node
		isVertical = True #vertical means: 1 inside cell above, and one down of the splitcell.
		splitNodeNoneNeighbourOnPath = (splitCell[0],splitCell[1]+1)
		splitNodeNeighbourOnPathOfSplitcellNode = (splitCell[0]+1,splitCell[1])
		
		if  (splitNodeNoneNeighbourOnPath) in neighbours:
			splitNodeNeighbourOnPathOfSplitcellNode, splitNodeNoneNeighbourOnPath = splitNodeNoneNeighbourOnPath, splitNodeNeighbourOnPathOfSplitcellNode
			isVertical = False
		
		#do the splitting of the cycle at the split node
		twoPaths = split_path_at_two_neighbours(self.__cycle[:],splitCell,splitNodeNeighbourOnPathOfSplitcellNode)
		#if the two nodes are front and end, one of the returning paths will be empty.
		pathA = twoPaths[0]
		pathB = twoPaths[1]
		
		#make sure diagonal is in next path to split.
		if (diagonalSplitNode) not in pathA:
			pathA,pathB = pathB,pathA
		paths = split_path_at_two_neighbours(pathA[:],diagonalSplitNode,splitNodeNoneNeighbourOnPath)
		
		#all three pieces of the original path are collected in paths
		paths.append(pathB)
		
		#sort out the pieces
		startNodePath = None
		endNodePath = None
		middlePath = None
		for path in paths:
			if self.__cycle[0] in path: #start node
				startNodePath = path
			elif self.__cycle[-2] in path: #end node
				endNodePath = path
			else:
				middlePath = path
		
		#error checking should not happen here, but let it be for the moment being....
		if startNodePath is None or endNodePath is None or middlePath is None:
			print "ASSERT ERROR In paths"
		
		if len(startNodePath)==0  or len(endNodePath)==0 or len(middlePath)==0:
			print "ASSERT ERROR In paths"
			
		 
		# cells[splitCell] = CELL_SPLITPOINT
		
		#
		self.__splitPaths.append([middlePath + [middlePath[0]], endNodePath + startNodePath + [endNodePath[0]]])
		
		# if (len(middlePath) == 0):
		# #perfect split when splitnodes where first and last element on given path
			# return cells, [startNodePath + [startNodePath[0]], endNodePath + [endNodePath[0]]]
		# else:
			# return cells,[middlePath + [middlePath[0]], endNodePath + startNodePath + [endNodePath[0]]]
	
def split_path_at_two_neighbours(path,splitA,splitB):
	#old  __split_path_loop
	# ASSERT:splitNode and checkNode are neighbours (direction unknown)
	
	#make sure path is not a cycle
	if path[0] == path[-1]:
		path = path[:-1]
	
	indexA = path.index(splitA)
	indexB = path.index(splitB)
	
	if indexA > indexB : #sort indeces, indexA must be smallest one
		indexA,indexB = indexB,indexA
	
	if (indexA == 0 and indexB == len(path)-1):
		#exception case where "neighbours are the first and last element on the path
		pathA = path
		pathB = []
	else:
		pathA = path[:indexA+1] #include splitNode
		pathB = path[indexB:] #could be empty...
	
	return [pathA,pathB]
		
	
		
def standardized_hamilton_cycle(path):
	#assert hamilton cycle with (row,col) tuples, first and last tuple equal
	#so it always starts with (0,0) and ends with (0,0)
	
	if path[0] ==(0,0):
		return path
	else:
		path = path[:-1]
		i = path.index((0,0))
		return path[i:]+ path[:i] +[(0,0)]
		
		
if __name__== "__main__":
	path = [ (2, 1), (2, 0), (3, 0),(3,1),(3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2),(2,1)]    #hamilton cycle
	# path = [ (2, 1), (2, 0), (3, 0),(3,1)],[(3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2)]    #valid paths
	
	ROWS = 4
	COLS = 6
	cycle = HamiltonCycle(ROWS,COLS,path)
	cycle.create_cell_pattern_from_hamilton_cycle()
	
	cycle.find_all_neighbour_hamilton_cycles()
	# print vars(cycle)["_HamiltonCycle__lattice"]
	
	checkPath = latticeOperations.Lattice(ROWS,COLS)
	for splitPath in cycle.get_split_paths():
		checkPath.new_paths_to_lattice(splitPath)
		checkPath.string_from_paths()
		
		print splitPath
		
		
	print cycle