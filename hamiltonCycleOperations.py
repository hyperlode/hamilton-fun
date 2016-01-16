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
	
	cycle = HamiltonCycle(4,6,path)
	cycle.create_cell_pattern_from_hamilton_cycle()
	# print vars(cycle)["_HamiltonCycle__lattice"]
	cycle.print_cells_ASCII()
	print cycle