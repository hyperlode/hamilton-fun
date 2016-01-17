import latticeOperations

CELL_INSIDE = 0
CELL_OUTSIDE = 1 
CELL_SPLITPOINT_POTENTIAL = 2
CELL_SPLITPOINT = 3
CELL_RECOMBINATION_CANDIDATE = 4

class TwoCycles():
	def __init__(self, lattice_rows, lattice_cols, paths,cells,forbiddenRecombinator=None):
		#if forbiddenrecombinator is cell, assume that this was the split cell, and that the path ends lead to this cell.
		
		self.__lattice = latticeOperations.Lattice(lattice_rows, lattice_cols)
		self.__lattice.new_paths_to_lattice(paths)
		self.__paths = paths
		self.forbiddenRecombinationCell = forbiddenRecombinator
		self.__foundLoopConnectionCells = []
		self.__foundHamiltonCycles = []
		self.__cells_extra_info = cells
		self.__cells = {}
		
		
		for k,v in self.__cells_extra_info.items():
			if v == CELL_SPLITPOINT or v == CELL_RECOMBINATION_CANDIDATE or v == CELL_OUTSIDE:
				self.__cells[k] = CELL_OUTSIDE
			else:
				self.__cells[k] = v
		
	def print_cells_ASCII(self, withExtraInfo = True):
		
		if withExtraInfo:
			cells = self.__cells_extra_info
		else:
			cells = self.__cells
		
		printcells = ""
		for row in range(self.__lattice.rows() - 1):
			printrow = ""
			for col in range(self.__lattice.cols() - 1):
				if cells[(row,col)] ==CELL_INSIDE:
					printrow += "X"
				elif cells [(row,col)] == CELL_SPLITPOINT_POTENTIAL:
					printrow += "x"
				elif cells [(row,col)] == CELL_OUTSIDE:
					printrow += " "
				elif cells [(row,col)] == CELL_SPLITPOINT:
					printrow += "+"
				elif cells [(row,col)] == CELL_RECOMBINATION_CANDIDATE:
					printrow += "o"
				else:
					printrow += "?"
			printcells += printrow + "\n"
		print printcells
	
	#################################################
	#################################################
	#################################################
	
	def find_recombination_cells(self):
		#return all possible cells for recombination for two paths.
		
		#make sure shortest path comes first.
		shortestPath = self.__paths[0]
		if len(self.__paths[0])>len(self.__paths[1]):
			shortestPath =  self.__paths[1]
		
		#if the end of the path is included, the cell where the path was cut (cutting cell) will be discovered again as recombination candidate
		# end = -1
		# if self.forbiddenRecombinationCell is not None:
			# end= None
		end = None  #lets not play little games here...
		
		#find the adjecent cells to the path
		adjecentCells = []
		previousNode= shortestPath[0]
		for node in shortestPath[1:end]:
			#path is loop, so first node is repeated, so all node combinations (also the lats one) are found
			adjecentCells.extend(self.find_adjecentCells_from_two_path_coords(node,previousNode))
			previousNode = node
		#select only the "outside" cells
		outsideAdjecentCells = [cell for cell in adjecentCells if self.__cells[cell] == CELL_OUTSIDE]		
		
		
		#if a cell is twice more in the list, it means that it is bordering two cells on the path (inside of curve), this can never be a valid recombination cell, so it should be deleted!
		#this is an armpit situation
		allCellsOnce = set(outsideAdjecentCells) #first create list where all cells are exactly once
		for i in allCellsOnce:
			outsideAdjecentCells.remove(i) #remove the exactly once ones from the normal list, 
		
		for leftOversAreRepeaters in outsideAdjecentCells: #those that remain are the once that should be removed completly
			while leftOversAreRepeaters in allCellsOnce:
				#make sure all are removed.
				allCellsOnce.remove(leftOversAreRepeaters)
		
		outsideAdjecentCells = allCellsOnce
		
		#if an outside cell has three outside neighbours, it is not valid either (there is no border from the other path to connect to, only corners)
		approvedOutsiders = []
		for outsider in outsideAdjecentCells:
			outsideNeighboursCount = 0
			for neighbour in self.find_neighbour_cells(outsider):
				try:
					if self.__cells[neighbour] == CELL_OUTSIDE or self.__cells[neighbour] == CELL_RECOMBINATION_CANDIDATE:
						outsideNeighboursCount +=1
				except:
					print "the cells"
					
					raise
			if outsideNeighboursCount<3:
				approvedOutsiders.append(outsider)
			
		if self.forbiddenRecombinationCell is not None:
			try:
				approvedOutsiders.remove(self.forbiddenRecombinationCell)
			except:
				print approvedOutsiders
				print self.forbiddenRecombinationCell
				print Exception("the original splitpoint is not found in neighbours, or possible recombinators, this is shocking!")
		
		
		#write down in lattice
		for cell in approvedOutsiders:
			self.__cells_extra_info[cell]= CELL_RECOMBINATION_CANDIDATE
		
		self.__foundLoopConnectionCells  =approvedOutsiders
	
	
	def find_neighbour_cells(self, cell):
		neighbours = []
		#test horizontally
		if cell[1] < self.__lattice.cols() -2:
			#east
			neighbours.append((cell[0],cell[1]+1) )
		
		if cell[1] > 0:
			#West
			neighbours.append((cell[0],cell[1]-1) )
		
		if cell[0] < self.__lattice.rows()-2:
			#South
			neighbours.append((cell[0]+1,cell[1]) )
		
		if cell[0] > 0:
			#north
			neighbours.append((cell[0]-1,cell[1]) )	
			
		return neighbours
	
	
	def find_adjecentCells_from_two_path_coords(self,coord1, coord2):
		#ASSERT coord1 and coord2 are orthogonal neighbours on lattice
		#sort coords
		adjecent = []
		
		if coord1[0] == coord2[0]:
			#horizontal (row the same)
			#check sequence and make sure it is normalized (eliminate path direction)
			if coord1[1] >coord2[1]:
				coord1,coord2 = coord2,coord1
			try:
				self.__cells[coord1] 
				adjecent.append(coord1)
			except:
				#ASSERT cell outside  boundaries
				pass
				
			try:
				self.__cells[(coord1[0]-1,coord1[1])] 
				adjecent.append((coord1[0]-1,coord1[1]))
			except:
				#ASSERT cell outside  boundaries
				pass
			
		else:
			#vertical node on path
			if coord1[0] > coord2[0]:
				coord1,coord2 = coord2,coord1
			
			try:
				self.__cells[coord1] 
				adjecent.append(coord1)
			except:
				#ASSERT cell outside  boundaries
				pass
				
			try:
				self.__cells[(coord1[0],coord1[1]-1)] 
				adjecent.append((coord1[0],coord1[1]-1))
			except:
				#ASSERT cell outside  boundaries
				pass
			
			
		return adjecent
	
	def find_all_hamilton_neighbours(self):
		self.find_recombination_cells()
		for cell in self.__foundLoopConnectionCells:
			self.__foundHamiltonCycles.append(self.convert_cycles_to_hamilton_cycle(cell))
	
	def get_hamilton_cycles(self):
		return self.__foundHamiltonCycles
	
	
	def convert_cycles_to_hamilton_cycle(self,cell):
		#assert two closed paths in arg. paths
		#connects paths at given cell
		
		#check if node (corresponding to cell (top left node of cell that is) is in first path, if not , switch
		#assume cellId node is in pathA
		pathA = self.__paths[0]
		pathB = self.__paths[1]
		if cell not in pathA:
			pathB,pathA = pathA, pathB
		
		neighboursA = self.neighbours_of_node(pathA, cell)
		
		#check whether the reconnection cell will have horizontal or vertical neighbours.
		connectionCellHasVerticalNeighbourCells = True
		neighbourToReconnectCellNodeWith = (cell[0],cell[1]+1) #assumehorizontal connection
		neighbourOnSameSide = (cell[0]+1,cell[1])
		if neighbourToReconnectCellNodeWith in neighboursA: #if path is horizontal, there is a vertical reconnection
			neighbourToReconnectCellNodeWith, neighbourOnSameSide = neighbourOnSameSide, neighbourToReconnectCellNodeWith 
			connectionCellHasVerticalNeighbourCells = False
			
		pathA1, pathA2 = split_path_at_two_neighbours(pathA, cell, neighbourOnSameSide )
		pathA = pathA2 + pathA1
		
		pathB1, pathB2 = split_path_at_two_neighbours(pathB, (cell[0]+1, cell[1]+1), neighbourToReconnectCellNodeWith)
		pathB = pathB2 + pathB1

		if pathA[-1] not in self.__lattice.find_neighbour_nodes(pathB[0]):
		#check if neighbours on lattice.
			print pathA
			print pathB
			raise Exception("paths matched wrongly!")
		
		#generate new cells from path
		recombined = pathA + pathB + [pathA[0]]
		return recombined
		# cells = self.create_cell_pattern_from_hamilton_cycle(recombined)
		
			
	def neighbours_of_node(self, path, node):
		#from old neighbours_on_cycle in first program.
		
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
		
		
		