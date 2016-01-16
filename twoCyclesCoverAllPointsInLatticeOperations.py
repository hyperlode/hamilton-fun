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
		self.__potentialRecombinationCells = []
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
		return approvedOutsiders
	
	
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