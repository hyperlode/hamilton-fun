import latticeOperations

CELL_INSIDE = 0
CELL_OUTSIDE = 1 
CELL_SPLITPOINT_POTENTIAL = 2
CELL_SPLITPOINT = 3
CELL_RECOMBINATION_CANDIDATE = 4

class TwoCycles():
	def __init__(self, lattice_rows, lattice_cols, paths,cells,forbiddenRecombinator=None):
		
		self.__lattice = latticeOperations.Lattice(lattice_rows, lattice_cols)
		self.__lattice.new_paths_to_lattice(paths)
		
		self.forbiddenRecombinationCell = forbiddenRecombinator
		self.__potentialRecombinationCells = []
		self.__cells_extra_info = cells
		self.__cells = self.__cells_extra_info.copy()
		
		if forbiddenRecombinator is None:
			for k,v in self.__cells:
				if v == CELL_SPLITPOINT:
					v=CELL_OUTSIDE
	
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
	
	def find_recombination_cells(self, paths, cells,includeOriginalCutCell = False):
		#return all possible cells for recombination for two paths.
		# possibleRecombinationCells= []
		
		#paths is two cycles that, together all nodes in the lattice are taken.
		
		
		#make sure shortest path comes first.
		shortestPath = paths[0]
		if len(paths[0])>len(paths[1]):
			shortestPath =  paths[1]
		
		#if the end of the path is included, the cell where the path was cut (cutting cell) will be discovered again as recombination candidate
		end = -1
		if includeOriginalCutCell:
			end= None
		
		#find the adjecent cells to the path
		adjecentCells = []
		previousNode= shortestPath[0]
		for node in shortestPath[1:end]:
			adjecentCells.extend(self.find_adjecentCells_from_two_path_coords(node,previousNode,cells))
			previousNode = node
			
		#select only the "outside" cells
		#ERROR this is not enough, armpit situations should be filtered out!!! (when the cells make a 90 degree turn, the adject cell is 
		outsideAdjecentCells = [cell for cell in adjecentCells if cells[cell] == CELL_OUTSIDE]		
		print"jdijei"
		print outsideAdjecentCells
		
		#if a cell is twice more in the list, it means that it is bordering two cells on the path (inside of curve), this can never be a valid recombination cell, so it should be deleted!
		allCellsOnce = set(outsideAdjecentCells)
		for i in allCellsOnce:
			outsideAdjecentCells.remove(i)
		
		print allCellsOnce
		print outsideAdjecentCells
		for leftOversAreRepeaters in outsideAdjecentCells:
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
					if cells[neighbour] == CELL_OUTSIDE or cells[neighbour] == CELL_RECOMBINATION_CANDIDATE:
						outsideNeighboursCount +=1
				except:
					print "the cells"
					print cells
					raise
			if outsideNeighboursCount<3:
				approvedOutsiders.append(outsider)
		
		# print "outside adjectenss:"
		# print allCellsOnce
		
		#write down in lattice
		for cell in allCellsOnce:
			cells[cell]= CELL_RECOMBINATION_CANDIDATE
		return approvedOutsiders,cells
	