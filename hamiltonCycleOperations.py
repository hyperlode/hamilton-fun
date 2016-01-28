import latticeOperations
import twoCyclesCoverAllPointsInLatticeOperations
import random 
import graphs
import fileOperations

CELL_INSIDE = 0
CELL_OUTSIDE = 1 
CELL_SPLITPOINT_POTENTIAL = 2
CELL_SPLITPOINT = 3
CELL_RECOMBINATION_CANDIDATE = 4

class HamiltonCycle():
	
	def __init__(self, lattice_rows, lattice_cols, cycle,cycleByName = None):
		self.__lattice = latticeOperations.Lattice(lattice_rows, lattice_cols)
		
		if cycleByName is not None:
			#generate cycle from given name
			
			cycle = self.initFromName(cycleByName)
			
			
		#if no cycle provided, make one warning: super slow for bigger lattices, better just provide one!!!
		if cycle is None:
			# cycle = self.__lattice.find_hamilton_cycle()
			# print cycle
			cycle = self.do_manual_infill()
		
		#let path start from (0,0)  here we don't yet really check for hamilton, but if error, we know it is wrong already
		try:
			cycle = standardized_hamilton_cycle(cycle)
		except:
			print "make sure to provide a single hamilton cycle"
			raise 
			
		self.addCycle(cycle) 
		#from here cycle is asserted to be hamilton
		self.__cells_extra_info = {}
		self.__cells = {}
		
		self.__cycle = cycle 
		
		self.__splitCells = [] #all the potential splitting points.
		self.__splitPathsData = [] 
		self.__neighbourHamiltonCycles = []
		
		self.create_cell_pattern_from_hamilton_cycle()
	
	def initFromName(self, cycleByName):
		
		self.__cells = self.cells_from_nameString(self.__lattice.rows(),self.__lattice.cols(),cycleByName)
		
		return self.__cycle_from_cells()
		
	
	def __str__(self):
		self.__lattice.string_from_paths()
		return str(self.__lattice)
	
	
	def is_hamilton_cycle_existing(self):
		#only cycle possible if even number of nodes in lattice.
		return self.__lattice.rows() * self.__lattice.cols() % 2 == 0
	
	def do_manual_infill(self):
		if not self.is_hamilton_cycle_existing():
			raise Exception("no hamilton cycle possible")
		
		path = self.__infillRecursive([(0,0)])
		path = path + [path[0]]
		return path
	def get_cycle_as_tuple(self):
		return tuple(self.__cycle)
		
	def __infillRecursive(self,path):
		
		n = path[-1]
		if self.__lattice.cols()%2 != 0:
			directions = [(n[0]+1,n[1]),(n[0],n[1]+1),(n[0],n[1]-1),(n[0]-1,n[1])]
		else:
			directions = [(n[0],n[1]+1),(n[0]+1,n[1]),(n[0]-1,n[1]),(n[0],n[1]-1)]
			
		orthoLatticeNodeNeighbours= self.__lattice.find_neighbour_nodes(n)
		if directions[0] in orthoLatticeNodeNeighbours and directions[0] not in path:
			path.append(directions[0])
			self.__infillRecursive(path)
		elif directions[1] in orthoLatticeNodeNeighbours and directions[1] not in path:
			path.append(directions[1])
			self.__infillRecursive(path)
		elif directions[2] in orthoLatticeNodeNeighbours and directions[2] not in path:
			path.append(directions[2])
			self.__infillRecursive(path)
		elif directions[3] in orthoLatticeNodeNeighbours and directions[3] not in path:
			path.append(directions[3])
			self.__infillRecursive(path)
		# elif len(path) == self.__rows * self.__cols:
			# return
		else:
			pass
			
		return path	
			
			
		
		
	def addCycle(self,path):
		#add hamilton cycle, and check if it is one
		self.__lattice.new_paths_to_lattice(path)
		stats = self.__lattice.stats()
		if not stats["hamilton"] or stats["cycles"]!=1:
			print self
			raise Exception("no hamilton path provided")
	
			
	def get_split_pathsData(self):
		return self.__splitPathsData
		
	def create_cell_pattern_from_hamilton_cycle(self):
		#if a cycle is drawn, there is an inside and an outside. 
		#each four points of the lattice define a cell. This cell is in or outside the cycle.
		
		#default all cells are outside path
		self.__cells = {(r,c):CELL_OUTSIDE for r in range(self.__lattice.rows()-1) for c in range(self.__lattice.cols()-1)}
		try:
			self.__explore_inside_of_hamilton_cycle((0,0)) 
		except:
			print "creation of cells failed"
			raise
		self.__cells_extra_info = self.__cells.copy()
	
	def __explore_inside_of_hamilton_cycle (self, nodeCell):
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
				self.__explore_inside_of_hamilton_cycle( nextNodeCell )

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
	
	def cells_to_string(self, rowDivider = None, withExtraInfo  = False, asListOfStrings = False):
		#rowDivider is a string one or more character to put between rows, or value is ROWNUMBER, it will then output the rowNumber as divider
		
		newLineChar  = rowDivider
		if withExtraInfo:
			cells = self.__cells_extra_info
		else:
			cells = self.__cells
		
		if asListOfStrings:
			printcells = []
		else:
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
					
			if asListOfStrings:
				printcells.append(printrow)
			else:
				if rowDivider == "ROWNUMBER":
					# newLineChar = str(row + 1)
					newLineChar = '{number:0{width}d}'.format(width=len(str(self.__lattice.rows() - 2)), number=row+1)
				if newLineChar is not None:
					printcells += printrow + newLineChar
				else:
					printcells += printrow
		
		if not asListOfStrings and newLineChar is not None and len(newLineChar)>0:
			#no divider on end of string.
			printcells = printcells[:-len(newLineChar)]
			
		return printcells
		
		
	def print_cells_ASCII(self, withExtraInfo = True):
		
		print self.cells_to_string("\n", withExtraInfo)
	
	# def get_cycle_as_nameString(self, rowDivider = "_",withExtraInfo = False):
	def get_cycle_as_nameString(self, rowDivider = "", withExtraInfo = False):
		# rowDivider ="ROWNUMBER"   #my little baby, but unnecessary overkill!
		return self.cells_to_string(rowDivider = rowDivider, withExtraInfo = withExtraInfo)
	
	def get_isoMorphs_as_nameString(self, noFlippedOnes = False, includeOriginalName = False):
		baseNameAsList = self.cells_to_string( rowDivider = None, withExtraInfo  = False, asListOfStrings = True)
		# baseNameAsList = ["123","456"]
		
		baseNameString = "".join(baseNameAsList)
		
		
		
		# print baseNameString
		
		#180deg
		rot180 = baseNameString[::-1]
		# print rot180
		# isoMorphs.append(rot180)
		#flipHori
		flipHori = "".join([ row[::-1] for row in baseNameAsList ])
		# print flipHori
		# isoMorphs.append(flipHori)
		#flipVert
		flipVert = "".join(baseNameAsList[::-1])
		# print flipVert
		# isoMorphs.append(flipVert)
		#90deg
		rot90 = []
		# print self.__lattice.cols()
		for col in range(self.__lattice.cols() -1):
		# for col in range(3):
			# print "lode"
			transposedRow = ""
			# print self.__lattice.rows()
			for row in range(self.__lattice.rows() - 1):
			# for row in range(2):
				# print "brecht"
				transposedRow += baseNameAsList[row][col]
			rot90.append(transposedRow)
		
		rot90String = "".join(rot90)
		# print rot90String
		# isoMorphs.append(rot90String)
		#rot 90 hori flip
		rot90HoriFlip = "".join([row[::-1] for row in rot90])
		# print rot90HoriFlip
		# isoMorphs.append(rot90HoriFlip)
		#rot270
		rot270String = rot90String[::-1]
		# print rot270String
		# isoMorphs.append(rot270String)
		#rot 90 flip vert
		rot90FlipVert = "".join(rot90[::-1])
		# print rot90FlipVert 
		# isoMorphs.append(rot90FlipVert)
		
		isoMorphs = [] 
		
		#add original path
		if includeOriginalName:
			isoMorphs.append(baseNameString)
		
		#add rotations
		isoMorphs.extend([rot90HoriFlip,rot180,rot90FlipVert])
		
		#add flipped rotations
		if not noFlippedOnes:
			isoMorphs.extend( [ flipHori,rot270String,flipVert,rot90String])
		
		
		return isoMorphs
		
	def cells_from_nameString(self,rows,cols,nameString):
		#dangerous! check result!
		#create all outside cells for this situation
		#default all cells are outside path
		cells = {(r, c):CELL_OUTSIDE for r in range(rows - 1) for c in range(cols - 1)}
		
	
		#autoRowDividerCheck
		extraLength = len(nameString) - (rows-1) * (cols-1)
		newLineLength = extraLength / (rows - 2) #check length of one newline spacer
		
		# nameString = nameString[cols-1:cols + newLineLength+1]
		cleanedNameString = ""
		
		for r in range(rows-1): 
			# print nameString [  r * (cols + newLineLength  -1) : r * (cols + newLineLength -1) + cols - 1 ]
			# cleanedNameString += nameString [  r * (cols + newLineLength + 1 ) : r * (cols + newLineLength + 1) + cols - 1 ]
			cleanedNameString += nameString [  r * (cols + newLineLength -1 ) : r * (cols + newLineLength - 1) + cols - 1 ]
		
		index = 0
		for row in range(rows-1):
			for col in range(cols-1):	
				
				if cleanedNameString[index] == "X":
					cells[(row,col)] = CELL_INSIDE
				elif  cleanedNameString[index] == " ":
					cells[(row,col)] = CELL_OUTSIDE
				else:
					print cleanedNameString
					raise Exception("strange values in string, only inside and outside allowed")
				
				
				index += 1
		return cells
	
	def __cycle_from_cells(self):
		#only used for alternative init with name
		#dangerous! check result!
		#def getPostitBoundary(self):
		'''
		numbers= postit corners (defined)
		
		
		'''
		'''
		get circumference of postit notes (if cannot be closed as one path (= start point = endpoint), illegal!
		
		walk from posititcorner to postit corner
		'''
		
		#first elements are always the same
		boundary = [(0,0), (0,1)]
		
		lengthBoundary  = (self.__lattice.rows() ) * (self.__lattice.cols() ) #defined, lenght of the path is equal to "cornerpositions" as each cornerposition is always taken.
		# print "jijijijijiji"
		nodeBackup = None
		# print self.__cells
		# self.print)
		index = 0
		while len(boundary) < lengthBoundary:
			# index += 1
			# if index > 50:
				# raise
			#print len(boundary)
			#get coord of active cell 
			# row,col =  
			activeNode = boundary[-1] 
			
			#get ortho neighbour coords:
			neighbourNodes = self.__lattice.find_neighbour_nodes(activeNode)
			
			#freeNode 
			nodesToCheck = [node for node in neighbourNodes if node not in boundary ]
			addedNode = False
			
			for node in nodesToCheck:
				# print "========"
				# print node
				adjecentCells = self.find_adjecentCells_from_two_path_coords(activeNode, node)
				# print adjecentCells
				
				if len(adjecentCells)== 2:
					if self.__cells[adjecentCells[0]] != self.__cells[adjecentCells[1]]:
						#different (one outside, one inside cell) cell types
						addedNode = True
						boundary.append(node)
						# print "this is the node"
						break
				elif len(adjecentCells)== 1:
					nodeBackup = node
					# boundary.append(node)
					# print "this is the nodeeee"
					# break
			if not addedNode:
				# print nodeBackup
				# print "added as backup"
				boundary.append(nodeBackup)
				
		return boundary	+ [boundary[0]]
			
			
		
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
	
	######################################
	#neighbour finder
	######################################
	
	def find_all_neighbour_hamilton_cycles(self):
		self.find_split_cells()
		
		#find cells where cycle can be split
		for cell in self.__splitCells:
			self.split(cell)
		
		#for every splitcell there are two paths, for every two paths get all the recombinations. 
		for splitPath in self.get_split_pathsData():
			twoLoops = twoCyclesCoverAllPointsInLatticeOperations.TwoCycles(self.__lattice.rows(), self.__lattice.cols(),splitPath["paths"],splitPath["cells"],splitPath["splitCell"])
			twoLoops.find_all_hamilton_neighbours() #search the neighbours
			self.__neighbourHamiltonCycles.extend(twoLoops.get_hamilton_cycles()) #get the neighbours
		
		#generalize the found cycles (all paths should start with 0,0)
		self.__neighbourHamiltonCycles = [standardized_hamilton_cycle(cycle) for cycle in self.__neighbourHamiltonCycles]
		
			
	def get_hamilton_cycle_neighbours(self):
		return self.__neighbourHamiltonCycles
			
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
			self.__cells_extra_info[cell] = CELL_SPLITPOINT_POTENTIAL 
		
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
			
		
		cellsWithSplitpoint = self.__cells.copy()
		cellsWithSplitpoint[splitCell] = CELL_SPLITPOINT

		self.__splitPathsData.append({"cells": cellsWithSplitpoint, "splitCell": splitCell, "paths":[middlePath + [middlePath[0]], endNodePath + startNodePath + [endNodePath[0]]]})
	
	def get_cycle_as_detailed_nameString(self):
		'''
		name that has meaning for sorting and classification. According to cells.
		
		3 elements:
		-every ring has a letter: outer ring = A , then B, ...
		-every side has an orientation: upper side = North, right side = East,...
		-every cell per ring per side has a number. starting from zero,... 
		
		Corners are part of the "most clockwise" orientation.
		
		if rows and cols are equal, every cycle has four 90deg rotations, all mirrored = 8 cycles for one
		if rows and cols not equal: 4 cycles for one.
		
		
		
		'''
		
		cellRows = self.__lattice.rows() - 1
		cellCols = self.__lattice.cols() - 1
		
		#check amount of rings:
		#select highest number, if it is odd, add one, divide it by two.
		keep,b = self.__lattice.rows(), self.__lattice.cols()
		if keep<b:
			keep,b = b,keep
		
		rings = keep/2
		
		# #if a single cell is in the center, it will not be detected in the normal program flow, we will define it a North
		# if cellRows == cellCols:
			# centerCellAvailable = True
		# else:
			# centerCellAvailable = False
			
		if rings>26:
			print "ASSERT ERROR: max 26 rings implemented..."
			raise
		
		nameString = "HamiltonCycle_r{}c{}_openCells:".format(self.__lattice.rows(), self.__lattice.cols())
		
		#name
		for ringI, ring in enumerate(["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","X","Y","Z"][:rings]):
			# print "ring:{}".format(ring)
			
			#NORTH
			for i in range(cellCols - 1 - 2*ringI):
				if self.__cells[(ringI,i + ringI)] == CELL_OUTSIDE:
					nameString += ring+"1N"+ str(i) + "-"
				# print "{},{}".format(ringI,i + ringI)
			
			#EAST
			for i in range(cellRows - 1 - 2*ringI):
				if self.__cells[(i + ringI, cellCols - 1 - ringI)] == CELL_OUTSIDE:
					nameString += ring+"2E"+ str(i) + "-"
				# print "{},{}".format(i + ringI, cellCols - 1 - ringI)
			
			#SOUTH
			for i in range(cellCols - 1 - 2*ringI):
				if self.__cells[(cellRows - 1 - ringI, cellCols - 1 -  i - ringI)] == CELL_OUTSIDE:
					nameString += ring+"3S"+ str(i) + "-"
				# print "{},{}".format(cellRows - 1 - ringI, cellCols - 1 -  i - ringI)
				
			#WEST
			for i in range(cellRows - 1 - 2*ringI):
				if self.__cells[( cellRows - 1 - i - ringI, ringI)] == CELL_OUTSIDE:
					nameString += ring+"4W"+ str(i) + "-"
				# print "{},{}".format( cellRows - 1 - i - ringI, ringI)
			
		#the center cell in rows==cols lattices is an exception, and defined as North.
		if cellRows == cellCols:
			if self.__cells[(ringI, ringI)] == CELL_OUTSIDE:
				nameString += ring+"1N0" + "-"
			# print "{},{}".format(  ringI, ringI)
			
			
		return nameString[:-1]
		# # # #name
		# # # for ringI, ring in enumerate(["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","X","Y","Z"][:rings]):
			# # # print "ring:{}".format(ring)
			
			# # # #NORTH
			# # # for i in range(cellCols - 1 - 2*ringI):
				# # # print ring+"N"+ str(i)
				# # # print "{},{}".format(ringI,i + ringI)
			
			# # # #EAST
			# # # for i in range(cellRows - 1 - 2*ringI):
				# # # print ring+"E"+ str(i)
				# # # print "{},{}".format(i + ringI, cellCols - 1 - ringI)
			
			
			# # # #SOUTH
			# # # for i in range(cellCols - 1 - 2*ringI):
				# # # print ring+"S"+ str(i)
				# # # print "{},{}".format(cellRows - 1 - ringI, cellCols - 1 -  i - ringI)
			
			# # # #WEST
			# # # for i in range(cellRows - 1 - 2*ringI):
				# # # print ring+"W"+ str(i)
				# # # print "{},{}".format( cellRows - 1 - i - ringI, ringI)
			
		# # # #the center cell in rows==cols lattices is an exception, and defined as North.
		# # # if cellRows == cellCols:
			# # # print ring+"N"+ str(0)
			# # # print "{},{}".format(  ringI, ringI)
			

		
		#go over each cell

		
		#if empty cell, define it and add to string.
		
		#return string
		
		pass
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
		
def getNeighbourCycles(rows, cols, path):
	cycle = HamiltonCycle(ROWS,COLS,path)
	# cycle.print_cells_ASCII(True)
	nstr = cycle.get_cycle_as_nameString("")
	
	cycle.find_all_neighbour_hamilton_cycles()
	neighbourCycles = cycle.get_hamilton_cycle_neighbours()
	# print nstr
	# retest = HamiltonCycle(ROWS,COLS,None,nstr)
	# print "=======!!!!!===="
	# retest.print_cells_ASCII(True)
	
	# print retest.get_cycle_as_nameString()
	# print "==========="
	return neighbourCycles

def getNeighbourCyclesAsNames(rows, cols, cycleName, includeIsoMorphs = False):
	#if includeIsoMorphs==True will also include all rotations and symmetricals from every neighbour
	
	cycle = HamiltonCycle(rows,cols,None,cycleName)
	
	cycle.find_all_neighbour_hamilton_cycles()
	
	neighbourCycleNames = []
	neighbourCycleNamesIsoMorphs_set = set([])
	
	neighbourCycles = cycle.get_hamilton_cycle_neighbours()
	for cycle in neighbourCycles:
		nameCycle = HamiltonCycle(ROWS,COLS,cycle)
		neighbourCycleNames.append(nameCycle.get_cycle_as_nameString(""))
		
		if includeIsoMorphs:
			neighbourCycleNamesIsoMorphs_set |= set(nameCycle.get_isoMorphs_as_nameString(noFlippedOnes = False, includeOriginalName = False))
		
	if includeIsoMorphs:
		return neighbourCycleNames, neighbourCycleNamesIsoMorphs_set
	else:
		return neighbourCycleNames

def getNeighbourCycles(rows, cols, cycleData):
	cycle = HamiltonCycle(rows,cols,cycleData)
	cycle.find_all_neighbour_hamilton_cycles()
	neighbourCyclesData = cycle.get_hamilton_cycle_neighbours()
	
	return neighbourCyclesData
	
def getAllPossibilities(rows, cols):
	initCycle = HamiltonCycle(rows,cols,None)
	initCycleName = initCycle.get_cycle_as_nameString("")
	cyclesToInvestigate = set([initCycleName])
	all = set([])
	i = 0
	while len(cyclesToInvestigate)>0:
		
		checkCycleName = cyclesToInvestigate.pop() #get the cycle to investigate
		all.add(checkCycleName) #add cycle to found ones
		# checkCycle = HamiltonCycle(rows,cols,None,checkCycleName)
		
		checkCycleNeighbourNames = getNeighbourCyclesAsNames(rows, cols, checkCycleName)
		
		for cycle in checkCycleNeighbourNames:
			if cycle not in all:
				cyclesToInvestigate.add(cycle)
		if i % 200 == 0:
			print "-------stats----"
			print len(cyclesToInvestigate)
			print len(all)
		if i% 500000 == 0 and i > 400000:
			fileOperations.linesToFile( r"c:\temp\foundHamiltonCyclesFor{}rows_{}cols_step{}.txt".format(ROWS,COLS,i),list(all)+ list(cyclesToInvestigate))
		i += 1
	return all
	
	
def getAllPossibilities_fast(rows, cols):
	initCycle = HamiltonCycle(rows,cols,None)
	initCycleData = initCycle.get_cycle_as_tuple()
	
	# fileOperations.dumpDataPickle(tmp, r"c:\temp\tmpdata.txt")
	# print fileOperations.retrieveDataPickle(r"c:\temp\tmpdata.txt")
	# initCycleName = initCycle.get_cycle_as_nameString("")
	cyclesToInvestigate = set( [])
	cyclesToInvestigate.add(initCycleData)
	
	all = set([])
	i = 0
	while len(cyclesToInvestigate)>0:
		
		checkCycleData = cyclesToInvestigate.pop() #get the cycle to investigate
		all.add(checkCycleData) #add cycle to found ones
		checkCycleNeighboursData = getNeighbourCycles(rows, cols, list(checkCycleData))
		
		for cycle in checkCycleNeighboursData:
			cycle_as_tuple = tuple(cycle)
			if cycle_as_tuple not in all:
				cyclesToInvestigate.add(cycle_as_tuple)
		if i % 200 == 0:
			print "-------stats----"
			print len(cyclesToInvestigate)
			print len(all)
		# if i% 5 == 0 and i > 4:
		if i% 500000 == 0 and i > 400000:
			fileOperations.dumpDataPickle(all, r"c:\temp\foundHamiltonCyclesFor{}rows_{}cols_neighboursFound_step{}.txt".format(ROWS,COLS,i))
			fileOperations.dumpDataPickle(cyclesToInvestigate, r"c:\temp\foundHamiltonCyclesFor{}rows_{}cols_neighboursToInvestigate_step{}.txt".format(ROWS,COLS,i))
		i += 1
	
	return all

def getAllPossibilities_super_fast(rows, cols):
	#from every neighbour, also takes isomorphs (rotations and symmetries) right away.
	
	initCycle = HamiltonCycle(rows,cols,None)
	initCycleName = initCycle.get_cycle_as_nameString("")
	cyclesToInvestigate = set([initCycleName])
	all = set([])
	i = 0
	while len(cyclesToInvestigate)>0:
		
		checkCycleName = cyclesToInvestigate.pop() #get the cycle to investigate
		 #add cycle to found ones
		
		checkCycleNeighbourNames, isoMorphs_set = getNeighbourCyclesAsNames(rows, cols, checkCycleName, True)
		# checkCycleNeighbourNames = getNeighbourCyclesAsNames(rows, cols, checkCycleName, False)
		
		
		
		
		
		for cycle in checkCycleNeighbourNames:
			if cycle not in all:
				cyclesToInvestigate.add(cycle)
		
		all |= isoMorphs_set
		all.add(checkCycleName)

		
		if i % 200 == 0:
			print "-------stats----"
			print len(cyclesToInvestigate)
			print len(all)
			# print checkCycleNeighbourNames
			
		if i% 500000 == 0 and i > 400000:
			fileOperations.linesToFile( r"c:\temp\foundHamiltonCyclesFor{}rows_{}cols_step{}.txt".format(ROWS,COLS,i),list(all)+ list(cyclesToInvestigate))
		i += 1
	return all
	
def convertCyclesDataToNames(rows, cols, cyclesData):
	cyclesNames = []
	for cycleData in cyclesData:
		cycle = HamiltonCycle(rows,cols,list(cycleData))
		cyclesNames.append(cycle.get_cycle_as_nameString(""))
	return cyclesNames

def convertCyclesDataToDetailedNames(rows, cols, cyclesData):
	detailedNames = []
	for cycleData in cyclesData:
		cycle = HamiltonCycle(rows,cols,list(cycleData))
		detailedNames.append(cycle.get_cycle_as_detailed_nameString())
	return detailedNames

def convertCyclesNamesToDetailedNames(rows, cols, cyclesNames):
	detailedNames = []
	for name in cyclesNames:
		cycle = HamiltonCycle(rows,cols,None, name)
		detailedNames.append(cycle.get_cycle_as_detailed_nameString())
	return detailedNames


def printIsoMorphs(rows, cols, withFlippedOnes = True, cycleData = None, addRandomness = True):
	if cycleData is None:
		print "creating own cycledata with some iterations..."
		startCycle = None
		for i in range(100):
			initCycle = HamiltonCycle(ROWS,COLS,startCycle)
			initCycle.find_all_neighbour_hamilton_cycles()
			if addRandomness:
				startCycle = random.choice(initCycle.get_hamilton_cycle_neighbours())
			else:
				startCycle = initCycle.get_hamilton_cycle_neighbours()[1]
		
	else:
		initCycle = HamiltonCycle(ROWS,COLS,cycleData)
	print "original cycle:"
	initCycle.print_cells_ASCII(False)	
	print "all isomorphs:"
	isoMorphs =  initCycle.get_isoMorphs_as_nameString()
		
	
	for isoMorph in isoMorphs:
		initCycle = HamiltonCycle(ROWS,COLS,None,isoMorph)
		initCycle.print_cells_ASCII(False)	
		print ""
	
if __name__== "__main__":
	path = [ (2, 1), (2, 0), (3, 0),(3,1),(3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2),(2,1)]    #hamilton cycle
	# path = [ (2, 1), (2, 0), (3, 0),(3,1)],[(3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2)]    #valid paths
	path = None
	ROWS = 8
	COLS = 8
	print "ROWS:{}, COLS:{}".format(str(ROWS),str(COLS))
	
	
	
	
	timePointAnchor = fileOperations.getTime()
	 
	# allCyclesData = getAllPossibilities_fast(ROWS,COLS)
	# allCycleNames_A = convertCyclesDataToDetailedNames(ROWS, COLS, allCyclesData)
	# print len(allCycleNames_A)
	# fileOperations.linesToFile( r"c:\temp\foundHamiltonCyclesFor{}rows_{}cols_detailedNameString.txt".format(ROWS,COLS),list(allCycleNames_A))
	# deltaT = fileOperations.getTime() - timePointAnchor
	# timePointAnchor = fileOperations.getTime() - timePointAnchor
	# print deltaT
	
	# allCyclesData = getAllPossibilities(ROWS,COLS)
	# allCycleNames_A = convertCyclesDataToDetailedNames(ROWS, COLS, allCyclesData)
	# print len(allCycleNames_A)
	# deltaT = fileOperations.getTime() - timePointAnchor
	# timePointAnchor = fileOperations.getTime() - timePointAnchor
	# print deltaT
	
	allCycleNames = getAllPossibilities_super_fast(ROWS,COLS)
	allCycleNames_A = convertCyclesNamesToDetailedNames(ROWS, COLS, allCycleNames)
	fileOperations.linesToFile( r"c:\temp\foundHamiltonCyclesFor{}rows_{}cols_detailedNameString.txt".format(ROWS,COLS),allCycleNames_A)
	print len(allCycleNames)
	deltaT = fileOperations.getTime() - timePointAnchor
	timePointAnchor = fileOperations.getTime() - timePointAnchor
	print deltaT
	
	# fileOperations.linesToFile( r"c:\temp\foundHamiltonCyclesFor{}rows_{}cols_detailedNameString.txt".format(ROWS,COLS),list(allCycleNames_A))
	# fileOperations.linesToFile( r"c:\temp\foundHamiltonCyclesFor{}rows_{}cols_BBBB.txt".format(ROWS,COLS),list(allCycleNames))
	print "done"
	
	# initCycle = HamiltonCycle(ROWS,COLS,None)
	# initCycle.print_cells_ASCII(False)
	# neighbours = initCycle.getNeighbourCycles()
	
	
	# printIsoMorphs(ROWS, COLS, True, None, True)