import graphs
import fileOperations
import random
from collections import defaultdict

CELL_INSIDE = 0
CELL_OUTSIDE = 1 
CELL_SPLITPOINT_POTENTIAL = 2
CELL_SPLITPOINT = 3
CELL_RECOMBINATION_CANDIDATE = 4

class HamiltonFun:
	def __init__(self, lattice_rows, lattice_cols):
		""" initializes a graph object """
		self.__rows = lattice_rows
		self.__cols = lattice_cols
		self.__rowCells = self.__rows-1
		self.__colCells = self.__cols-1
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
		
	def is_pattern_hamilton_cycle(self, cells):
		#number of inside cells must be 
		#((n-2)**2/2) bramz formula: empty spaces (n = points!!!! no postits (points=  corners of postits i.e.= 5x5 postits = 6x6 points
		#self.emptySpaces = (rows + 1-2) *  (cols + 1-2) / 2
		
		emptySpaces = (self.__rowCells + 1-2) *  (self.__colCells + 1-2) / 2
		
		#chapter COUNTING in https://codefisher.org/catch/blog/2015/04/22/python-how-group-and-count-dictionaries/ 
		d = defaultdict(int)
		for value in cells.values():
			d[value] += 1	
		
		
		# 
		# print "0-0-0--"
		if (emptySpaces != d[CELL_OUTSIDE] +  d[CELL_RECOMBINATION_CANDIDATE]):
			print "is hamiltoncycle?:"
			print self.print_cells_ASCII(cells)
			print emptySpaces == d[CELL_OUTSIDE] +  d[CELL_RECOMBINATION_CANDIDATE]
			print d
			print d[CELL_OUTSIDE]
			print d[CELL_RECOMBINATION_CANDIDATE]
			print cells
			
			raise "ASSERT ERROR empty spaces in hamilton not OK"
		return emptySpaces == d[CELL_OUTSIDE] +  d[CELL_RECOMBINATION_CANDIDATE]
		
	def create_cell_pattern_from_hamilton_cycle(self, path):
		#if a cycle is drawn, there is an inside and an outside. 
		#each four points of the lattice define a cell. This cell is in or outside the cycle.
		
		#ASSERT hamilton cycle
		
		#default all cells are outside path
		cells = {(r,c):CELL_OUTSIDE for r in range(self.__rows-1) for c in range(self.__cols-1)}
		try:
			return self.__explore_inside_of__hamilton_cycle(path, cells, (0,0)) 
		except:
			print 'foiajef;oija;oiefj;iajse;fiajsef'
			print path
			raise
	
	
	
	def __explore_inside_of__hamilton_cycle (self, path, cells, nodeCell):
		#every cell has four directions, N, E, S, W . inside cells are all connected. Go over inside recursively
		#add nodeCell as inside in cells
		cells[nodeCell] = CELL_INSIDE
		
		#four directions, (common nodes with adjecent cell and adjecent cell )
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
			pathNeighboursNodeA = self.neighbours_on_cycle(path, nodeA)
			# pathNeighboursNodeB = self.neighbours_on_cycle(path, nodeB)
			
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
				neighboursOnPathTopLeftNode =  self.neighbours_on_cycle(path,(row,col))
				
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
				neighboursOnPathBottomRightNode =  self.neighbours_on_cycle(path,(row+1,col+1))
				
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
			cells[cell] = CELL_SPLITPOINT_POTENTIAL 
		
		return cells, splitCells
	
	def split(self, path, splitCell, cells):
		#split the path in two at splitcell
		nodeTopLeftPositionInPath = path.index(splitCell)
		nodeBottomRightPositionInPath = path.index((splitCell[0]+1, splitCell[1]+1))
		
		
		#  drA, dcA,drB,dcB,drC,dcC in [(0,0,0,1,-1,0),(0,1,1,1,0,1),(1,1,1,0,1,0),(1,0,0,0,0,-1)]:
		#determine if vertical or horizontal
		neighbours = self.neighbours_on_cycle(path, splitCell)
		#orthogonal neighbours of first node E and S
		#orthoNodesInLattice =  {(row,col+1),(row+1,col)} #create as set
		
		# print "neighs:"
		# print neighbours
		#preserve start and end node of path
		startNode = path[0]
		endNode = path[-2] #-2 because of loop (repetition)
		
		#save state of split cell
		isVertical = True
		splitNodeOpposite = (splitCell[0],splitCell[1]+1)
		splitNodeNeighbourOnPathOfSplitcellNode = (splitCell[0]+1,splitCell[1])
		
		if  (splitNodeOpposite) in neighbours:
			#the connection is horizonal (col +1)
			splitNodeNeighbourOnPathOfSplitcellNode, splitNodeOpposite = splitNodeOpposite, splitNodeNeighbourOnPathOfSplitcellNode
			isVertical = False
		# print "isHorizontally connected?  =   :{}".format(str(not isVertical))
		# print splitNodeOpposite
		# print splitNodeNeighbourOnPathOfSplitcellNode
		# print "kdkdk"
		# print "papapapa"
		# print path
		
		twoPaths = self.__split_path_loop(path[:-1],splitCell,splitNodeNeighbourOnPathOfSplitcellNode)
		#if the two nodes are front and end, one of the returning paths will be empty.
		pathA = twoPaths[0]
		pathB = twoPaths[1]
		
		# print "fijasiejf"
		# print twoPaths
		
		diagonalSplitNode = (splitCell[0]+1,splitCell[1]+1)
		
		if (diagonalSplitNode) not in pathA:
			pathA,pathB = pathB,pathA
		
		# print "pathA"
		# print pathA
		
		try:
			twoPaths = self.__split_path_loop(pathA,diagonalSplitNode,splitNodeOpposite)
		except:
			print pathA
			print pathB
			
			print self.print_paths_ASCII([pathA,pathB])
			print "ifjieiejfiief"
			raise 'ffefe'
			
		pathA = twoPaths[0]
		pathC = twoPaths[1]
		
		paths = [pathA,pathB,pathC]
		
		startNodePath = None
		endNodePath = None
		middlePath = None
		for path in paths:
			if startNode in path:
				startNodePath = path
			elif endNode in path:
				endNodePath = path
			else:
				middlePath = path
				
		if startNodePath is None or endNodePath is None or middlePath is None:
			print "ASSERT ERROR In paths"
		cells[splitCell] = CELL_SPLITPOINT
		
		if (len(middlePath) == 0):
		#perfect split when splitnodes where first and last element on given path
			return cells, [startNodePath + [startNodePath[0]], endNodePath + [endNodePath[0]]]
		else:
			return cells,[middlePath + [middlePath[0]], endNodePath + startNodePath + [endNodePath[0]]]
		
		#concatenate the two pieces that are one path
		
		
		# twoPaths = [twoPaths[0] + [twoPaths[0][0]] , twoPaths[1] + [twoPaths[1][0]]]
		
		# 
		
		# return cells, twoPaths
		# else:
			# #vertical to horizontal
	
	def find_adjecentCells_from_two_path_coords(self,coord1, coord2,cells):
		#ASSERT coord1 and coord2 are orthogonal neighbours on lattice
		#sort coords
		adjecent = []
		
		if coord1[0] == coord2[0]:
			
			#horizontal (row the same)
			
			#check sequence and make sure it is normalized (eliminate path direction)
			if coord1[1] >coord2[1]:
				coord1,coord2 = coord2,coord1
			# print "horizzz"	
			# print coord1, coord2
			
			try:
				cells[coord1] 
				adjecent.append(coord1)
			except:
				#ASSERT cell outside  boundaries
				pass
				
			try:
				cells[(coord1[0]-1,coord1[1])] 
				adjecent.append((coord1[0]-1,coord1[1]))
			except:
				#ASSERT cell outside  boundaries
				pass
			
			# if (coord1[0]-1) <0:
				# #edge of lattice
				# return [coord1]
			# else:
				# return [coord1, (coord1[0]-1,coord1[1])]
				
		else:
			#vertical
			#coord1[1] == coord2[1]:
			#check sequence and make sure it is normalized (eliminate path direction)
			
			if coord1[0] > coord2[0]:
				coord1,coord2 = coord2,coord1
			# print "verticalll"
			# print coord1 , coord2
			
			
			try:
				cells[coord1] 
				adjecent.append(coord1)
			except:
				#ASSERT cell outside  boundaries
				pass
				
			try:
				cells[(coord1[0],coord1[1]-1)] 
				adjecent.append((coord1[0],coord1[1]-1))
			except:
				#ASSERT cell outside  boundaries
				pass
			
			# if coord1[1]-1 < 0:
				# return [coord1]
			# else:
				# return [coord1, (coord1[0],coord1[1]-1)]
		return adjecent
	def find_neighbour_cells(self, cell):
		neighbours = []
		#test horizontally
		if cell[1] < self.__cols -2:
			#east
			neighbours.append((cell[0],cell[1]+1) )
		
		if cell[1] > 0:
			#West
			neighbours.append((cell[0],cell[1]-1) )
		
		if cell[0] < self.__rows -2:
			#South
			neighbours.append((cell[0]+1,cell[1]) )
		
		if cell[0] > 0:
			#north
			neighbours.append((cell[0]-1,cell[1]) )	
			
		return neighbours
		
		
	def find_recombination_cells(self, paths, cells,includeOriginalCutCell = False):
		#return all possible cells for recombination for two paths.
		
		possibleRecombinationCells= []
		
		shortestPath = paths[0]
		#make sure shortest path comes first.
		if len(paths[0])>len(paths[1]):
			shortestPath =  paths[1]
		
		end = -1
		if includeOriginalCutCell:
			#if the end of the path is included, the cell where the path was cut (cutting cell) will be discovered again as recombination candidate
			end= None
		
		#find the adjecent cells to the path
		adjecentCells = []
		previousNode= shortestPath[0]
		for node in shortestPath[1:end]:
			adjecentCells.extend(self.find_adjecentCells_from_two_path_coords(node,previousNode,cells))
			# print node
			# print previousNode
			# print '========'
			previousNode = node
			# pass
		
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
	
	def __split_path_loop(self,path,splitA,splitB):
		# ASSERT:splitNode and checkNode are neighbours (direction unknown)
		# path = path[:-1]
		
		# print "faefaestartsplitpath"
		# print path
		# print splitA
		try:
			indexA = path.index(splitA)
			# print "indeces:"
			# print indexA
		except:
			print splitA
			print splitB
			print path
			print '====='
			raise 'fff'
			
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
		
		
	def recombine(self,paths, cell, cells):
		#assert two closed paths in arg. paths
		#connects paths at given cell
		
		#check if node (corresponding to cell (top left node of cell that is) is in first path, if not , switch
		#assume cellId node is in pathA
		pathA = paths[0]
		pathB = paths[1]
		if cell not in pathA:
			pathB,pathA = pathA, pathB
		
		#
		
		#get position of nodes in paths
		#nodeTopLeftPositionInPath = pathA.index(cell)
		#nodeBottomRightPositionInPath = pathB.index((cell[0]+1, cell[1]+1))
		
		neighboursA = self.neighbours_on_cycle(pathA, cell)
		# print "recombinecells:"
		# print neighboursA
		# print cell
		# print "pathA:"
		# print pathA
		# print pathB
		# print cells
		# self.print_cells_ASCII(cells)
		
		isHorizontalReconnection = True
		neighbourToReconnectCellNodeWith = (cell[0],cell[1]+1) #assumehorizontal connection
		neighbourOnSameSide = (cell[0]+1,cell[1])
		if (cell[0],cell[1] +1 ) in neighboursA: #if path is horizontal, there is a vertical reconnection
			neighbourToReconnectCellNodeWith = (cell[0]+1,cell[1])
			neighbourOnSameSide = (cell[0],cell[1]+1)
			isHorizontalReconnection = False
		# print "is Horizontal?:{}".format (str(isHorizontalReconnection))
			
		pathA1, pathA2 = self.__split_path_loop(pathA[:-1], cell, neighbourOnSameSide )
		
		# print "splitfirst loop:"
		# print pathA1
		# print pathA2
		
		pathA = pathA2 + pathA1
		# print pathA
		
		# if isHorizontalReconnection:
			# pathB1, pathB2 = self.__split_path_loop(pathB[:-1], cell, (cell[0]+1, cell[1]))
			
		# print "neighbourToReconnectCellNodeWith"
		# print neighbourToReconnectCellNodeWith
		# print (cell[0]+1,cell[1]+1)
		try:
			pathB1, pathB2 = self.__split_path_loop(pathB[:-1], (cell[0]+1, cell[1]+1), neighbourToReconnectCellNodeWith)
		except:
			print "is Horizontal?:{}".format (str(isHorizontalReconnection))
			print cells
			print cell
			print paths
			self.print_cells_ASCII(cells)
			raise 'wwew'
			
		# else:
			# # print pathB
			# # print (cell[0]+1, cell[1]+1)
			# # print (cell[0], cell[1]+1)
			# pathB1, pathB2 = self.__split_path_loop(pathB[:-1], (cell[0]+1, cell[1]+1), (cell[0]+1, cell[1]))
		
		pathB = pathB2 + pathB1
		
		# print "second path loop:"
		# print pathB1
		# print pathB2
		# print pathB
		
		#generate new cells from path
		recombined = pathA + pathB + [pathA[0]]
		try:
			cells = self.create_cell_pattern_from_hamilton_cycle(recombined)
		except:
			print paths
			print cell
			self.print_cells_ASCII(cells)
			print "eijfijeifjiefj"
		
		# if len(recombined) != 25:
			# raise Exception
		# else:
			# print "0000000000000000000000000"
			# print len(recombined)
		return recombined, cells
			
		
		
		
	# def split_search_return_one_neighbour_or_none(self,path, node):
		# pass
	
	
	
	def neighbours_on_cycle(self,path, node):
		#assert path is cycle (last element equals first in path)
		
		#get index of node in path
		try:
			pos = path.index(node)
		except:
			print "tesijtes"
			print path
			
			raise
		
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
		
	def print_path_ASCII(self,path):
		
		self.print_paths_ASCII([path])
	

	def print_paths_ASCII(self,paths):
		#from a list of vertices, and rows and cols, print path on screen
		#paths as in multiple pats in one lattice! (so, if more than one path, cant be hamilton cycle anymore)
		
		#CREATE EMPTY LATTICE
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
		#pathsDoubled = []
		
		paths = [[(2*r,2*c) for r,c in path] for path in paths]
		
		for path in paths:
			#fill lattice with path
			previousNode = path[0]
			for node in path[1:]:
			
				rowCoord = (previousNode[0] + node[0])/2
				colCoord = (previousNode[1] + node[1])/2
				latticeMinimalCoords[rowCoord][colCoord] = "X"
				previousNode = node
				
		#print lattice
		for printrow in latticeMinimalCoords:
			print "".join(printrow)

def standardized_hamilton_cycle(path):
	#assert hamilton cycle with (row,col) tuples, first and last tuple equal
	#so it always starts with (0,0) and ends with (0,0)
	
	if path[0] ==(0,0):
		return path
	else:
		path = path[:-1]
		i = path.index((0,0))
		return path[i:]+ path[:i] +[(0,0)]
		
		
		

	
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

		
def getNeighbours( path,rows, cols,verbose = False):
	lattice = HamiltonFun(rows, cols)
	tttttcells,splitCells= lattice.find_split_cells(path)
	if verbose:
		print "with potentials:"
		lattice.print_cells_ASCII(tttttcells)

	neighbours = []
	
	for splitCell in splitCells:
		if verbose:
			print "splitcell and lattice to split: {}".format(splitCell)
			lattice.print_cells_ASCII(tttttcells)
			
		cells,splitPaths =  lattice.split(tttttpath[:], splitCell ,tttttcells.copy())
		
		if verbose:
			print "split done:"
			lattice.print_cells_ASCII(cells)
			lattice.print_paths_ASCII(splitPaths)
			print splitPaths
			print "origpath"
			print path
		try:
			recombinationCandidates, cellsb = lattice.find_recombination_cells(splitPaths, cells)
		except:
			print "failed recomb"
			
			lattice.print_cells_ASCII(cells)
			raise
		
		if verbose:
			# print "recomb candidates:"
			# print recombinationCandidates
			print "after recomb cand:"
			lattice.print_cells_ASCII(cellsb)
			
			# print cellsb
			# print lattice.is_pattern_hamilton_cycle(cellsb)
			
		
		for cand in recombinationCandidates:
			# print "candidate {} :".format(str(cand))
			# print splitPaths
			try:
				newPath, cellsc =  lattice.recombine(splitPaths, cand, cells.copy())
				
				# print "aidjfiajsdfjja PATH LENGTH"
				# print len(pathb)
				
			except:
				lattice.print_cells_ASCII(cells)
				raise
			
				
			neighbours.append(standardized_hamilton_cycle(newPath))
			
			if verbose:
				print "recombined:"
				# print path
				lattice.print_cells_ASCII(cellsc)
	
	return neighbours
	
if __name__ == "__main__":
	ROWS = 4
	COLS = 6
	ITERATIONS = 100
	lattice = HamiltonFun(ROWS, COLS)
	
	# dictje = lattice.latticeGraphDictGenerator()
	# hoi = graphs.Graph(dictje)
	# print hoi
	# print hoi.vertices()
	
	
	# print lattice.is_hamilton_cycle_existing()
	# tttttpath = lattice.set_up_first_cycle()
	# [[(2, 1), (2, 1), (2, 0), (3, 0), (3, 1), (2, 1)], [(3, 2), (3, 3), (3, 4), (3,
# 5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (2, 3), (1, 3), (0,
# 3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2), (3, 2)]]
	# tttttpath = [(2, 1), (2, 1), (2, 0), (3, 0), (3, 1), (2, 1),(3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (2, 3), (1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2), (3, 2)]
	tttttpath = [ (2, 1), (2, 0), (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2), (2,1)]
	
	print lattice.print_path_ASCII(tttttpath)
	
	neighbours = []
	
	for iter in range(ITERATIONS):
		try:
			n = getNeighbours(tttttpath,ROWS,COLS, False)
			neighbours.extend(n)
			
			
		except:
			print "wrong wrong wrong"
			# print lattice.print_cells_ASCII(cellsc)
			print lattice.print_path_ASCII(tttttpath)
			raise 
	
		
		print "----------------------------"
		
		tttttpath = random.choice(neighbours)
		#tttttpath = neighbours[0]
		# try:
			# tttttpath = neighbours[0]
		# except:
			# print neighbours
			# raise
	print len(neighbours)
	neighbour_paths_no_doubles = [list(x) for x in set(tuple(x) for x in neighbours)]
	print len(neighbour_paths_no_doubles)
	for neighbour in neighbour_paths_no_doubles:
		lattice.print_path_ASCII(neighbour)
		print "\n"
		
		# print standardized_hamilton_cycle(neighbour	)