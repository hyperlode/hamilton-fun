import graphs
import fileOperations
import random
from collections import defaultdict




class Lattice:
	def __init__(self, lattice_rows, lattice_cols):
		""" initializes a graph object """
		self.__rows = lattice_rows
		self.__cols = lattice_cols
		# self.__rowCells = self.__rows-1
		# self.__colCells = self.__cols-1
		self.__latticeGraphDict = {}
		self.__create_lattice_graph_dict()
		self.__latticeGraph = graphs.Graph(self.__latticeGraphDict)
		
		self.__paths = []  #list of paths in the lattice. 
		self.__classification = {}
		# self.__initCycle = []
		# self.__activePath = []
		# self.__activeCells = []
		# self.__splitPoints = []
		# self.__neighbourCycles = []
		# self.__splitCells = []
	
	def rows(self):
		return self.__rows
	def cols(self):
		return self.__cols
	def __create_lattice_graph_dict(self):
		#def latticeGraphDictGenerator(self):
		#of this type, where the node name is a coordinate: (row,col)
		# g = { "a" : ["d","f"],
		   # "b" : ["c","b"],
		   # "c" : ["b", "c", "d", "e"],
		   # "d" : ["a", "c"],
		   # "e" : ["c"],
		   # "f" : ["a"]
		# }

		#lattice: every node is connected (edges) to its neighbours
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
				self.__latticeGraphDict[(row,col)]=neighbours
		
		return True
		
	def find_neighbour_nodes(self, node):
		#neighbours in lattice as list
		return self.__latticeGraphDict[node]
	
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
	def stats(self):
		return self.__classification
		
	def new_paths_to_lattice(self,paths):
		#feed with single path as a list, or a list of paths.  paths can be cycles or paths.
		self.__paths = []
		
		#convert to double list if single path given.
		if not isinstance(paths[0], list):
			paths = [paths]
			
		if self.checkAndClassifyPaths(paths):
			self.__paths = paths
		else:
			#given path not valid
			raise
		return 
	
	def checkAndClassifyPaths(self, paths):
		
		#not tested for this lattice, only if it "could" fit in this lattice. 
		#multiple paths should not overlap. paths can be path or cycle. individual path should not overlap. every element of path should go to neighbour.
		
		classification = {"cycles":0, "numberOfPaths":0, "hamilton":False,"error":False,"visitedNodes":0,"allNodesVisited":False}
		
		visitedNodes = []
		
		#number of paths
		classification["numberOfPaths"]= len(paths)
		
		#number of loops
		for path in paths:
		
			
			#check if paths are loops
			if path[0] == path[-1]:
				classification["cycles"] += 1
				classification["visitedNodes"] += len(path)-1 #cycle has repetition of beginning node on end (don't count twice)
				visitedNodes.extend(path[:-1])
			else:
				classification["visitedNodes"] += len(path)  
				visitedNodes.extend(path)		#add all visited nodes
			
			#check if valid neighbours in node string
			previousNode = path[0]
			
			for node in path[1:]:
				
				if node not in self.find_neighbour_nodes(previousNode):
					classification["error"] = "node {} not neighbour of {} in path".format(str(node),str(previousNode))
					print self.find_neighbour_nodes(previousNode)
					print classification
					print path
					
					raise Exception("no valid neighbour, does node exists anyways? see error in classification dict")
				previousNode = node
				#add 1 to number of nodes visited
				
				
			
		#check if nodes are visited multiple times
		if len(visitedNodes) != len(set(visitedNodes)):
			classification["error"] = "nodes are visited multiple times"
			print paths
			print classification
			raise Exception ("nodes visited multiple times")
		
		#check if not too many nodes are visited
		if classification["visitedNodes"] > self.__rows * self.__cols :
			classification["error"] = "too many nodes visited number of nodes visited: {}".format(classification["visitedNodes"])
			print paths
			print classification
			raise Exception ("too many nodes visited")
			
		#check if all nodes are visited
		if classification["visitedNodes"] == self.__rows * self.__cols:
			classification["allNodesVisited"] = True
		
		#check if Hamilton path or Cycle
		if classification["allNodesVisited"] and classification["numberOfPaths"] == 1:
			classification["hamilton"] = True
			
		
		
		self.__classification = classification

		return not classification["error"]
	
	def __str__(self):
		return fileOperations.dictToString(self.__classification)
		# return str(self.__classification)
	
	def string_from_paths(self):
	
		#from a list of vertices, and rows and cols, print path on screen
		#paths as in multiple pats in one lattice! (so, if more than one path, cant be hamilton cycle anymore)
		
		paths = self.__paths
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
			# print path
			#fill lattice with path
			previousNode = path[0]
			# print previousNode
			for node in path[1:]:
			
				rowCoord = (previousNode[0] + node[0])/2
				colCoord = (previousNode[1] + node[1])/2
				latticeMinimalCoords[rowCoord][colCoord] = "X"
				previousNode = node
		#print lattice
		for printrow in latticeMinimalCoords:
			print "".join(printrow)
		print "\n"

# class HamiltonCycle():
	
	


			
if __name__== "__main__":
	test = Lattice(4,6)
	# path = [ (2, 1), (2, 0), (3, 0),(3, 1), (3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2), (2,1)]
	path = [ (2, 1), (2, 0), (3, 0),(2,1)],[(3, 1), (3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2), (2,1),(3,1)] #(2,1) node at end multiple visits
	path = [ (2, 1), (2, 0), (3, 0),(2,1)],[(3, 1), (3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2),(3,1)]    #invalid path, has diagonal in loops...
	path = [ (2, 1), (2, 0), (3, 0),(3,1),(2,1)],[(3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2),(3,2)]    #valid path
	path = [ (2, 1), (2, 0), (3, 0),(3,1),(3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2),(2,1)]    #hamilton cycle
	path = [ (2, 1), (2, 0), (3, 0),(3,1),(3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3), (1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2)]    #hamilton path
	path = [ (2, 1), (2, 0), (3, 0),(3,1),(3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3)],[(1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2)]    #two paths
	path = [ (3, 0), (2, 0), (2, 1) ],[(3,1),(3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3)],[(1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2)]    #three paths, one inverted direction
	path = [ (3, 0), (2, 0), (2, 1),(3,1),(3,0)],[(3, 2), (3, 3), (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3)],[(1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2)]    #two paths , one cycle
	path = [ (3, 0), (2, 0), (2, 1),(3,1),(3,0)],[ (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3)],[(1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2)]    #two paths , one cycle   not all nodes visited
	path = [ (3, 0), (2, 0), (2, 1),(3,1),(3,0)],[ (3, 4), (3,5), (2, 5), (1, 5), (0, 5),  (1, 4), (2, 4), (2, 3)],[(1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (1, 1), (1, 2), (2, 2)]    #two paths , one cycle   not all nodes visited   neighbour missing
	path = [ (3, 0), (2, 0), (2, 1),(3,1),(3,0)],[ (3, 4), (3,5), (2, 5), (1, 5), (0, 5), (0, 4), (1, 4), (2, 4), (2, 3)],[(1, 3), (0,3), (0, 2), (0, 1), (0, 0), (1, 0), (10, 1), (1, 2), (2, 2)]    #two paths , one cycle   not all nodes visited, non existing node
	
	print len(path)
	print test.new_paths_to_lattice(path)
	print test.string_from_paths()
	print str(test)