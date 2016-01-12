import graphs
import fileOperations

class HamiltonFun:
	def __init__(self, lattice_rows, lattice_cols):
		""" initializes a graph object """
		self.__rows = lattice_rows
		self.__cols = lattice_cols
		self.__latticeGraphDict = self.latticeGraphDictGenerator()
		self.__latticeGraph = graphs.Graph(self.__latticeGraphDict)
	
	
	
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
	
	def find_all_Hamilton_cycles(self):
		paths = self.__latticeGraph.find_all_paths((0,0),(0,1)) 
		self.__hamiltonCycles= []
		for path in paths:
			if len(path) == len(self.__latticeGraph.vertices()):
				path.append(path[0]) #close the path
				self.__hamiltonCycles.append( path)
		return self.__hamiltonCycles
	
		
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

if __name__ == "__main__":
	rows = 4
	cols = 4

	lattice = HamiltonFun(rows,cols)
	cycles = lattice.find_all_Hamilton_cycles()
	lattice.print_path_ASCII(cycles[0])