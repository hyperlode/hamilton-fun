import graphs
import fileOperations

# class HamiltonFun:
	 # def __init__(self, lattice_rows, lattice_cols):
		
        # """ initializes a graph object """
        # self.__graph_dict = graph_dict
	

def latticeGraphDictGenerator(rows, cols):
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
	for row in range(rows):
		for col in range(cols):
			neighbours = []
			north = row -1
			if north >= 0:
				neighbours.append((north,col))
			south = row +1
			if south < rows:
				neighbours.append((south,col))
			east = col + 1
			if east < cols:
				neighbours.append((row,east))
			west = col - 1
			if west >= 0:
				neighbours.append((row,west))
			
				
			latticeGraph[(row,col)]=neighbours
	
	return latticeGraph	

def find_path(graph):
	print graph.find_path((0,0),(0,1))
	
def find_all_paths(graph):
	paths = graph.find_all_paths((0,0),(0,1)) 
	hamiltonPaths= []
	for path in paths:
		if len(path) == len(graph.vertices()):
			hamiltonPaths.append( path)
	return hamiltonPaths
	
def print_path(path,rows, cols):
	#from a list of vertices, and rows and cols, print path on screen
	
	#create row with points on each node
	print_empty_row = cols*"X"
	print_empty_row = list(" ".join(print_empty_row))
	
	#create completely blank row
	blankRowPrint =  list(((2*cols-1)*" "))
	
	#add all rows to lattice
	latticeMinimalCoords = []
	for row in range(rows):
		latticeMinimalCoords.append(print_empty_row[:])
		latticeMinimalCoords.append(blankRowPrint[:])
	latticeMinimalCoords.pop()
	
	#make closed path by repeating first element as last position
	path.append(path[0])
	
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
	rows = 6
	cols = 5
	lattice =  latticeGraphDictGenerator(rows,cols)
	lode = graphs.Graph(lattice)
	# print lode
	# find_path(lode)
	paths = find_all_paths(lode)
	# print paths[0]
	# print_path(paths[0],rows,cols)
	
	
	for path in paths:
		print_path(path,rows,cols)
		print "\n"