import graphs
import fileOperations

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

if __name__ == "__main__":

	print latticeGraphDictGenerator(3,3)
	# lode = graphs.Graph()