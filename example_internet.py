from itertools import combinations as comb
# http://arthurvause.blogspot.ca/2013/05/hamiltonian-cycles-on-square-grid-part-2.html  
def graph(m): # define a rectangular grid graph m by m
  g = { n:set([n+1,n-m,n+m]) & set(xrange(m*m)) for n in xrange(0,m*m,m)}
  g.update( { n:set([n-1,n+1,n-m,n+m]) & set(xrange(m*m))
            for n in xrange(m*m) if 0 < n%m < m-1  })
  g.update({ n:set([n-1,n-m,n+m]) & set(xrange(m*m))
            for n in xrange(m-1,m*m,m)})
  return g

def printgrid(outside,M):
  s="\n"
  for i in xrange(M):
    for j in xrange(M):
      s+= "  " if M*i+j in outside else "* "
    s+="\n"
  print s

def find_inside_connected(graph, maxnode, start, pathset,grid):
  pathset.add(start)
  for node in graph[start]:
    if node < maxnode and grid[node] and node not in pathset:
      find_inside_connected(graph, maxnode, node, pathset, grid)

def find_outside_connected(graph, start, pathset,grid):
  pathset.add(start)
  for node in graph[start]:
    if not(grid[node]) and node not in pathset:
      find_outside_connected(graph, node, pathset,grid)

def valid(gr,M,grid): # check that there are (M-1)^2 / 2 outside cells connected to perimeter
  global perimeter
  outside_connected=set()
  for s in perimeter:
    if not(grid[s]):
      find_outside_connected(gr, s, outside_connected, grid)

  return len(outside_connected)==((M-1)*(M-1))/2

def hamiltonian(outside,M):
  inside = set(range(M*M))-set(outside)
  # convert to an M+1 by M+1 array to trace the Hamiltonian
  converted = set([x+x/M for x in inside])
  N=M+1
  ham=[0,1]
  pos=1
  move=1
  while pos != 0:

    if move==1: # going right, 
      if pos-N in converted: move=-N
      elif pos in converted: move=+1
      else: move=+N  
        
    elif move==N: # going down, 
      if pos in converted: move=1
      elif pos-1 in converted: move=+N
      else: move=-1  

    elif move== -N: # going up, 
      if pos-N-1 in converted: move=-1
      elif pos-N in converted: move=-N
      else: move=1  

    elif move== -1: # going left, 
      if pos-1 in converted: move=N
      elif pos-N-1 in converted: move=-1
      else: move=-N  
        
    else:
      raise Exception("no more cases")
     
    pos+=move  
    ham.append(pos)

  if len(ham) != N*N+1 or len(set(ham)) != N*N:
    print ham
    raise Exception("Invalid Hamiltonian")

  return ham

# powerset of an iterable, up to L elements (adapted from itertools recipes)
from itertools import chain, combinations, product
def powerset(iterable, L=None): 
  if L == None: L = len(iterable)
  s = list(iterable)
  return chain.from_iterable(combinations(s, r) for r in range(L+1))

def findHamiltonians(gr, grid,M, x, numoutside ):
  global MAXOUTSIDE, MAXINSIDE, BOTTOMLEFT, BOTTOMRIGHT

  if x-numoutside > MAXINSIDE:
    return
  if numoutside  > MAXOUTSIDE:
    return
  if x==M*M: # all cells populated, check for a valid configuration
    if valid(gr,M,grid):
      outside = set([i for i in xrange(M*M) if not grid[i]])
      yield hamiltonian(outside,M)
    return

  if x == BOTTOMLEFT: # bottom corners have to be inside
    options=[True]
  elif x==BOTTOMRIGHT: # bottom corners have to be inside
    if grid[x-M-1] and grid[x-1]==grid[x-M]:
      return
    else:
      options=[True]
  elif x > BOTTOMLEFT and not grid[x-1]: # no adjacent outside cells on bottom row
    if grid[x-M-1] and not(grid[x-M]):
      return
    else:
      options=[True]
  elif x%M == 0: # start of a new row

    # check that all inside cells are connected to last completed row
    inside_connected = set(i for i in xrange(x-M,x) if grid[i])
    for s in xrange(x-M,x):
      if grid[s]:
        find_inside_connected(gr, x, s, inside_connected, grid)
    if inside_connected <>  set(i for i in xrange(x) if grid[i]):
      return

    if not(grid[x-M]) : # no adjacent outside cells on left column
      options=[True]
    else:
      options=[True,False]

  #elif x%M == M-1 and not(grid[x-M]) and not(grid[x-1]) and grid[x-M-1] :
  #    return

  else:  # avoid an invalid 2 by 2 configuration   
    if grid[x-1]^grid[x-M]:
      options=[True,False]
    else:
      options = [ not(grid[x-M-1]) ]

  for op in options:
    grid[x]=op
    for f in findHamiltonians(gr,grid,M, x+1, numoutside + (0 if grid[x] else 1)):
      yield f       

perimeter=set()

def Hamiltonians(N):
  global perimeter, MAXOUTSIDE, MAXINSIDE, BOTTOMLEFT, BOTTOMRIGHT
  M=N-1
  MAXOUTSIDE = (M-1)*(M-1)/2
  MAXINSIDE = (M*M+2*M-1)/2
  BOTTOMLEFT = M*(M-1)
  BOTTOMRIGHT = M*M-1
  gr=graph(M)
  noncorners = set(range(M*M))-set((0,M-1,M*(M-1),M*M-1))
  allnodes=set(range(M*M))
  centre = set([x for x in xrange(M*M) if x/M not in (0,M-1) and x%M not in (0,M-1)])
  perimeter = allnodes-centre

  for p in powerset( range(1,M-1), (M-1)/2 ):
    if all( p[i+1] !=p[i]+1 for i in xrange(len(p)-1)):
      grid = [True]*(M*M)
      for i in p:
        grid[i]=False
      for f in findHamiltonians(gr,grid,M,M, sum([1 for x in xrange(0,M) if not(grid[x])]) ):
        yield f

if __name__=="__main__":
  from time import clock
  start= clock()
  for N in (4,6):
    start=clock()
    solutions = [h for h in Hamiltonians(N)]
    interval=clock()-start
    print "N=",N,":", len(solutions),"solutions,",interval, "sec"
    print "N=",N,":", len(set([tuple(x) for x in solutions]))  , "distinct solutions" 

  start=clock()
  counter=0
  for h in Hamiltonians(8):
    counter+=1
    if counter%10000 == 0:
      print counter,"Hamiltonians found,", clock()-start,"sec"
  print counter," Hamiltonians for N=8,", clock()-start, "sec"