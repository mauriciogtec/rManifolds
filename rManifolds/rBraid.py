import  random, operator, math, spherogram, numpy
class rBraid(): 
    ""
    def __init__(self,n,k,filename=""):	
      	# Author: Mauricio B. Garcia Tec
        # Date: February 2014
        # In:
        #     n.- the number of strands in the braid group.
        #     k.- the length of the random walk.
        # Out:
        #     A knot or link which is the closure of the random braid   
        random.seed()

        const = max(1,k/25.0) # a posteriori resizes the knot conveniently, controls the height/width ratio.
        
        # Bottom of the braid.
        strand_pos = [range(n)]+[range(n)] # the vertex position of the i-th strand.
        swap, swapinv = range(n), range(n) # auxiliaries for computations with permutations
        height = [[ (s+1) for s in range(n)]] + [[n+1 for s in range(n)]] 
        crossings = [] # Pairs indicating the positions of the n crossings. Left edge passes under the right edge.
        cross_sph_list = [ [] for s in range(n)] # This is the memory intensive part. We will create a list of length n. 
			   # Each element is a sublist of triples [a,b,c] where a \in {0,1,2,...,n-1} indicates the vertices touched by each strand and 
                           # and b \in {-1,1} indicates whether the flow of the strand at that vertex is above (1) or below (-1).
        cross_orientation = []  # vector with values {0,1} indication type of crossing orientation, right over left (0), otherwise (1)
        
        # The random walk.
        for i in range(k): # the random walk.
            gen = random.randint(1,2*(n-1)) # each number represents with equal probability a generator of the braid group.
                                            # each pair is a generator and its inverse.
            pos = (gen-1)/2 # connects the pos-th and (pos+1)-th strands.
            above = (gen-1)%2 # tells which strand passes above. If 1 then left over right, 0 otherwise.
            height += [[(n + 2 + i) for s in range(n)]] # for simplicity in the resulting code, we will generate 2*k levels during the walk.
            ### strand position (efficiency to be improved).
            swap_aux = swap[pos]
            swap[pos] = swap[pos+1]
            swap[pos+1] = swap_aux
            swap_inv = [[swap[j],j] for j in range(n)]
            swap_inv.sort(key=operator.itemgetter(0))
            swap_inv = [row[1] for row in swap_inv]
            ### 
            strand_pos += [[s for s in swap_inv]]
            cross1 = swap[pos]
            cross2 = swap[pos+1]
            if above == 0:
                crossings += [[cross1,cross2]]
                cross_sph_list[swap[pos]] += [[i,-1]]
                cross_sph_list[swap[pos+1]] += [[i,1]]
            else:
                crossings += [[cross2,cross1]]
                cross_sph_list[swap[pos]] += [[i,1]]
                cross_sph_list[swap[pos+1]] += [[i,-1]]
            cross_orientation += [above]

        # Top of the braid.
        strand_pos += [swap_inv] 
        height += [[(2*n+k+1-swap_inv[s]) for s in range(n)]]      
            
        # Detects the number of components factoring the last permutation; each cycle stands for a component. Could also be improved.
        remaining = {s for s in range(n)}
        perm = [swap_inv[i] for i in range(n)]
        cycles = [] 
        subcycle = []
        place,smallest = 0,0
        while remaining!=set():
            if place in remaining:
                remaining.remove(place)
                auxplace = place
                place = perm[auxplace]
                subcycle += [place]
            else:
                cycles += [[s for s in subcycle]]
                smallest += 1
                while not smallest in remaining and smallest < n-1:
                    smallest += 1                    
                place = smallest
                subcycle = []
        cycles += [[s for s in subcycle]] # is the cycle factorization of the permutation.
        linkcomponents = len(cycles)
        self.num_components = linkcomponents
        
        # Note on vertex counting: to each strand we will add 2 points to close the strand. Each strand has
        #   k+3 vertices as we defined them, that makes k+5 vertices.
            
        # With the generated braid we will I) build a .lnk file a store it there II) create an object using the spherogram module

	# I) We will print the object into 'filename'.lnk for the PLink editor. We store the string into 'output'.
        output = '% Link Projection \n'+str(linkcomponents)+'\n'
        output += '   0 0\n'
        cycle_pos = [len(cycles[0])]
        
        for t in range(linkcomponents-1): 
            next = t+1
            cycle_pos += [cycle_pos[t] + len(cycles[next])]
            
        for comp in range(linkcomponents-1): # link component closure and opening.
            output += '   '+str(cycle_pos[comp]*(k+5))+' '+str(cycle_pos[comp]*(k+5))+'\n'
                
        # Printing the vertices          
        vertexnumber = n*(k+5)
        output += str(vertexnumber)+'\n'
        perm = [s for sublist in cycles for s in sublist] # we print vertices in cycle order, in reverse order as swap.
        
        for i in range(n):
            x = [const*row[perm[i]] for row in strand_pos]
            y = [row[perm[i]] for row in height]
            lastx = x[k+2]/const # scaling for improving visualization
            # ---- now we add two points for closing the braid so that we can close the braid.
            x += [2*const*n-lastx] + [2*const*n-lastx] 
            y += [2*n+k+1-lastx] + [lastx+1]
            for j in range(k+5):
                output += '   '+str(x[j])+' '+str(-y[j])+'\n'
                    
        # Printing the edges          
        output += str(vertexnumber)+'\n'
        cycle_pos = [0] + [(k+5)*s for s in cycle_pos] 
        
        for s in range(linkcomponents):
            for t in range(cycle_pos[s],cycle_pos[s+1]):
                output += '   '+str(t)+' '+str((t+1)%cycle_pos[s+1]+cycle_pos[s]*((t+1)/cycle_pos[s+1]))+'\n'
                    
        # Printing the crosses
        output += str(k)+'\n'
        # Luckily, we know that there are exactly k crossings and where they take place.         
        perm_inv = [[perm[i],i] for i in range(n)]
        perm_inv.sort(key=operator.itemgetter(0))
        perm_inv = [row[1] for row in perm_inv] 
        # Invert permutation to rewrite crossings in terms of edge numbers. For example, if 
        # there is a cross between 2 and 1 but 2 was printed first, then we want pinv so that pinv(2)=1.    
        crossings = [[(k+5)*perm_inv[s[0]],(k+5)*perm_inv[s[1]]] for s in crossings]
        crossings = [[crossings[i][0]+i+1,crossings[i][1]+i+1] for i in range(k)]

        for s in range(k):
            output += '   '+str(crossings[s][0])+' '+str(crossings[s][1])+'\n'

        output += str(-1)	
    
        if (filename!=""):
            namelen = len(filename)
            if (namelen >4):
                if (filename[(namelen-4):namelen]!=".lnk"):
                    filename = filename + ".lnk"
            else:
                filename = filename + ".lnk"
            self.filename = filename
            with open(filename,'wt') as file_output:
                file_output.write(output)
            print('A file named "'+filename+'" was created in the current working directory.')

        # II) To create an object directly into SnapPy we use the spheregram module convention.
        Sph_Crossings = [spherogram.links.links.Crossing(j) for j in range(k)]
        # we will glue cross_sph_list according to observed cycles.
        for i in range(linkcomponents):
            current_cycle = cycles[i]
            len_current_cycle = len(current_cycle)
            cross_list = []
            for s in range(len_current_cycle):
                cross_list += cross_sph_list[current_cycle[s]]
            N = len(cross_list)
            for j in range(N): # we have to proceed by cases.
                vertex = cross_list[j]
                if vertex[1] == -1: # if present vertex flows above or below.
                    prev = 0
                    proc = 2
                else:
                    if cross_orientation[vertex[0]] == 0: # left over right type
                        prev = 3
                        proc = 1
                    else:
                        prev = 1
                        proc = 3
                # 1.- Connect to previous in cycle.
                prev_vertex = cross_list[(j-1)%N]
                if prev_vertex[1] == -1: # if previous flows below, it is easy, it connects to 2.
                    Sph_Crossings[vertex[0]][prev] = Sph_Crossings[prev_vertex[0]][2]
                else: # we have to detect orientation of previous.
                    if cross_orientation[prev_vertex[0]] == 0: # previous is left over right
                        Sph_Crossings[vertex[0]][prev] = Sph_Crossings[prev_vertex[0]][1]
                    else:
                        Sph_Crossings[vertex[0]][prev] = Sph_Crossings[prev_vertex[0]][3]
                # 2.- Connect to proceeding in cycle.
                next_vertex = cross_list[(j+1)%N]
                if  next_vertex[1] == -1: # if proceeding flows below, it is easy, it connects to 2.
                    Sph_Crossings[vertex[0]][proc] = Sph_Crossings[next_vertex[0]][0]
                else: # we have to detect orientation of previous.
                    if cross_orientation[next_vertex[0]] == 0: # proceeding is right over left
                        Sph_Crossings[vertex[0]][proc] = Sph_Crossings[next_vertex[0]][3]
                    else:
                        Sph_Crossings[vertex[0]][proc] = Sph_Crossings[next_vertex[0]][1]                
                            
        # We know store the complement of the braid into a manifold object.
        self.manifold = spherogram.links.links.Link(Sph_Crossings).exterior()
