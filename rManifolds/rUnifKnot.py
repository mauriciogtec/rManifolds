import  random, operator, math, spherogram, numpy
class rUnifKnot(): 
    ""
    def __init__(self,n,filename=""):
    	# Author: Mauricio B. Garcia Tec.
    	# Date: Februrary 2014.
    	# This function generates a random knot in SnapPy using the module Spherogram

        # In: int n.-the number of desired random vertices.
        # Out: a object of type Link in SnapPy
        random.seed() # iniatilize random number generator
        self.num_components = 1
    
        # Generation of n random vertices.
        const = 500 # So that the knot can be visualized in the PLink editor witout automatic scaling
        x = [const*random.random() for i in range(n)]	
        y = [const*random.random() for i in range(n)]
        z = [const*random.random() for i in range(n)]
    
        # 1. This information shows that a knot is to be saved (a one-component link).
        output = '% Link Projection \n1 \n   0    0'	
    
        # 2. Prints vertices information.
        output += '\n'+str(n)
        for i in range(n):
            output += '\n   '+str(x[i])+'  '+str(y[i])
        # 3. Connectivity information in standard form.
        output += '\n'+str(n)
        for i in range(n):
            output += '\n    '+ str((i+1) % n) +'  '+str(i)
        
        # 1. Crossings: 
        # CHALLENGE: change from vertex information to crossings information.
        # We shall add vertex by vertex, looking for crossings with the ones already added.
        crossinfo = [] # used for storing crossing information. A list with 6-tuples containing
		       # (a) an edge, (b) another edge (each crossing pair appears twice)
		       # (c) 1 if first edge goes above, -1 if second index goes above
		       # (d) relative distance from the crossing point to the origin of the first-edge
		       # (e) orientation of the crossing: 1 if the edge passing above is going to the right else -1
		       # (f) (pot) the coordinates of the point
		       # (g) crossid (see below)
        crossid = 0 # an id for crossings since each crossing appears twice.
    
        # The following two will be used when printing the .lnk file
        under = [] # 
        above = [] # 
    
        for i in range(n):
            # We add the i-th edge and then check for crossings with the ones already added. We use %n% since it is a circular list.
            a0,b0,c0,a1,b1,c1 = x[i],y[i],z[i],x[(i+1)%n],y[(i+1)%n],z[(i+1)%n] 
            for j in range(i-1):
                # Projection onto the x an y coordinates and use a parametric model to verify there is a crossing.
                # The formula in use comes from solving the system of two linear equations on t and u:  
                #	                (x0,y0)+t(x1-x0,y1-y0) = (a0,b0)+u(a1-a0,b1-b0).
                # If 0 <= t,u <=1, then there is a crossing between two line segements.
                x0,y0,z0,x1,y1,z1 = x[j],y[j],z[j],x[(j+1)%n],y[(j+1)%n],z[(j+1)%n]  
                A,B,C,D = a1-a0,-(x1-x0),b1-b0,-(y1-y0)
                det = A*D-B*C
                t = (D*(x0-a0)-B*(y0-b0))/det
                u = (-C*(x0-a0)+A*(y0-b0))/det	
                if 0<t and t<1 and 0<u and u<1:			
                    if c0+t*(c1-c0) > z0+u*(z1-z0):
                        aux = 1 # first index goes above.
                        under += [j]
                        above += [i]
                    else: 	
                        aux = -1 # second-index goes above.
                        under += [i]
                        above += [j]
                    # To know whether the edge above is going to the right or to the left with respect to edge below
                    # we will use the determinant.
                    # convention: 1 = the above edge is to the right with respect to the edge below	
                    #            -1 = the above edge is to the left with respect to the edge below 
                    coords = [x0+u*(x1-x0),y0+u*(y1-y0),z0+u*(z1-z0)] # crossing coordinates (if crossing)
                    signum = math.copysign(1,(x1-x0)*(b1-b0)-(a1-a0)*(y1-y0))
                    crossinfo += [(j,i,aux,u,aux*signum,coords,crossid)] # the smaller the u, the closest it is to the vertex vec[j]
                    crossinfo += [(i,j,-aux,t,aux*signum,coords,crossid)] # the smaller the t, the closest it is to the vertex vec[i]
                    # Saving the t and u allows us to recreate the crossing information while going around the knot.
                    crossid += 1 # one more crossing		
        N = len(crossinfo) # vector length.
        totcross = N/2 # total number of crossings found.	
        # We want to sort the crossings as if we were traveling the from one vertex to the next in order and going around the knot.
        crossinfo.sort(key=operator.itemgetter(0,3)) #  (Nice stuff!).
        crossid = [row[6] for row in crossinfo] # we extract the crossid column for simplicit
    
        # Print the .lnk file for the PLink editor.
        output += '\n'+str(totcross)
        for i in range(totcross):
            output += '\n   '+str(under[i])+'    '+str(above[i])
        output += '\n'+str(-1)	
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
    
        # Use the Spherogram module instead.
        crossings = [spherogram.links.links.Crossing(j) for j in range(totcross)] # to start a variable for crossing conventions
        # We need to mark each crossing
        k = 0
        while k < N:
		# We will proceed as follows:
		# For each crossing in visit order (crossinfo) we will have two cases: whether the flows is going above or below.
		# In both cases we look at the previous and the following knot and with that determine how to tie the crossings.
		# The convention is that each crossing has entries 0,1,2,3 where 0 corresponds to the edge incoming from below
		# having the lowest index. Entries 1,2,3 and assigned moving counterclockwise.
		# There will be 3 interesting cases. If the flow goes below, if the flow goes above in standard direction (the above
		# edge passes to the right) and if it goes above in skew direction (to the left).
		# FIRST CASE	
		# We know in this case that the previous crossing in the list will be connected to  crossings[k][0] and 
		# the one that follows to crossings[k][2]. We need now to check the crossing type of the previous and
		# the following to know to which of their entries to connect.
            if crossinfo[k][2]==-1:
                prev = 0
                proc = 2
                # SECOND CASE	
                # If above then the sign will determine which entry to connect.
            else: 
                if crossinfo[k][4]==1:
                    prev = 3			
                    proc = 1
                    # THIRD CASE
                else:
                    prev = 1
                    proc = 3
            # Now we connect to the previous and next crossings again depending on their type.
            if crossinfo[(k-1)%N][2]==-1: # Case a) The previous crossing's flow is below.
                crossings[crossid[k]][prev] = crossings[crossid[(k-1)%N]][2] 		
            else: # Case b) If the flow is going above, we need to know we must look at the signum (left or right).
                if crossinfo[(k-1)%N][4] == 1: #  Case i) standard direction, the above edge passes to the right.
                    crossings[crossid[k]][prev] = crossings[crossid[(k-1)%N]][1] 
                else: # Case ii) skew direction
                    crossings[crossid[k]][prev] = crossings[crossid[(k-1)%N]][3]
            # The instructions are similar for the following cross.
            if crossinfo[(k+1)%N][2]==-1: # Case a) The follwing crossing's flow is below.
                crossings[crossid[k]][proc] = crossings[crossid[(k+1)%N]][0] 		
            else: # Case b) Flow goes above, we look at signum
                if crossinfo[(k+1)%N][4] == 1: # Case i) standard direction, the above edge passes to the right.
                    crossings[crossid[k]][proc] = crossings[crossid[(k+1)%N]][3] 
                else: # Case ii) skew direction
                    crossings[crossid[k]][proc] = crossings[crossid[(k+1)%N]][1]						
            k += 1	
        # Now we craft the knot.
        self.manifold = spherogram.links.links.Link(crossings).exterior()
