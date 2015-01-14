import  random, operator, math, spherogram, numpy
class rGaussianKnot(): 
    ""
    def __init__(self, n, filename="", sigma=1):
    	# Author: Mauricio B. Garcia Tec.
    	# Date: Februrary 2014.
    	# This function generates a random knot in SnapPy using the module Spherogram

        # In: int n.-the number of desired random vertices.
        # Out: a object of type Link in SnapPy

        # Generation of n random points following a gaussian distribution and given that X_0=X_(n+1)=0
        x,y,z = [0],[0],[0]        
        for i in range(1,n):
            cond_sd = sigma*math.sqrt((n-i)/(1.0*(n+1-i))) 
            # Random walk given that X_n = 0. Conditional variance is (s/t)*(t-s) where t is the time at 
            # the known future state and current time. By past independence, we can set t = n+1-i and s = 1
            cond_mean_x = -x[i-1]/(1.0*(n+1-i))
            cond_mean_y = -y[i-1]/(1.0*(n+1-i))
            cond_mean_z = -x[i-1]/(1.0*(n+1-i))
            x += [x[i-1] + numpy.random.normal(cond_mean_x,cond_sd,1)[0]]	
            y += [y[i-1] + numpy.random.normal(cond_mean_y,cond_sd,1)[0]]
            z += [z[i-1] + numpy.random.normal(cond_mean_z,cond_sd,1)[0]]

        # The rest of the algorithm is the same as in the uniform case. I will code comments.
        output = '% Link Projection \n1 \n   0    0'	  
        output += '\n'+str(n)
        for i in range(n):
            output += '\n   '+str(x[i])+'  '+str(y[i])
        output += '\n'+str(n)
        for i in range(n):
            output += '\n    '+ str((i+1) % n) +'  '+str(i)
        crossinfo = [] 
        crossid = 0 
        under = [] 
        above = []  
        for i in range(n):
            a0,b0,c0,a1,b1,c1 = x[i],y[i],z[i],x[(i+1)%n],y[(i+1)%n],z[(i+1)%n] 
            for j in range(i-1):

                x0,y0,z0,x1,y1,z1 = x[j],y[j],z[j],x[(j+1)%n],y[(j+1)%n],z[(j+1)%n]  
                A,B,C,D = a1-a0,-(x1-x0),b1-b0,-(y1-y0)
                det = A*D-B*C
                t = (D*(x0-a0)-B*(y0-b0))/det
                u = (-C*(x0-a0)+A*(y0-b0))/det	
                if 0<t and t<1 and 0<u and u<1:			
                    if c0+t*(c1-c0) > z0+u*(z1-z0):
                        aux = 1 
                        under += [j]
                        above += [i]
                    else: 	
                        aux = -1
                        under += [i]
                        above += [j]
                    coords = [x0+u*(x1-x0),y0+u*(y1-y0),z0+u*(z1-z0)] 
                    signum = math.copysign(1,(x1-x0)*(b1-b0)-(a1-a0)*(y1-y0))
                    crossinfo += [(j,i,aux,u,aux*signum,coords,crossid)] 
                    crossinfo += [(i,j,-aux,t,aux*signum,coords,crossid)] 
                    crossid += 1 		
        N = len(crossinfo) 
        totcross = N/2 	
        crossinfo.sort(key=operator.itemgetter(0,3)) 
        crossid = [row[6] for row in crossinfo] 

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
            with open(filename,'wt') as file_output:
                file_output.write(output)
	    self.filename = filename
            print('A file named "'+filename+'" was created in the current working directory.')
    
        crossings = [spherogram.links.links.Crossing(j) for j in range(totcross)] # to start a variable for crossing conventions
        k = 0
        while k < N:
            if crossinfo[k][2]==-1:
                prev = 0
                proc = 2
            else: 
                if crossinfo[k][4]==1:
                    prev = 3			
                    proc = 1
                else:
                    prev = 1
                    proc = 3
            if crossinfo[(k-1)%N][2]==-1: 
                crossings[crossid[k]][prev] = crossings[crossid[(k-1)%N]][2] 		
            else: 
                if crossinfo[(k-1)%N][4] == 1:
                    crossings[crossid[k]][prev] = crossings[crossid[(k-1)%N]][1] 
                else:
                    crossings[crossid[k]][prev] = crossings[crossid[(k-1)%N]][3]
            if crossinfo[(k+1)%N][2]==-1: 
                crossings[crossid[k]][proc] = crossings[crossid[(k+1)%N]][0] 		
            else: 
                if crossinfo[(k+1)%N][4] == 1: 
                    crossings[crossid[k]][proc] = crossings[crossid[(k+1)%N]][3] 
                else:
                    crossings[crossid[k]][proc] = crossings[crossid[(k+1)%N]][1]						
            k += 1	
        self.manifold = spherogram.links.links.Link(crossings).exterior()
