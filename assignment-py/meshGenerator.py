#   Generate list of position of nodes according to a geometric series
#   for assignement in "Nonlinear Finite Element Methods" 
#   in summer term 2020
#   lecturer in charge: Dr. Geralf H�tter
from inputParams import INNER_RADIUS,OUTER_RADIUS, NUMBER_OF_ELEMENTS
from matplotlib.pyplot import plot,subplots,savefig,figtext,show
from numpy import zeros, array,transpose
from datetime import datetime
## Input parameters
b = OUTER_RADIUS # Outer radius
a = INNER_RADIUS # Inner radius
nelem = NUMBER_OF_ELEMENTS # number of elements
meshrefinementfactor = 2 # ratio of element sizes at outer and inner radius
# ratio between element sizes of subsequent elements for a geometric series
q = meshrefinementfactor**(1.0/(nelem-1))
# size of first interval
dr = (b-a)*(1-q)/(1-meshrefinementfactor*q)
rnode = a
rnodes = [a,]
# loop over all elements
for i in range(nelem):
    rnode = rnode + dr
    rnodes.append(rnode)
    dr = dr*q
# visulaizing location of nodes
rnodes = array(rnodes)
fig,ax = subplots()
ax.plot(rnodes,zeros(nelem+1),'x',label="nodes")
figtext(0.01,0.01,"Mesh generated at "+str(datetime.today()),fontsize='small')
ax.set(xlabel='r (mm)')
ax.legend()
savefig("generatedMesh.png",dpi=600,metadata={'Author':"Venkata Mukund Kashyap Yedunuthala"})