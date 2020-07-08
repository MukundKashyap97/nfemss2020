## File that generates analytical solution in linear-elastic 
# case of the given investigation.




from inputParams import *
from meshGenerator import a,b,rnodes,nelem
from numpy import array
def exact_sol(r):
    '''Function to generate the analytical solution of the given investigation
    Required parameters:
        r: radius as input
        UPSILON: poisson ration
        P_MAX: maximum applicable pressure
        a,b: radii
        E: young's modulus
    Output generated: an array called analytical_solution that consists of values.'''
    u = (1+UPSILON)*(P_MAX/E)*(a**2/(b**2-a**2))*(((1-(2*UPSILON))*r)+(b**2/r))
    return u
analytical_solution = []
r = 0 
for i in range(nelem+1):
    r = rnodes[i]
    analytical_solution.append(exact_sol(r)) 
analytical_solution = array(analytical_solution)