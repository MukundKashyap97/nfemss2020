## File consisting of all input parameters required
# part of the assignment of non linear finite element methods
# summer term 2020
# Author: Venkata Mukund Kashyap Yedunuthala


# To import required packages from numpy library
from numpy import array,zeros,arange

# Initializing all required parameters
INNER_RADIUS = 60 # internal radius in mm
OUTER_RADIUS = 120 # external radius in mm
E = 7e4 # Young's modulus in MPa
UPSILON = 0.30 # poisson ration
Q = 35000 # modulus Q in MPa
TIME_SCALE = 3 # characteristic time scale in s
P_MAX = 50 # maximum applicable pressure in MPa
START_TIME = 0 # time at which solver starts
LOADING_TIME = 6 # time till loading stops
END_TIME = 30 # time at which solver stops 
TIME_STEP = 0.01 # incrememnt in time
MAX_ITERATIONS = 20 # maximum number of iterations in NR scheme
GAUSS_POINT = 0 # gauss point as per quadrature scheme used
NUMBER_OF_ELEMENTS = 10 # number of elements in the mesh

# the following generates the material matrix C
C = (E/((1+UPSILON)*(1-2*UPSILON)))*array([[(1-UPSILON),UPSILON],[UPSILON,(1-UPSILON)]])

# the following generates the array of the times
TIME_SPAN = arange(START_TIME,END_TIME,TIME_STEP)