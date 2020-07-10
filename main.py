## NON LINEAR FINITE ELEMENTS 
## Assignment - Summer term 2020
# Author: Venkata Mukund Kashyap Yedunuthala
# Computational Materials Science
# This file consists of main solver, element routine, and material routine
# User is expected to execute this file AFTER updating inputParams.py


# importing required packages from python libraries
import numpy as np 
import time
import matplotlib.pyplot as plt
from datetime import datetime
# importing parameters and functions from other files
from inputParams import *
from meshGenerator import a,b,rnodes,nelem
from analytical_sol import analytical_solution
from convergenceStudy import convergenceStudy_elem,convergenceStudy_dt

def mat_routine(strain,del_strain,prev_sigma_ov,time_scale,Q,timeStep,C):
    '''Function for material routine employed within element subroutine
    Input parameters:
        strain          : strain values updated in the element routine
        del_strain      : differentiated strain values, also updated in the element routine
        prev_sigma_ov   : overstress values of the previous timestep
        time_scale      : characteristic time scale value, updated in inputParams.py
        Q               : modulus defined in inputParams.py
        timeStep        : increment in time in seconds
        C               : material matrix for a plain-strain element, updated in inputParams.py
    Output parameters:
        sigma           : updated stress values
        sigma_ov        : updated overstress values
        Ct              : material tangent matrix'''
    dt = timeStep
    T = time_scale
    # the deviator of strain rate is calculated
    dev = del_strain - (np.sum(del_strain)/3) 
    # which is used to update overstress using Euler Modified method
    sigma_ov = (1/(1+(dt/(2*T)))) * ((prev_sigma_ov*(1-(dt/(2*T)))) + Q*dev)
    # the material tangent matrix is thus calculated
    Ct = np.add(C,(Q/(1+(dt/(2*T))))*np.array([[2/3, -1/3],[-1/3, 2/3]]))
    # which are used to update stress
    sigma = np.matmul(C,strain) + sigma_ov
    return sigma, sigma_ov, Ct


def element_routine(rnodes,global_u,del_u,xi,T,Q,sigma_ov,timeStep,time_scale,C):
    ''' Function that describes the element subroutine that is employed for each element at every timestep
    Input parameters:
        rnodes      : the nodes at which element routine is performed
        global_u    : the global displacements corresponding to the element under consideration
        del_u       : the displacement increments from the previous iteration
        xi          : value of gauss point according to quadrature scheme
        T           : value of time at that particular timestep
        Q           : modulus value as defined in inputParams.py
        sigma_ov    : overstress values (of the previous timestep)
        timeStep    : time increment value
        time_scale  : characteristic time scale as defined in inputParams.py
        C           : material law for plain strain element
    Output parameters:
        elem_K      : element stiffness matrix
        elem_F      : nodal forces
        sigma       : updated stresses
        sigma_ov    : updated overstresses'''
    r1 = rnodes[0]
    r2 = rnodes[1]
    elem_u = global_u # assigning values to local variables for ease
    # the overstress value of previous timestep is stored to use in EM scheme
    prev_sigma_ov = sigma_ov
    # the jacobian for this case is calculated
    Jacobian = (r2-r1)/2
    # shape function N
    N = np.array(([0.5*(1-xi),0.5*(1+xi)]))
    # strain displacement matrix B
    B = np.array(([[-1/(2*Jacobian),1/(2*Jacobian)],[1/(r1+r2),1/(r1+r2)]]))
    B_t = np.transpose(B)
    # calculating strains
    strain = np.matmul(B,elem_u)
    del_strain = np.matmul(B,del_u)
    # stress update using material routine
    sigma, sigma_ov, Ct = mat_routine(strain,del_strain,prev_sigma_ov,time_scale,Q,timeStep,C)
    # calculating element stiffness matrix, and nodal forces
    elem_F = 2*np.matmul(B_t,sigma)*np.matmul(np.transpose(N),rnodes)*Jacobian
    elem_K = 2*np.matmul(B_t,np.matmul(Ct,B))*np.matmul(N,rnodes)*Jacobian
    # return statement for the values
    return elem_K, elem_F, sigma, sigma_ov

# main program begins here, with allocation of space for some required variables
sigma_ov = zeros((2,1))
sigma = zeros((2,1))
global_u = zeros((nelem+1,1))
placeholder_del_u = zeros((nelem+1,1))
placeholder_sigma_ov = zeros((2,1)) # allocating space and initializing to zero, required in first iteration
disp_evolution_last = [] # list to contain time history of displacements


## SOLVER
# =========
start_time = time.time() # making note of time at start to calculate total time taken by solver
## time loop
for T in TIME_SPAN:
    # load scaling
    if T <= LOADING_TIME:
        load_scale = (1/LOADING_TIME)*T
        P = P_MAX*a*load_scale
    else:
        P = P_MAX*a
    it = 0
    # newton-raphson loop
    while it <= MAX_ITERATIONS:
        # allocating space and initializing to zero at each iteration
        elem_K = np.zeros(2) # element stiffness matrix
        elem_F = np.zeros((2,1)) # nodal forces
        global_K = np.zeros((nelem+1,nelem+1)) # global stiffness matrix
        global_int_F = np.zeros((nelem+1,1)) # global internal forces
        global_ext_F = np.zeros((nelem+1,1)) # global external forces
        del_u = np.zeros((nelem+1,1)) # increment in displacement
        global_ext_F[0] = P # assigning value to global external forces
        sigma_evolution = [] # list to store values of stress at each time step
        # element routine loop and assembly of elemental matrices
        for k in range(nelem):
            elem_K, elem_F, sigma, sigma_ov = element_routine(rnodes[k:k+2],global_u[k:k+2],placeholder_del_u[k:k+2],GAUSS_POINT,T,Q,placeholder_sigma_ov,TIME_STEP,TIME_SCALE,C)
            sigma_evolution.append(sigma)
            global_K[k:k+2,k:k+2] = global_K[k:k+2,k:k+2] + elem_K
            global_int_F[k:k+2] = global_int_F[k:k+2] + elem_F
        # solving non-linear equation for displacement increment
        del_u = np.linalg.solve(global_K,(global_ext_F-global_int_F))
        # updating displacements
        global_u = global_u + del_u
        placeholder_del_u = del_u[:]
        # calculation of residual forces
        R = global_int_F - global_ext_F
        # newton-raphson convergence criteria
        if np.amax(np.abs(R)) < 0.005*(np.amax(np.abs(global_int_F))) and np.amax(np.abs(del_u)) < 0.005*np.amax(np.abs(global_u)):
            break
        else:
            it+=1
    disp_evolution_last.append(global_u[-1]) # for time history 
    placeholder_sigma_ov = sigma_ov # to carry forward to next time step
end_time =time.time() # time at the end of execution of solver
sigma_evolution = np.array(sigma_evolution) # converting list to a NumPy array for ease

# CONVERGENCE STUDY
# Uncomment the following lines to perform convergence study
# Generates plots that are saved directly to the current directory
# REMOVE ANY .JSON FILES PRESENT FROM PRIOR STUDY BEFORE EXECUTING
# =====
# for convergence study with respect to number of elements while keeping timestep constant, uncomment below line
# convergenceStudy_elem(nelem,rnodes,global_u)

# for convergence study with respect to time step keeping number of elements constant, uncomment below line
# convergenceStudy_dt(TIME_STEP,rnodes,global_u)
# =====

## Post-processing section
print("Convergence achieved in "+str(end_time-start_time)+" seconds")
# plot for displacement distribution
fig,ax1 = plt.subplots()
ax1.plot(rnodes,global_u,c='C0',marker='.',lw=0.7,label='Displacement')
ax1.set(
    xlabel = 'r (mm)',
    ylabel = 'Displacement '+r' $u_r$'+' (mm)',
    title = 'Distribution of '+r' $u_r$'+' at last time step'
 )
ax1.legend()
plt.figtext(0.01,0.01,"Plot generated on "+str(datetime.today()),fontsize='small')
fig.savefig("./displacementDistribution.png",dpi=600,metadata={'Author':"Venkata Mukund Kashyap Yedunuthala"})

# # plot for stress distributions at the last timestep
fig,ax2 = plt.subplots(nrows=2)
ax2[0].plot(rnodes[1:],sigma_evolution[:,0],c='C2',marker='.',label=r'$\sigma_{rr}$')
ax2[1].plot(rnodes[1:],sigma_evolution[:,1],c='C3',marker='.',label=r'$\sigma_{\phi\phi}$')
ax2[0].set(
    title = 'Distribution of '+r'$\sigma_{rr}$ '+'and '+r'$\sigma_{\phi\phi}$ '+'at last time step',
    ylabel = 'Stress '+r'$\sigma_{rr}$' +' (MPa)'
)
ax2[1].set(
    xlabel = 'r'+' (mm)',
    ylabel = 'Stress '+r'$\sigma_{\phi\phi}$'+' (MPa)'
)
ax2[0].legend()
ax2[1].legend()
plt.figtext(0.01,0.01,"Plot generated at "+str(datetime.today()),fontsize='small')
fig.savefig("./stressDistribution.png",dpi=600,metadata={'Author':"Venkata Mukund Kashyap Yedunuthala"})


# # plot for comparison of displacements obtained with analytical, for verification purposes
fig,ax3 = plt.subplots()
ax3.plot(rnodes,global_u,c='C0',marker='.',lw=1,label='Obtained solution')
ax3.plot(rnodes,analytical_solution,'--',c='C3',marker='*',lw=0.5,label='Analytical solution')
ax3.set(
    xlabel = 'r'+' (mm)',
    ylabel = 'Displacement '+r' $u_r$'+' (mm)',
    title = 'Comparison of computed and analytical '+r' $u_r$'+' at last time step'
 )
ax3.legend()
plt.figtext(0.01,0.01,"Plot generated on "+str(datetime.today()),fontsize='small')
fig.savefig("./comparison.png",dpi=600,metadata={'Author':"Venkata Mukund Kashyap Yedunuthala"})


# # plot for the time history of displacement 
fig,ax4 = plt.subplots()
ax4.plot(TIME_SPAN,disp_evolution_last,label='Displacement')
ax4.set(
    xlabel = 'Time t (s)',
    ylabel = 'Displacement '+r' $u_r$' +' (mm)',
    title = 'Evolution of displacement '+r' $u_r$ (r=b)'+' over time'
 )
ax4.legend(loc=(0.6,0.07))
plt.figtext(0.01,0.01,"Plot generated on "+str(datetime.today()),fontsize='small')
fig.savefig("./displacement_evolution.png",dpi=600,metadata={'Author':"Venkata Mukund Kashyap Yedunuthala"})