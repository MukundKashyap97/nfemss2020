import json
from json.decoder import JSONDecodeError
import matplotlib.pyplot as plt
from datetime import datetime
from numpy import transpose,array

def convergenceStudy_elem(nelem,rnodes,global_u):
    """Function to generate json files required for convergence study
    with respect to number of elements 
    This function is called in main file. 
    Input parameters:
        nelem       : number of elements
        rnodes      : coordinates of nodes
        global_u    : global displacements at nodes
    Output parameters:
    Returns None, generates json file consisting of displacements
    depending upon the number of elements. """

    study_u_nelem = {}
    study_u_nelem[nelem] = [rnodes.tolist(),transpose(global_u).tolist()]
    try:
        with open("convergenceStudy_nelem.json", 'r') as infile:
            old_data = json.load(infile)
            data = old_data + [study_u_nelem,]
        with open("convergenceStudy_nelem.json", 'w') as outfile:
            json.dump(data,outfile,sort_keys=True,indent=4)
    except:
        with open("convergenceStudy_nelem.json",'w') as outfile:
            json.dump([study_u_nelem,],outfile,indent=4)
    
def convergenceStudy_dt(TIME_STEP,rnodes,global_u):
    """Function to generate json files required for convergence study
    with variation in time step. 
    This function is called in main file. 
    Input parameters:
        TIME_STEP   : increment in time
        rnodes      : coordinates of nodes
        global_u    : global displacements at nodes
    Output parameters:
    Returns None, generates json file consisting of displacements
    depending upon the increment in time. """
    study_u_dt = {}
    study_u_dt[TIME_STEP] = [rnodes.tolist(),transpose(global_u).tolist()]
    try:
        with open("convergenceStudy_dt.json", 'r') as infile:
            old_data = json.load(infile)
            data = old_data + [study_u_dt,]
        with open("convergenceStudy_dt.json", 'w') as outfile:
            json.dump(data,outfile,sort_keys=True,indent=4)
    except:
        with open("convergenceStudy_dt.json",'w') as outfile:
            json.dump([study_u_dt,],outfile,indent=4)

def convergenceStudy_elem_publish():
    """ Function to generate plot of convergence study
    with respect to number of elements. 
    REQUIRES JSON FILES TO BE PRESENT AS GENERATED USING ABOVE FUNCTIONS
    ---
    Execution in this file, see below.
    Generates convergenceStudy_u_nelem.png file """
    with open("convergenceStudy_nelem.json","r") as f:
        data = json.load(f)
    fig,ax1 = plt.subplots()
    for d in data:
        for values in d.values():
            ax1.plot(values[0],values[1][0],lw=0.7,label=str(len(values[0])-1)+" elements")
    ax1.set(
        xlabel = 'r'+' (mm)',
        ylabel = 'Displacement '+r' $u_r$'+' (mm)',
        title = 'Convergence study of '+r' $u_r$'+' with respect to number of elements'
    )
    ax1.legend()
    plt.figtext(0.01,0.01,"Plot generated on "+str(datetime.today()),fontsize='small')
    fig.savefig("./convergenceStudy_u_nelem.png",dpi=600,metadata={'Author':"Venkata Mukund Kashyap Yedunuthala"})

def convergenceStudy_dt_publish():
    """ Function to generate plot of convergence study
    with respect to increment in time. 
    REQUIRES JSON FILES TO BE PRESENT AS GENERATED USING ABOVE FUNCTIONS
    ---
    Execution in this file, see below.
    Generates convergenceStudy_u_dt.png file """
    with open("convergenceStudy_dt.json","r") as f:
        data = json.load(f)
    fig,ax2 = plt.subplots()
    for d in data:
        for key, values in zip(d.keys(),d.values()):
            ax2.plot(values[0],values[1][0],lw=0.7,label='Time step = '+str(key))
    ax2.set(
        xlabel = 'r (mm)',
        ylabel = 'Displacement '+r' $u_r$'+' (mm)',
        title = 'Convergence study of '+r' $u_r$'+' with respect to number of elements'
    )
    ax2.legend()
    plt.figtext(0.01,0.01,"Plot generated on "+str(datetime.today()),fontsize='small')
    fig.savefig("./convergenceStudy_u_dt.png",dpi=600,metadata={'Author':"Venkata Mukund Kashyap Yedunuthala"})

## Uncomment the following to generate plots
# ====    
# convergenceStudy_dt_publish()
# convergenceStudy_elem_publish()
# ====