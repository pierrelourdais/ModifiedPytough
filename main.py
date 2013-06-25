#!/usr/bin/python

'''
Created on 22 mai 2013
@author: Pierre
---------------------------------
Purpose: 
enables all the TOUGH2 simulation
---------------------------------
Date:   \ Comments:
06/12/13\3D Mesh
        \
---------------------------------

---------------------------------
Date:   \ Modifications:
        \
        \
---------------------------------
'''
#Packages
import math
import numpy
from mulgrids import *
from t2data import *
from t2grids import *
from t2incons import *
from t2listing import *
from t2thermo import *
from copy import copy
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import vtk

'''
SECONDARY FUNCTIONS 
* initial_interp: set, for each block, the initial temperature thanks to the thermocouples (wall + vessel)
    datas.
ERROR:
    - need scipy
    - interp2 need to be tested

* setinitialT: extract datas from experiments and enable to fill thermocouples[i][2]
ERROR: not tested

* setinitialTw: the same, withe the wall thermocouples
ERROR: not tested

'''
def initial_interp(position,thermocouples):
#    Ecrire une fonction qui prend en argument thermocouples[i][2]
#    Qui effectue une interpolation logarithmique radiale
#    Qui effectue une interpolation verticale lineaire
#    Qui retourne la temperarure interpolee
    positionw=numpy.zeros(number_thermocouplesw+2)
    positionw[0]=length_bottom
    positionw[-1]=length_bottom+length_vessel
    temperaturew=numpy.zeros(number_thermocouplesw+2)
    temperaturew[0]=T_bottom
    temperaturew[-1]=T_top
    for i in range(1,number_thermocouplesw):
        positionw[i]=thermocouplesw[i][0]
        temperaturew[i]=thermocouplesw[i][1]
    
    ri=(radii_sand*0,radii_sand*1./3,radii_sand*1./2,radii_sand*2./3,radii_sand*1)
    zi=(length_bottom+length_vessel*(0/5),length_bottom+length_vessel*(1./5),length_bottom+length_vessel*(2./5),length_bottom+length_vessel*(3./5),length_bottom+length_vessel*(4./5),length_bottom+length_vessel*(5./5))    
#    Assumption: only conduction, no convection
#     Thermal Diffusion: T(r)= A ln (r) + B
#     It makes sense between the two closest points
    a1=(thermocouples[3][2]-thermocouples[2][2])/math.log(2)*math.log(2./3)+thermocouples[2][2]
    a2=(thermocouples[3][2]-thermocouples[2][2])/math.log(2)*math.log(4./3)+thermocouples[2][2]
    a3=thermocouples[5][2]+(thermocouples[6][2]-thermocouples[5][2])/math.log(2)*math.log(3./2)
    a4=thermocouples[9][2]+(thermocouples[10][2]-thermocouples[9][2])/math.log(2)*math.log(3./2)
    a5=thermocouples[13][2]+(thermocouples[14][2]-thermocouples[13][2])/math.log(2)*math.log(3./2)
    a6=thermocouples[18][2]+(thermocouples[19][2]-thermocouples[18][2])/math.log(2)*math.log(3./2)
    a7=(thermocouples[18][2]-thermocouples[17][2])/math.log(2)*math.log(2./3)+thermocouples[17][2]
    a8=0
    a9=0
    Ti=[[thermocouples[1][2],a1,thermocouples[2][2],a2,thermocouples[3][2]],
        [thermocouples[4][2],thermocouples[5][2],a3,thermocouples[6][2],thermocouples[7][2]],
        [thermocouples[8][2],thermocouples[9][2],a4,thermocouples[10][2],thermocouples[11][2]],
        [thermocouples[12][2],thermocouples[13][2],a5,thermocouples[14][2],thermocouples[15][2]],
        [thermocouples[17][2],thermocouples[18][2],a6,thermocouples[19][2],thermocouples[20][2]],
        [thermocouples[17][2],a7,thermocouples[18][2],a8,a9]
        ]
    import scipy.interpolate 
#     rii=range(0,radii_sand,radialelements)
#     zii=range(length_bottom,length_bottom+length_vessel,verticalelements)
    newfunc=scipy.interpolate.interp2d(ri,zi,Ti,kind='cubic')
    
    if length_bottom<=position[2]<length_bottom+length_vessel :
        if 0<=position[0]<radii_sand: interpolation=newfunc(position[0],position[2])
        else: interpolation=numpy.interp(position[2], positionw, temperaturew)
    elif 0<= position[2] < length_bottom : interpolation= T_bottom
    else: interpolation = T_top
    
    return interpolation



def initial_interpw(position,thermocouplesw):
    z=position[2]
    import numpy.interp
    zi=[]
    Ti=[]
    for i in range(number_thermocouplesw):
        zi.append(thermocouplesw[i][0])
        Ti.append(thermocouplesw[i][1])
    interpolation=numpy.interp(z,zi,Ti)
    return interpolation
    

    
def exportnumericaldatas(infile, outfile):
    '''
Enable to get the results.
resultat['time']= number of time steps
resultat['elements']=number of elements
resultat[i][j]['options']=value of options for time i and elements j
--options:
    T
    P
    name
    DL
    DG
    xco2aq

'''
    outfile=outfile
    f=open(infile, 'r')
    content=f.readlines()
    separateur=[]
    time=[]
    
    str=content[0][:9]
    compteur=0
    while str!=' MESH HAS': 
        compteur+=1
        str=content[compteur][:9]    
    elements=int(content[compteur].split()[2])
    
    
    
    
    
    for i in range(len(content)):
        if content[i][:27]=='          OUTPUT DATA AFTER':
            separateur.append(i)
            time.append(float(content[i+4][1:13]))
    separateur.append(len(content))
    resultat={}
    resultat['time']=len(separateur)-1
    resultat['elements']=elements
    print resultat['time']
    
    import re
    for i in range(len(separateur)-1):
        resultat[i]={}
        resultat[i]['time']=time[i]
        element=0
        for j in range(separateur[i],separateur[i+1]):
            line=content[j]
            matchObj=re.search(r' (\S)(\S)(\S) (\S) ..(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)',line)
            if matchObj!=None:
                element+=1
                a=line.split()
                name=a[0]+a[1]
                resultat[i][element]={}
                resultat[i][element]['name']=name
                resultat[i][element]['T']=float(a[4])
                resultat[i][element]['P']=float(a[3])
                resultat[i][element]['DL']=float(a[12])
                resultat[i][element]['DG']=float(a[11])
                resultat[i][element]['xco2aq']=float(a[8])
    import csv
    w = csv.writer(open(outfile, "w"))
    for key, val in resultat.items():
        w.writerow([key, val])
    
    return resultat
            
    
#    for i in range(len(separateur)):
#         difference=0
#         while (content[separateur[i]+difference]+'None').split()[0]!='0b0':
#             difference+=1
        
        
#         resultat[i]={}
#         resultat[i]['time']=time[i]
#         for j in range(elements):
#             blank=0
#             test=content[separateur[i]+difference+j]
#             
#             test1=test.split()
#             
#             while (len(test1)==0 or test1[-1]!='0.00'):
#                 blank+=1
#                 test=content[separateur[i]+difference+j+blank]
#                 test1=test.split()
#                 
#     
#                 
#                 
#             resultat[i][j]={}
#             a=content[separateur[i]+j+blank+difference].split()
#             resultat[i][j]['T']=float(a[4])
#             resultat[i][j]['P']=float(a[3])    
#             resultat[i][j]['name']=a[0]+a[1] 
#             resultat[i][j]['DL']=float(a[12]) 
#             resultat[i][j]['DG']=float(a[11])
#             resultat[i][j]['xco2aq']=float(a[8])
#             
#         import csv
#         w = csv.writer(open(outfile, "w"))
#         for key, val in resultat.items():
#             w.writerow([key, val])
#     
#     return resultat

def importnumericaldatas(infile):

    infile=infile
    import csv
    dict = {}
    for key, val in csv.reader(open(infile)):
        dict[key] = val
    
    dic={}
    
    
    a=dict['0']
    b=a.split('}, ')
    dic['time']=int(dict['time'])
    dic['elements']=int(dict['elements'])
    for i in range(len(dict)-2):
        dic[i]={}
        a=dict[str(i)].split('}, ')
        c=a[len(a)-1].split()
        dic[i]['time']=float(c[1][0:-1])
        for j in range(len(a)-2):
            dic[i][j]={}
            temp=a[j].split()
            dic[i][j]['T']=float(temp[12])
            dic[i][j]['P']=float(temp[10][0:-1])
            dic[i][j]['DL']=float(temp[2][0:-1])
            dic[i][j]['name']=temp[4][1:-1]
            dic[i][j]['DG']=float(temp[6][0:-1])
            dic[i][j]['xco2aq']=float(temp[8][0:-1])
    
    return dic

def findindiceblock(position, gridblock):
    #find the nearest block
    #return indice in blocklist
    distance_min_square=length_vessel**2
    indice=None
    for i in range(len(gridblock.blocklist)):
        distance_square=0
        for j in range(3):
            distance_square+=(position[j]-gridblock.blocklist[i].centre[j])**2
        if distance_square<distance_min_square:
            indice=i
            distance_min_square=distance_square
    return indice
        
        

'''INPUT DATA'''
#Simulation name
name='vessel'
filename='vessel.core'
output='vessel.listing'
outfile='output.csv'
import_file_path=''
workfile='6.18.2013_11.07 AM.txt'


#Dimensions: (meters)
length_top= 0.1
length_bottom= 0.1
length_vessel= 0.508

radii_sand=0.08
radii_vessel=0.01
radii_pipe=0.01

#Thermocouple datas, create a dictionnary
# 'number':[radial,vertical_position_in_the_vessel, initial_temperature]
T_bottom=60
T_top=80
number_thermocouples=22
T=setinitialT()
thermocouples={1:[0*        radii_sand,0.0,   0*   length_vessel/5 +length_bottom,    T[0]],
               2:[1./2*  radii_sand,0.0,   0*   length_vessel/5 +length_bottom,    T[1]],
               3:[1*    radii_sand,0.0,   0*   length_vessel/5 +length_bottom,    T[2]],
               4:[0*    radii_sand,0.0,   1*   length_vessel/5 +length_bottom,    T[3]],
               5:[1./3*  radii_sand,0.0,   1*   length_vessel/5 +length_bottom,    T[4]],
               6:[2./3*  radii_sand,0.0,   1*   length_vessel/5 +length_bottom,    T[5]],
               7:[1*    radii_sand,0.0,   1*   length_vessel/5 +length_bottom,    T[6]],
               8:[0*    radii_sand,0.0,   2*   length_vessel/5 +length_bottom,    T[7]],
               9:[1./3*  radii_sand,0.0,   2*   length_vessel/5 +length_bottom,    T[8]],
               10:[2./3* radii_sand,0.0,   2*   length_vessel/5 +length_bottom,    T[9]],
               11:[1*   radii_sand,0.0,   2*   length_vessel/5 +length_bottom,    T[10]],
               12:[0*   radii_sand,0.0,   3*   length_vessel/5 +length_bottom,    T[11]],
               13:[1./3* radii_sand,0.0,   3*   length_vessel/5 +length_bottom,    T[12]],
               14:[2./3* radii_sand,0.0,   3*   length_vessel/5 +length_bottom,    T[13]],
               15:[1*   radii_sand,0.0,   3*   length_vessel/5 +length_bottom,    T[14]],
               16:[-2./3*radii_sand,0.0,   3*   length_vessel/5 +length_bottom,    T[15]],
               17:[0*   radii_sand,0.0,   4*   length_vessel/5 +length_bottom,    T[16]],
               18:[1./3* radii_sand,0.0,   4*   length_vessel/5 +length_bottom,    T[17]],
               19:[2./3* radii_sand,0.0,   4*   length_vessel/5 +length_bottom,    T[18]],
               20:[1*   radii_sand,0.0,   4*   length_vessel/5 +length_bottom,    T[19]],
               21:[0*   radii_sand,0.0,   5*   length_vessel/5 +length_bottom,    T[20]],
               22:[1./2* radii_sand,0.0,   5*   length_vessel/5 +length_bottom,    T[21]]
               }

number_thermocouplesw=9;
Tw=setinitialTw()
thermocouplesw={1:[length_bottom+0          ,Tw[0]],
                2:[length_bottom+0.053      ,Tw[1]],
                3:[length_bottom+0.097      ,Tw[2]],
                4:[length_bottom+0.147      ,Tw[3]],
                5:[length_bottom+0.198      ,Tw[4]],
                6:[length_bottom+0.296      ,Tw[5]],
                7:[length_bottom+0.395      ,Tw[6]],
                8:[length_bottom+0.505      ,Tw[7]],
                9:[length_bottom+0.508      ,Tw[8]]
                }
'''MESHFILE'''
#Grid realisation
#######################################################################################
#Geometry
hpipei=0.01
hpipe_=0.1
hsand=0.507
hpipeo=hpipei
hbottom=0.03
hbottom3=hpipei+hpipe_-hbottom
htop=hbottom
htop3=hbottom3

rpipe=0.01
rsand=0.045
rwall=0.05
rcap=0.07

#######################################################################################
#Mesh
theta_step=10
layer_step_middle=10
layer_step_top=5
layer_step_bottom=5
middle_radii_step=10
top_radii_step=int(middle_radii_step/rwall*rcap)
bottom_radii_step=int(middle_radii_step/rwall*rcap)


#######################################################################################
vessel=t2data()
#######################################################################################
vessel_grid=t2grid()
vessel_grid.add_rocktype(rocktype('SAND1',nad=2,permeability=[9.3e-12]*3,density=2600.0,porosity=0.3942,conductivity=0.5,specific_heat=830,compressibility=4.5e-10))
vessel_grid.rocktypelist[0].relative_permeability={'type':7, 'parameters':[.457,.3,1,.05]}
vessel_grid.rocktypelist[0].capillarity={'type':8, 'parameters':[.457,0,5.1e-5,1e7,.999]}

vessel_grid.add_rocktype(rocktype('STEEL',permeability=[0.0]*3,density=8000.0,porosity=1.0e-9,conductivity=5.3,specific_heat=300))
vessel_grid.add_rocktype(rocktype('PIPE_',permeability=[1.0e-9]*3,density=8000.0,porosity=0.99,conductivity=6.689e-02,specific_heat=500))
vessel_grid.add_rocktype(rocktype('PIPEI',permeability=[1.0e-9]*3,density=2600.0e40,porosity=0.11,conductivity=16.3e9,specific_heat=1.0e4))
vessel_grid.add_rocktype(rocktype('PIPEO',permeability=[1.0e-9]*3,density=8000.0,porosity=0.99,conductivity=6.689e-02,specific_heat=500))

#######################################################################################


############################
#MIDDLE
middle=mulgrid()
middle.create_cylinder(hmin=hbottom+hbottom3,length=hsand,radii=rwall,theta_step=theta_step,layer_step=layer_step_middle,radii_step=middle_radii_step,nameofcylinder='m',filename='middle')
middle_grid=t2grid()
middle_grid.fromgeo(middle)
middle_grid.get_block_centres_defined()
for blk in middle_grid.blocklist[1:]:
    radii_blk=math.sqrt(blk.centre[0]**2+blk.centre[1]**2)
    if radii_blk<rsand:blk.rocktype=vessel_grid.rocktype['SAND1']
    else:blk.rocktype=vessel_grid.rocktype['STEEL']
############################

############################
#top
top=mulgrid()
top.create_cylinder(hmin=hbottom+hbottom3+hsand,length=htop3+htop,radii=rcap,theta_step=theta_step,layer_step=layer_step_top,radii_step=top_radii_step,nameofcylinder='t',filename='top')
top_grid=t2grid()
top_grid.add_rocktype(rocktype('PIPE_',permeability=[1.0e-9]*3,density=8000.0,porosity=0.99,conductivity=6.689e-02,specific_heat=500))
top_grid.fromgeo(top)
top_grid.get_block_centres_defined()
for blk in top_grid.blocklist[1:]:
    radii_blk=math.sqrt(blk.centre[0]**2+blk.centre[1]**2)
    if blk.centre[2]<hbottom+hbottom3+hsand+hpipe_:
        if radii_blk<rpipe:blk.rocktype=vessel_grid.rocktype['PIPE_']
        else:blk.rocktype=vessel_grid.rocktype['STEEL']
    else:
        if radii_blk<rpipe:blk.rocktype=vessel_grid.rocktype['PIPEO']
        else:blk.rocktype=vessel_grid.rocktype['STEEL']
############################

############################
#bottom
bottom=mulgrid()
bottom.create_cylinder(hmin=0,length=hbottom3+hbottom,radii=rcap,theta_step=theta_step,layer_step=layer_step_bottom,radii_step=bottom_radii_step,nameofcylinder='b',filename='bottom')
bottom_grid=t2grid()
bottom_grid.fromgeo(bottom)
bottom_grid.get_block_centres_defined()
for blk in bottom_grid.blocklist[1:]:
    radii_blk=math.sqrt(blk.centre[0]**2+blk.centre[1]**2)
    if blk.centre[2]>hpipei:
        if radii_blk<rpipe:blk.rocktype=vessel_grid.rocktype['PIPE_']
        else:blk.rocktype=vessel_grid.rocktype['STEEL']
    else:
        if radii_blk<rpipe:blk.rocktype=vessel_grid.rocktype['PIPEI']
        else:blk.rocktype=vessel_grid.rocktype['STEEL']
############################

############################
#MESH GENERAL
vessel_grid=bottom_grid+top_grid#+middle_grid
vessel_grid.get_block_centres_defined()
############################
#CONNECTION
vessel_grid.create_vertical_connection(middle_grid, rwall/middle_radii_step, bottom_grid, rcap/bottom_radii_step, theta_step)
vessel_grid.create_vertical_connection(top_grid, rcap/top_radii_step, middle_grid, rwall/middle_radii_step, theta_step)

###########################
for blk in vessel_grid.blocklist[:]:
    blk.ahtx=0
#


'''INITIAL CONDITIONS'''
incon=t2incon()
vessel_grid.blocklist[0].centre=[0,0,length_top]
for blk in vessel_grid.blocklist[1:]:
    if blk.rocktype=='SAND1':
        T=initial_interp(blk.centre,thermocouples)
        t=t2blockincon(block=blk.name, porosity=0,variable=[.1379e8,0,1,T])
        incon.add_incon(t)
    
    elif blk.centre[2]<length_bottom:
        T=thermocouplesw[1][1]
        t=t2blockincon(block=blk.name, porosity=0,variable=[.1379e8,0,1,T])
        incon.add_incon(t)
    elif blk.centre[2]>length_bottom +length_vessel:  
        T=thermocouplesw[number_thermocouplesw][1]
        t=t2blockincon(block=blk.name, porosity=0,variable=[.1379e8,0,1,T])
        incon.add_incon(t) 
    else:    
        T=initial_interpw(blk.centre,thermocouplesw)
        t=t2blockincon(block=blk.name, porosity=0,variable=[.1379e8,0,1,T])
        incon.add_incon(t)
        
    
'''TIME'''
#set times
times=[]
    
    
    
    
'''GENERATION '''
#set GENER
rate_inj=[1.0e-6,1.0e-4,1.6e-3]
time_inj=[0,0.01,0.02]
for i in range(1,200):
    time_inj.append(i)
    rate_inj.append(1.6e-3)


injection=t2generator(name='inj 1',block='at 0',type='COM3',gx=None,itab='',ltab=4,rate=rate_inj,time=time_inj,ex=None,fg=None,hg=None)




'''GENERAL PARAMETERS'''
#set MULTI
multi={'num_components':3, 'num_equations':4, 'num_phases':3, 'num_secondary_parameters':6}

#set SELEC
selec={'integer':[1,None,None,None,None, None,None,None,None,0,0,0,0,0,0,0], 'float':[.8,.8,None,None,None,None,None,None,None,None,None,None,None,None,]}

#set SOLVR
solver={'closure':1e-7, 'relative_max_iterations':.8, 'type':4, 'o_precond':'O0','z_precond':'Z1'}

#set COFT
hconnection=[]

#set FOFT
hblock=[]

#set INDOM
indom={'BSINK':[1.3790e07,0,1,10], 'PIPEI':[1.3790e07,0,1,15]}
#Set PARAM
parameters={'max_iterations':None, 'print_level':2, 'max_timesteps':2000, 'max_duration':None, 'print_interval':1, 
                    '_option_str':'','option':(1,1,0,0,0,0,0,0,0,0,0,2,0,0,0,0,4,0,0,0,0,4,0,0,0), 'diff0':None, 'texp':None, 'tstart':0.0, 'tstop':20000,
                    'const_timestep':-1.0,'timestep':[1], 'max_timestep':None, 'print_block':None, 'gravity':9.81,
                    'timestep_reduction':None, 'scale':None, 'relative_error':1.0e-5, 'absolute_error':None, 'pivot':None,
                    'upstream_weight':None, 'newton_weight':None, 'derivative_increment':1e-8, 'default_incons':[1.3790e07,0,1,95]}

'''.CORE CREATION'''
vessel=t2data()

vessel.title=name
vessel.filename=filename
vessel.grid=vessel_grid

vessel.grid.delete_rocktype('dfalt')
vessel.grid.add_rocktype(rocktype('SAND1',nad=2,permeability=[9.3e-12]*3,density=2600.0,porosity=0.3942,conductivity=0.5,specific_heat=830,compressibility=4.5e-10))
vessel.grid.rocktypelist[0].relative_permeability={'type':7, 'parameters':[.457,.3,1,.05]}
vessel.grid.rocktypelist[0].capillarity={'type':8, 'parameters':[.457,0,5.1e-5,1e7,.999]}
vessel.grid.add_rocktype(rocktype('STEEL',permeability=[0.0]*3,density=8000.0,porosity=1.0e-9,conductivity=5.3,specific_heat=300))
vessel.grid.add_rocktype(rocktype('PIPE_',permeability=[1.0e-9]*3,density=8000.0,porosity=0.99,conductivity=6.689e-02,specific_heat=500))
vessel.grid.add_rocktype(rocktype('PIPEI',permeability=[1.0e-9]*3,density=2600.0e40,porosity=0.11,conductivity=16.3e9,specific_heat=1.0e4))
vessel.grid.add_rocktype(rocktype('PIPEO',permeability=[1.0e-9]*3,density=8000.0,porosity=0.99,conductivity=6.689e-02,specific_heat=500))
vessel.grid.delete_block('ATM00')


vessel.generator=injection
vessel.generatorlist=[injection]
vessel.incon=incon
vessel.multi=multi
vessel.selection=selec
vessel.start=True 
vessel.solver=solver
vessel.history_connection=hconnection
vessel.history_block=hblock
vessel.indom=indom
vessel.parameter=copy(parameters)
vessel.output_times=times
vessel.indom=indom

vessel.write()

'''VTK'''
print top.num_atmosphere_blocks
print top_grid.rocktypelist
top_grid.add_rocktype(rocktype('SAND1',nad=2,permeability=[9.3e-12]*3,density=2600.0,porosity=0.3942,conductivity=0.5,specific_heat=830,compressibility=4.5e-10))
top_grid.rocktypelist[0].relative_permeability={'type':7, 'parameters':[.457,.3,1,.05]}
top_grid.rocktypelist[0].capillarity={'type':8, 'parameters':[.457,0,5.1e-5,1e7,.999]}

top_grid.add_rocktype(rocktype('STEEL',permeability=[0.0]*3,density=8000.0,porosity=1.0e-9,conductivity=5.3,specific_heat=300))
top_grid.add_rocktype(rocktype('PIPE_',permeability=[1.0e-9]*3,density=8000.0,porosity=0.99,conductivity=6.689e-02,specific_heat=500))
top_grid.add_rocktype(rocktype('PIPEI',permeability=[1.0e-9]*3,density=2600.0e40,porosity=0.11,conductivity=16.3e9,specific_heat=1.0e4))
top_grid.add_rocktype(rocktype('PIPEO',permeability=[1.0e-9]*3,density=8000.0,porosity=0.99,conductivity=6.689e-02,specific_heat=500))


print top_grid.rocktypelist

top_grid.write_vtk(geo=top, filename='vtk', wells=False)





'''RUN SIMULATION OR IMPORT SIMULATION'''
#If import simualtion, run=false

run=True
if run:vessel.run(incon_filename='vessel.core',simulator='zco2n    ', silent=False,save_filename='out')



'''POST-PROCESSING'''
if run:resultat=exportnumericaldatas(infile='vessel.listing',outfile=outfile)
else: resultat=importnumericaldatas(infile=outfile)

'''PLOTTING'''
Texperiment=setexperimentmesure()
color_level=['g','r','c','b','y','m']
stage=[0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5]
for i in range(1,number_thermocouples+1):
    
    x1=[]
    y1=[]
    y2=[]
    positionthermocouples=[]
    for k in range(3):
        positionthermocouples.append(thermocouples[i][k])
    indice=findindiceblock(position=positionthermocouples,gridblock=vessel.grid)    
    
    for j in range(resultat['time']):
        
        y1.append(resultat[j][indice]['T'])
        x1.append(resultat[j]['time'])
    plt.plot(x1,y1,'ro',color=color_level[stage[i]],linestyle="dashed", marker="o",)
    x2=range(len(Texperiment)-3)
    for h in range(len(Texperiment)-3):
        y2.append(Texperiment[h][i-1])
    plt.plot(x2,y2,'ro',color=color_level[stage[i]],linestyle="solid", marker="x",)
    xlabel='Time (s)'
    ylabel= 'Temperature (C)'
    name='tc'+str(i)+'.png'
    title='Thermocouple '+str(i)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel=(ylabel)
    plt.savefig(name,facecolor='w', edgecolor='w')
    plt.close()


















