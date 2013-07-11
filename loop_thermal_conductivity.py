#!/usr/bin/python

'''
Created on 18 juin 2013

@author: Pierre
'''
###################################
#MODULE
#Math: sin,cos,sqrt..
#Numpy: interp1D
#Mulgrids,t2data,t2grids,t2incons,t2listing,t2thermo: PyTOUGH (Modified)
#Copy: copy csv file
#Matplotlib: draw figures
#Scipy.interpolate: 2D interpolation
#os: crate directory, run tough2...
import math
import numpy
from mulgrids_l import *
from t2data_l import *
from t2grids_l import *
from t2incons_l import *
from t2listing_l import *
from t2thermo_l import *
from copy import copy
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
#import vtk
import scipy.interpolate 
import os


#Grid realisation
#######################################################################################
#Geometry
#You can see geometry: vessel.jpg
hpipei=0.05
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
theta_step=6
layer_step_middle=10
layer_step_top=6
layer_step_bottom=6
middle_radii_step=10
top_radii_step=int(middle_radii_step/rwall*rcap)
bottom_radii_step=int(middle_radii_step/rwall*rcap)


#######################################################################################
#Initialisation Thermocouples
#It is possible to add other thermocouples.

#VOLUME ONES
number_thermocouples=22
#i=1...number_thermocouples: name of the thermocouple. It is an integer (not a string).
#thermocouples[i][0]=radial position
#thermocouples[i][1]=theta position     Useless because everything is in a plane, but I needed it because of the 3D Mesh
#thermocouples[i][2]=vertical position
#thermocouples[i][3]=temperature
#thermocouples[i][4]=associated block
thermocouples={1:[0*     rsand,0.0,   0*   hsand/5 +hbottom+hbottom3,    0, None],
               2:[1./2*  rsand,0.0,   0*   hsand/5 +hbottom+hbottom3,    0, None],
               3:[1*     rsand,0.0,   0*   hsand/5 +hbottom+hbottom3,    0, None],
               4:[0*     rsand,0.0,   1*   hsand/5 +hbottom+hbottom3,    0, None],
               5:[1./3*  rsand,0.0,   1*   hsand/5 +hbottom+hbottom3,    0, None],
               6:[2./3*  rsand,0.0,   1*   hsand/5 +hbottom+hbottom3,    0, None],
               7:[1*     rsand,0.0,   1*   hsand/5 +hbottom+hbottom3,    0, None],
               8:[0*     rsand,0.0,   2*   hsand/5 +hbottom+hbottom3,    0, None],
               9:[1./3*  rsand,0.0,   2*   hsand/5 +hbottom+hbottom3,    0, None],
               10:[2./3* rsand,0.0,   2*   hsand/5 +hbottom+hbottom3,    0, None],
               11:[1*    rsand,0.0,   2*   hsand/5 +hbottom+hbottom3,    0, None],
               12:[0*    rsand,0.0,   3*   hsand/5 +hbottom+hbottom3,    0, None],
               13:[1./3* rsand,0.0,   3*   hsand/5 +hbottom+hbottom3,    0, None],
               14:[2./3* rsand,0.0,   3*   hsand/5 +hbottom+hbottom3,    0, None],
               15:[1*    rsand,0.0,   3*   hsand/5 +hbottom+hbottom3,    0, None],
               16:[-2./3*rsand,0.0,   3*   hsand/5 +hbottom+hbottom3,    0, None],
               17:[0*    rsand,0.0,   4*   hsand/5 +hbottom+hbottom3,    0, None],
               18:[1./3* rsand,0.0,   4*   hsand/5 +hbottom+hbottom3,    0, None],
               19:[2./3* rsand,0.0,   4*   hsand/5 +hbottom+hbottom3,    0, None],
               20:[1*    rsand,0.0,   4*   hsand/5 +hbottom+hbottom3,    0, None],
               21:[0*    rsand,0.0,   5*   hsand/5 +hbottom+hbottom3,    0, None],
               22:[1./2* rsand,0.0,   5*   hsand/5 +hbottom+hbottom3,    0, None]
               }



#SURFACE ONES
number_thermocouplesw=9;
#i= a,b, ... j: name of the thermocouples. They are strings
#BUT: f is not here, as it is not connected.
#thermocouplesw[i][0]=vertical position
#thermocouplesw[i][1]=temperature
#thermocouplesw[i][2]=indice block
thermocouplesw={'a':[hbottom+hbottom3+0          ,0,    None],
                'b':[hbottom+hbottom3+0.053      ,0,    None],
                'c':[hbottom+hbottom3+0.097      ,0,    None],
                'd':[hbottom+hbottom3+0.147      ,0,    None],
                'e':[hbottom+hbottom3+0.198      ,0,    None],
                'g':[hbottom+hbottom3+0.296      ,0,    None],
                'h':[hbottom+hbottom3+0.395      ,0,    None],
                'i':[hbottom+hbottom3+0.505      ,0,    None],
                'j':[hbottom+hbottom3+0.508      ,0,    None],
                }






def initial_interp(position,thermocouples):
#Initial temperatures inside the vessel are stored in thermocouples[i][3].
#This function executes interpolation to give initial temperature for each block of the mesh.
#position: the place you want to get temperature
#thermocouples: dictionnary of volume thermocouples
#return a float: temperature(position)
    positionradii=math.sqrt(position[0]**2+position[1]**2)
    ri=(rsand*0,rsand*1./3,rsand*1./2,rsand*2./3,rsand*1)
    zi=(hbottom3+hbottom+hsand*(0/5),hbottom3+hbottom+hsand*(1./5),hbottom3+hbottom+hsand*(2./5),hbottom3+hbottom+hsand*(3./5),hbottom3+hbottom+hsand*(4./5),hbottom3+hbottom+hsand*(5./5))    
#    Assumption: only conduction, no convection
#     Thermal Diffusion: T(r)= A ln (r) + B
#     It makes sense between the two closest points
    a1=(thermocouples[3][3]-thermocouples[2][3])/math.log(2)*math.log(2./3)+thermocouples[2][3]
    a2=(thermocouples[3][3]-thermocouples[2][3])/math.log(2)*math.log(4./3)+thermocouples[2][3]
    a3=thermocouples[5][3]+(thermocouples[6][3]-thermocouples[5][3])/math.log(2)*math.log(3./2)
    a4=thermocouples[9][3]+(thermocouples[10][3]-thermocouples[9][3])/math.log(2)*math.log(3./2)
    a5=thermocouples[13][3]+(thermocouples[14][3]-thermocouples[13][3])/math.log(2)*math.log(3./2)
    a6=thermocouples[18][3]+(thermocouples[19][3]-thermocouples[18][3])/math.log(2)*math.log(3./2)
    a7=(thermocouples[18][3]-thermocouples[17][3])/math.log(2)*math.log(2./3)+thermocouples[17][3]
    a8=0
    a9=0
    Ti=[[thermocouples[1][3],a1,thermocouples[2][3],a2,thermocouples[3][3]],
        [thermocouples[4][3],thermocouples[5][3],a3,thermocouples[6][3],thermocouples[7][3]],
        [thermocouples[8][3],thermocouples[9][3],a4,thermocouples[10][3],thermocouples[11][3]],
        [thermocouples[12][3],thermocouples[13][3],a5,thermocouples[14][3],thermocouples[15][3]],
        [thermocouples[17][3],thermocouples[18][3],a6,thermocouples[19][3],thermocouples[20][3]],
        [thermocouples[17][3],a7,thermocouples[18][3],a8,a9]
        ]
    
    newfunc=scipy.interpolate.interp2d(ri,zi,Ti,kind='cubic')
    interpolation=newfunc(positionradii,position[2])
    return interpolation



def initial_interpw(position,thermocouplesw):
#Initial temperatures around the vessel are stored in thermocouplesw[i][1].
#This function executes interpolation to give initial temperature for each block of the mesh.
#position: the place you want to get temperature
#thermocouplesw: dictionnary of surface thermocouples
#return a float: temperature(position)
    name=map(chr, range(97, 102))+map(chr,range(103,107))
    z=position[2]
    if z>hbottom3+hbottom+hsand:return thermocouplesw[name[number_thermocouples-1]][1]
    elif z<hbottom3+hbottom: return thermocouplesw[name[0]][1]
    else:
        import numpy.interp
        zi=[]
        Ti=[]
        for i in range(number_thermocouplesw):
            zi.append(thermocouplesw[name[i]][0])
            Ti.append(thermocouplesw[name[i]][1])
        interpolation=numpy.interp(z,zi,Ti)
        return interpolation
    
def findindiceblock(position, gridblock):
#find the nearest block of position place in gridblock mesh
# position: 3 elements list
#gridblock:t2block
#return indice in blocklist
#name of this block: gridblock.blocklist[indice].name
    distance_min_square=hsand**2
    indice=None
    for i in range(len(gridblock.blocklist)):
        distance_square=(position[0]-gridblock.blocklist[i].centre[0])**2+(position[1]-gridblock.blocklist[i].centre[1])**2+(position[2]-gridblock.blocklist[i].centre[2])**2
        if distance_square<distance_min_square:
            indice=i
            distance_min_square=distance_square
    return indice
#######################################################################################
#THERMAL CONDUCTIVITY
#To understand how it works: presentation.pdf, slide 19.
import math
from scipy.interpolate import interp1d
import scipy.io
from scipy.interpolate import interp2d

#We want to get k_fluid(T,P)
#Thanks to matlab structure: 'co2data3.mat'
# It return f(T,P)=k_fluid(T,P) 
mat = scipy.io.loadmat('co2data3.mat')
list1=[]
for i in range(600):
    list1.append(mat['P'][i][0])
    
list2=mat['T'][0]
f=interp2d(list2,list1,mat['tc'],)

#Density(T,P):
density=interp2d(list2,list1,mat['D'])



#X to x : image processing from pixel position(X) to real value (x)
#Y to y : image processing from pixel position(Y) to real value (y)
def Xtox(x): return math.pow(10.0,(7.0/(1217.0-64.0)*(float(x)-64.0))-2)

def Ytoy(y): return math.pow(10.0,((3.68/604.0)*(620.0-float(y)))-2.21)

phi1_X=[64,85,103,133,153,176,196,219,242,259,278,296,314,335,368,401,431,468,514,540,596,635,672,706,745,786,835,887,957,1006,1064,1116,1159,1191,1216]
phi1_Y=[76,100,116,144,165,187,206,224,245,259,274,286,296,308,320,335,345,359,373,381,398,408,415,424,432,440,451,455,467,476,482,487,493,497,499]

phi2_X=[64,82,102,122,141,163,185,213,236,257,281,305,331,352,381,409,438,466,493,520,546,576,602,631,668,697,731,771,812,853,890,937,979,1015,1059,1093,1143,1182,1216]
phi2_Y=[16,38,56,74,93,113,133,163,186,206,229,252,278,296,322,346,370,391,409,428,443,458,471,483,497,505,514,524,532,542,550,561,567,574,580,582,591,595,597]

phi1_x=[]
phi1_y=[]
for compteur in range(len(phi1_X)):
    phi1_x.append(Xtox(phi1_X[compteur]))
    phi1_y.append(Ytoy(phi1_Y[compteur]))

#phi1_interpolate(x) is phi1(ratio) function    
phi1_interpolate=interp1d(phi1_x,phi1_y)


phi2_x=[]
phi2_y=[]
for compteur in range(len(phi2_X)):
    phi2_x.append(Xtox(phi2_X[compteur]))
    phi2_y.append(Ytoy(phi2_Y[compteur]))

#phi2_interpolate(x) is phi2(ratio) function    
phi2_interpolate=interp1d(phi1_x,phi1_y)


epsilon=0.41 #Sand porosity (in the vessel)
k_solid=1 #W/(m.K)
#Pressure_=200 #Pressure you want to study.
Pressure_top=200 #Pressure at the top of the vessel

#I  assume pressure is nearly constant.
#It is possible to get a more realistic pressure (hydrostatic) using the following function:
#def Pressure(vertical_position,Pressure_top,temperature)
#It imports co2data3.mat (few lines before) in order to use a function:
#density=interp2(list2,list1,mat['D'])
#It uses a vertical discretisation: P(x-delta)=P(x)+density(temperature(x-delta),P(x))*g*delta
#It returns a more accurate pressure 

#Done:
def pressure(vertical_position,pressure_top):
    #need that thermocouple[i][2] is already full.
    delta=0.01
    step=int((hbottom+hbottom3+hsand-vertical_position)/delta)
    if vertical_position>hbottom+hbottom3+hsand:
        print 'error: not in the sand'
        return pressure_top
    elif vertical_position<hbottom+hbottom3:
        print 'error: not in the sand'
        return pressure_top
    else:
        pressure=pressure_top
        for i in range(step):
            temperature=initial_interp([0,0,vertical_position-i*step],thermocouples)
            pressure=pressure+density(temperature,pressure)*delta*9.81
    return pressure
        
            
    

def k_fluid(Temperature,Pressure):
    #Get k_fluid in operating conditions
    k_fluid=f(Temperature,Pressure)
    return k_fluid

def ratio(Temperature,Pressure):
    #Get ratio thanks to k_fluid(T,P) and constant k_solid
    ratio=k_fluid(Temperature,Pressure)/k_solid
    return ratio

def thermal_conductivity(Temperature,vertical_position,Pressure_top):
    Pressure=pressure(vertical_position,Pressure_top)
    ratio_k=ratio(Temperature,Pressure)
    phi1=phi1_interpolate(ratio_k)
     
    phi2=phi2_interpolate(ratio_k)
     
    phi=phi2+(phi1-phi2)*(epsilon-0.260)/0.216
     
    return float(phi)
    











#######################################################################################
#TIME DISCRETISAION
#Define a new class: step
#Step object own everything to launch a tough2 simulation
#It generates a .core file
#It runs the simulation and export a .listing file
#It exports a .csv file, essential to create a dictionary structure with all the results

#This class own several functions:
#initialisation:__init__
#setinitialT
#setinitialTw
#importpreviousnumericaldatas: useful if it is not the first step. It imports the previous .csv file and create a dictionary structure.
#mesh_definition: create MESH for .core file 
#vessel_definition: define vessel, a t2data() object from PyTOUGH class. This object is written as a .core file
#rocktype_definition: define all the rock_type. It creates the spatial dicretisation, using previous functions: density, thermal_conductivity
#create_sink: inactivate the rocktype_sink elements of the mesh
#parameters_definition: several parameters, used for writing .core file
#generation_definition: create a generation object, used for writing .core file
#incon definition:
#write_vtk: DONT WORK
#write: write .core file withe vessel, a t2data object (pytough)
#run: run TOUGH2 with .core file
#exportnumerical datas: import .listingfile and export a .csv file

#This object owns several variables:
#name:
#outfilecsv
#nameinfilecore
#nameoutfilelisting
#time
#vessel
#outfile
#infileexperiment
#resultat
#thermocouples
#thermocouplesw
#previousstepcsv
#mulgrid:


class step(object):
    def __init__(self,name,previousstepcsv,infileexperiment,time,outfile,outfilecsv):
        self.name=name
        self.outfilecsv=outfilecsv
        self.nameinfilecore=name+'.core'
        self.nameoutfilelisting=name+'.listing'
        self.time=time
        self.vessel=t2data()
        self.outfile=outfile
        self.infileexperiment=infileexperiment
        self.resultat={}
        self.thermocouples=thermocouples
        self.thermocouplesw=thermocouplesw
        self.previousstepcsv=previousstepcsv
        self.mulgrid=None
        
    
    def setinitialT(self):
        if self.infileexperiment!=None:
            document=str(self.infileexperiment)
            os.chdir(r'/Users/Pierre/Desktop/serveur/Serveur')
            f=open(document, 'r')
            content=f.readlines()
            data=content[6].split()
            permutation=[22,0,1,4,2,3]+range(5,22)
            for i in range(number_thermocouples):
                self.thermocouples[i+1][3]=float(data[permutation[i]])
        
    
    def setinitialTw(self):
        if self.infileexperiment!=None:
            document=str(self.infileexperiment)
            #os.chdir(r'/Users/Pierre/Desktop/serveur/Serveur')
            f=open(self.infileexperiment, 'r')
            content=f.readlines()
            data=content[6].split()
            permutation=[25,27,28,29,30,31,32,33,26]
            name=map(chr, range(97, 102))+map(chr,range(103,107))
            for i in range(number_thermocouplesw):
                self.thermocouplesw[name[i]][1]=float(data[permutation[i]])
    
    def importpreviousnumericaldatas(self):
        infile=self.previousstepcsv
        import csv
        dicti = {}
        for key, val in csv.reader(open(infile)):
            dicti[key] = val
        dic={}
        a=dicti['0']
        b=a.split('}, ')
        dic['time']=int(dicti['time'])
        dic['elements']=int(dicti['elements'])
        for i in range(len(dicti)-2):
            dic[i]={}
            a=dicti[str(i)].split('}, ')
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
        for i in range(number_thermocouples):
            indice_thermocouple=findindiceblock([self.thermocouples[i+1][0],self.thermocouples[i+1][1],self.thermocouples[i+1][2]],self.vessel.grid)
            self.thermocouples[i+1][3]=dic[len(dicti)-3][indice_thermocouple]['T']
            self.thermocouples[i+1][4]=indice_thermocouple
        
        name=map(chr, range(97, 102))+map(chr,range(103,107))
        for i in range(number_thermocouplesw):
            indice_thermocouplew=findindiceblock([rwall,0,self.thermocouplesw[name[i]][0]],self.vessel.grid)
            self.thermocouplesw[name[i]][1]=dic[len(dicti)-3][indice_thermocouplew]['T']
            self.thermocouplesw[name[i]][2]=indice_thermocouplew
            
    


    
    
    def mesh_definition(self):
        mesh=t2grid()
        
        middle=mulgrid()
        middle.create_cylinder(hmin=hbottom+hbottom3,length=hsand,radii=rwall,theta_step=theta_step,layer_step=layer_step_middle,radii_step=middle_radii_step,nameofcylinder='m',filename='middle')
        middle_grid=t2grid()
        middle_grid.fromgeo(middle)
        middle_grid.get_block_centres_defined()
        
        top=mulgrid()
        top.create_cylinder(hmin=hbottom+hbottom3+hsand,length=htop3+htop,radii=rcap,theta_step=theta_step,layer_step=layer_step_top,radii_step=top_radii_step,nameofcylinder='t',filename='top')
        top_grid=t2grid()
        top_grid.fromgeo(top)
        top_grid.get_block_centres_defined()
       
        bottom=mulgrid()
        bottom.create_cylinder(hmin=0,length=hbottom3+hbottom,radii=rcap,theta_step=theta_step,layer_step=layer_step_bottom,radii_step=bottom_radii_step,nameofcylinder='b',filename='bottom')
        bottom_grid=t2grid()
        bottom_grid.fromgeo(bottom)
        bottom_grid.get_block_centres_defined()
        
        self.mulgrid=bottom
        mesh=bottom_grid+middle_grid+top_grid
        
        mesh.create_vertical_connection(middle_grid, rwall/middle_radii_step, bottom_grid, rcap/bottom_radii_step, theta_step,rwall)
        mesh.create_vertical_connection(top_grid, rcap/top_radii_step, middle_grid, rwall/middle_radii_step, theta_step,rwall)
        
        for blk in mesh.blocklist[:]:
            blk.ahtx=0
        mesh.delete_block(mesh.blocklist[0].name)
        self.vessel.grid=mesh
        
        #print len(self.vessel.grid.blocklist)
        
        
        
    def vessel_definition(self):
        self.vessel.title=self.name
        self.vessel.filename='vessel'+self.name+'.core'
        
        
    
    def rocktype_definition(self):
        name_list=range(10)+map(chr, range(65, 76))
        self.vessel.grid.add_rocktype(rocktype('STEEL',permeability=[0.0]*3,density=8000.0,porosity=1.0e-9,conductivity=5.3,specific_heat=300))
        self.vessel.grid.add_rocktype(rocktype('PIPE_',permeability=[1.0e-9]*3,density=8000.0,porosity=0.99,conductivity=6.689e-02,specific_heat=500))
        self.vessel.grid.add_rocktype(rocktype('PIPEI',permeability=[1.0e-9]*3,density=2600.0e40,porosity=0.11,conductivity=16.3e9,specific_heat=1.0e4))
        self.vessel.grid.add_rocktype(rocktype('PIPEO',permeability=[1.0e-9]*3,density=8000.0,porosity=0.99,conductivity=6.689e-02,specific_heat=500))
         
        for i in range(len(name_list)):
                #Definir: qu'est-ce qui change?????
                #importer la temperature
            name_sand='SAND'+str(name_list[i])
            specific_heat=0
            
            conductivity=thermal_conductivity(thermocouples[i+1][3],thermocouples[i+1][2],Pressure_top)
            self.vessel.grid.add_rocktype(rocktype(name=name_sand,nad=2,permeability=[9.3e-12]*3,density=2600.0,porosity=0.3942,conductivity=conductivity,specific_heat=specific_heat,compressibility=4.5e-10))
            self.vessel.grid.rocktype[name_sand].relative_permeability={'type':7, 'parameters':[.457,.3,1,.05]}
            self.vessel.grid.rocktype[name_sand].capillarity={'type':8, 'parameters':[.457,0,5.1e-5,1e7,.999]}
        
            
        for blk in self.vessel.grid.blocklist[:]:
            blk_radii=math.sqrt(blk.centre[0]**2+blk.centre[1]**2)
            #print ' '
            #print rsand*1./6
            #print rsand*3./6
            #print rsand
            #print blk_radii
            #if blk_radii<rsand*1./6:print 'True'
            if blk.centre[2]<hpipei:
                if blk_radii<rpipe:blk.rocktype=self.vessel.grid.rocktype['PIPEI']
                else: blk.rocktype=self.vessel.grid.rocktype['STEEL']
            elif blk.centre[2]<hbottom3+hbottom:
                if blk_radii<rpipe:blk.rocktype=self.vessel.grid.rocktype['PIPE_']
                else: blk.rocktype=self.vessel.grid.rocktype['STEEL']
            elif blk.centre[2]<hbottom3+hbottom+hsand/5*(1/2):
                if blk_radii<rsand*1./4 : blk.rocktype=self.vessel.grid.rocktype['SAND0']
                if blk_radii<rsand*3./4  : blk.rocktype=self.vessel.grid.rocktype['SAND1']
                if blk_radii<rsand : blk.rocktype=self.vessel.grid.rocktype['SAND2']
                else:blk.rocktype=self.vessel.grid.rocktype['STEEL']
            elif blk.centre[2]<hbottom3+hbottom+hsand/5*(3/2):
                if blk_radii<rsand*1./6 : blk.rocktype=self.vessel.grid.rocktype['SAND3']
                if blk_radii<rsand*3./6 : blk.rocktype=self.vessel.grid.rocktype['SAND4']
                if blk_radii<rsand*5./6 : blk.rocktype=self.vessel.grid.rocktype['SAND5']
                if blk_radii<rsand : blk.rocktype=self.vessel.grid.rocktype['SAND6']
                else:blk.rocktype=self.vessel.grid.rocktype['STEEL']
            elif blk.centre[2]<hbottom3+hbottom+hsand/5*(5/2):
                if blk_radii<rsand*1./6 : blk.rocktype=self.vessel.grid.rocktype['SAND7']
                if blk_radii<rsand*3./6 : blk.rocktype=self.vessel.grid.rocktype['SAND8']
                if blk_radii<rsand*5./6 : blk.rocktype=self.vessel.grid.rocktype['SAND9']
                if blk_radii<rsand : blk.rocktype=self.vessel.grid.rocktype['SANDA']
                else:blk.rocktype=self.vessel.grid.rocktype['STEEL']
            elif blk.centre[2]<hbottom3+hbottom+hsand/5*(7/2):
                
                if blk_radii<rsand*1./6 : blk.rocktype=self.vessel.grid.rocktype['SANDB']
                if blk_radii<rsand*3./6 : blk.rocktype=self.vessel.grid.rocktype['SANDC']
                if blk_radii<rsand*5./6 : blk.rocktype=self.vessel.grid.rocktype['SANDD']
                if blk_radii<rsand : blk.rocktype=self.vessel.grid.rocktype['SANDE']
                else:blk.rocktype=self.vessel.grid.rocktype['STEEL']
            elif blk.centre[2]<hbottom3+hbottom+hsand/5*(9/2):
                if blk_radii<rsand*1./6 : blk.rocktype=self.vessel.grid.rocktype['SANDF']
                if blk_radii<rsand*3./6 : blk.rocktype=self.vessel.grid.rocktype['SANDG']
                if blk_radii<rsand*5./6 : blk.rocktype=self.vessel.grid.rocktype['SANDH']
                if blk_radii<rsand : blk.rocktype=self.vessel.grid.rocktype['SANDI']
                else:blk.rocktype=self.vessel.grid.rocktype['STEEL']
            elif blk.centre[2]<hbottom3+hbottom+hsand:
                if blk_radii< rsand*1./4: blk.rocktype=self.vessel.grid.rocktype['SANDJ']
                if blk_radii< rsand: blk.rocktype=self.vessel.grid.rocktype['SANDK']
                else:blk.rocktype=self.vessel.grid.rocktype['STEEL']
            elif blk.centre[2]<hbottom3+hbottom+hsand+hpipe_:
                if blk_radii<rpipe:blk.rocktype=self.vessel.grid.rocktype['PIPE_']
                else: blk.rocktype=self.vessel.grid.rocktype['STEEL']
            else:
                if blk_radii<rpipe:blk.rocktype=self.vessel.grid.rocktype['PIPEO']
                else: blk.rocktype=self.vessel.grid.rocktype['STEEL']
            #print blk.rocktype
    
    def create_sink(self,rocktype_sink):
        self.vessel.grid.add_block(t2block('0v001',0.0,self.vessel.grid.rocktypelist[0],centre=[0,0,0]))
        for blk in self.vessel.grid.blocklist[:]:
            if blk.rocktype.name==rocktype_sink:
                self.vessel.grid.blocklist.remove(blk)
                self.vessel.grid.blocklist.append(blk)
 
        
    
    def parameters_definition(self):
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
        
        #set TIMES
        times=[]
        
        #set INDOM
        indom={'BSINK':[1.3790e07,0,1,10], 'PIPEI':[1.3790e07,0,1,15]}
        #Set PARAM
        parameters={'max_iterations':None, 'print_level':2, 'max_timesteps':2000, 'max_duration':None, 'print_interval':1, 
                            '_option_str':'','option':(1,1,0,0,0,0,0,0,0,0,0,2,0,0,0,0,4,0,0,0,0,4,0,0,0), 'diff0':None, 'texp':None, 'tstart':0.0, 'tstop':self.time,
                            
                            'const_timestep':-1.0,'timestep':[1], 'max_timestep':None, 'print_block':None, 'gravity':9.81,
                            'timestep_reduction':None, 'scale':None, 'relative_error':1.0e-5, 'absolute_error':None, 'pivot':None,
                            'upstream_weight':None, 'newton_weight':None, 'derivative_increment':1e-8, 'default_incons':[1.3790e07,0,1,95]}
        self.vessel.multi=multi
        self.vessel.selection=selec
        self.vessel.solver=solver
        self.vessel.history_connection=hconnection
        self.vessel.history_block=hblock
        self.vessel.indom=indom
        self.vessel.start=True
        self.vessel.parameter=copy(parameters)
        self.vessel.output_times=times
                        
    def generation_definition(self,rate_inj,time_inj):
        #print self.vessel.grid.blocklist[0]
        injection=t2generator(name='inj 1',block='0b001',type='COM3',gx=None,itab='',ltab=4,rate=rate_inj,time=time_inj,ex=None,fg=None,hg=None)
        self.vessel.generator=injection
        self.vessel.generatorlist=[injection]
    
    
    def incon_definition(self):
        incon=t2incon()
        name=map(chr, range(97, 102))+map(chr,range(103,107))
        
        for blk in self.vessel.grid.blocklist[:]:
            blk_radii=math.sqrt(blk.centre[0]**2+blk.centre[1]**2)
            if str(blk.rocktype)[0:4]=='SAND':
                
                T=initial_interp(blk.centre,self.thermocouples)
                
                t=t2blockincon(block=blk.name, porosity=0.41,variable=[.1379e8,0,1,T])
                incon.add_incon(t)
            
            elif blk.centre[2]<hbottom+hbottom3:
                if blk_radii>rpipe:
                    T=thermocouplesw[name[0]][1]
                    t=t2blockincon(block=blk.name, porosity=0.41,variable=[.1379e8,0,1,T])
                    incon.add_incon(t)
                else:
                    T=Tinjection
                    t=t2blockincon(block=blk.name, porosity=0.41,variable=[.1379e8,0,1,T])
                    incon.add_incon(t)
            elif blk.centre[2]>hbottom+hbottom3 +hsand:  
                T=thermocouplesw[name[number_thermocouplesw-1]][1]
                t=t2blockincon(block=blk.name, porosity=0.41,variable=[.1379e8,0,1,T])
                incon.add_incon(t) 
            else:    
                T=initial_interpw(blk.centre,thermocouplesw)
                t=t2blockincon(block=blk.name, porosity=0.41,variable=[.1379e8,0,1,T])
                incon.add_incon(t)
        self.vessel.incon=incon
   
    def write_vtk(self,vtkfile):
        self.vessel.grid.write_vtk(geo=self.mulgrid,savefile=vtkfile,wells=False)
    
    
    def write(self):
        
        #self.vessel.grid.delete_block('0b001')
        
        self.vessel.write(filename=self.nameinfilecore)
        
            
    def run(self):
        self.vessel.run(incon_filename=self.nameinfilecore,simulator='zco2h',silent=False,save_filename=self.nameoutfilelisting)

        
    def exportnumericaldatas(self):
        outfilecsv=self.outfilecsv
        f=open(self.nameoutfilelisting, 'r')
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
        
         
        import re
        for i in range(len(separateur)-1):
            resultat[i]={}
            resultat[i]['time']=time[i]
            element=0
            for j in range(separateur[i],separateur[i+1]):
                line=content[j]
                #matchObj=re.search(r' (\S)(\S)(\S) (\S) ..(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  .(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)  (\S)(\S)(\S)(\S)(\S)(\S)',line)
                matchObj=re.search(r'(\s)[a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9](\s)[a-zA-Z0-9]....(\s)[0-9].[0-9][0-9][0-9][0-9][0-9][E]',line)
                
                if matchObj!=None and line[1]!='.':
                    element+=1
                    a=line.split()
                    if len(a)==14:
                        name=a[0]+a[1]
                        resultat[i][element]={}
                        resultat[i][element]['name']=name
                        resultat[i][element]['T']=float(a[4])
                        resultat[i][element]['P']=float(a[3])
                        resultat[i][element]['DL']=float(a[12])
                        resultat[i][element]['DG']=float(a[11])
                        resultat[i][element]['xco2aq']=float(a[8])
                    elif len(a)==13:
                        name=a[0]+a[1][0]
                        resultat[i][element]={}
                        resultat[i][element]['name']=name
                        resultat[i][element]['T']=float(a[3])
                        resultat[i][element]['P']=float(a[2])
                        resultat[i][element]['DL']=float(a[11])
                        resultat[i][element]['DG']=float(a[10])
                        resultat[i][element]['xco2aq']=float(a[7])
                    else:print 'erreur'
        import csv
        w = csv.writer(open(self.outfilecsv, "w"))
        for key, val in resultat.items():
            w.writerow([key, val])
        self.resultat=resultat
#END OF STEP OBJECT     
################################################################################


###############################################################################
#POSTPROCESSING

#use a new class of object: an EXPERIMENT

#Use several functions:
#create an object: __init__
#gatherdatas: gather all the step together, and create a whole .csv file: out.csv
#setexperimentmesure: import all the data from experiment file for post processing
#importnumericaldatas: import .csv file as a dictionnary for post-processing
#plot: an example of post-processing: figure out a comparison between experiment temperature and simulation temperature in each thermocouple place. Create a directory, with figures, out.csv and experimental file

#Use several variables:
#name: main name of the experiment
#Texperiment: Texperiment[i][j] get the experimental temperature of j thermocouple at i time
#experimentfile: the experiment file
#numericalfile: out.csv if gatherdata is used.
#dic: a dictionnary with all the simulation results:
#    examples of dic:
#        dic['time']: number of simulation time steps
#        dic['elements']:
#        dic[i]['time']:
#        dic[i][j]['T']:
#        dic[i][j]['xco2aq']:
#namedir: name of the directory created
        
class experiment(object):
    def __init__(self,name,experimentfile,namedir):
        self.name=name
        self.Texperiment={}
        self.experimentfile=experimentfile
        self.numericalfile=None
        self.dic={}
        self.namedir=namedir
        

            
            
    def gatherdatas(self, filename, steps):
        fout=open("out.csv","w")
        # first file:
        input1=filename+'0.csv'
        compteur=0
        for line in open(input1):
            if line[0]=='e':
                elements=int(line[9:])
            elif line[0]=='t':
                time=int(line[5:])
            else:
                fout.write(line)
                compteur=compteur+1
        # now the rest:
        
        digit=len(str(time-1))
        for num in range(2,steps):
            f = open(filename+str(num)+".csv")
            
            
            for line in f:
                number=time*(num-1)
                if line[0]!='e' and line[0]!='t':
                    if digit==1:
                        number_line=number+int(line[0])
                    else:
                        number_line=number+int(line[0:digit-1])
                    line=line[digit+1:]
                    
                    line=str(number_line)+','+line
                    fout.write(line)
                    compteur=compteur+1
                    
            f.close() # not really needed
        
        import csv
        times_total=compteur
        
        time_str='time'+','+str(times_total)
        elements_str='elements'+','+str(elements)
        fout.write(time_str)
        fout.write('\n')
        fout.write(elements_str)
        fout.close()
        
        fout=open("out.csv","r")
        
        self.numericalfile="out.csv"
        
        
        
        
        
                
    def setexperimentmesure(self):
        #Temperature
        #create a dictiionnary in TC order
        #Texperiment[i][j] get the j thermocouples at j time.
        #beginning for i and j at 0
        workfile=self.experimentfile
        f1=open(workfile, 'r')
        content=f1.readlines()
        permutation=[22,0,1,4,2,3]+range(5,22)
        Texperiment={}
        data=[]
        for i in range(6,len(content)):
            data.append(content[i].split())
            Texperiment[i-6]={}
            Texperiment[i-6]['time']=float(data[i-6][34])
            for j in range(number_thermocouples):
                Texperiment[i-6][j]=float(data[i-6][permutation[j]])
        
        self.Texperiment=Texperiment
            
        
    def importnumericaldatas(self):
        infile=self.numericalfile
        import csv
        dicti = {}
        for key, val in csv.reader(open(infile)):
            dicti[key] = val
        dic={}
        a=dicti['0']
        b=a.split('}, ')
        dic['time']=int(dicti['time'])
        dic['elements']=int(dicti['elements'])
        for i in range(len(dicti)-2):
            dic[i]={}
            a=dicti[str(i)].split('}, ')
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
        self.dic=dic
        
    def plot(self):
        import os
        
        if os.path.isdir(self.namedir):
            i=1
            while os.path.isdir(self.namedir+str(i)):
                i=i+1
            self.namedir=self.namedir+str(i)
        os.mkdir(self.namedir)
        import shutil
        shutil.copyfile(self.experimentfile, self.namedir+'/'+self.experimentfile)
        shutil.copyfile('out.csv', self.namedir+'/out.csv')
        os.chdir(self.namedir)
        self.setexperimentmesure()
        self.importnumericaldatas()
        color_level=['g','r','c','b','y','m']
        stage=[0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5]
        for i in range(1,number_thermocouples+1):
            xexperiment=[]
            xnumerical=[]
            yexperiment=[]
            ynumerical=[]
#             positionthermocouples=[]
#             for k in range(3):
#                 positionthermocouples.append(thermocouples[i][k])
#             indice=findindiceblock(position=positionthermocouples,gridblock=self.vessel.grid)    
            
            
            for j in range(self.dic['time']):
                
                ynumerical.append(self.dic[j][thermocouples[i][4]]['T'])
                if j==0:
                    xnumerical.append(self.dic[j]['time'])
                else:
                    delta_time=self.dic[j]['time']-self.dic[j-1]['time']
                    if delta_time>0:
                        xnumerical.append(xnumerical[-1]+delta_time)
                    else:
                        xnumerical.append(xnumerical[-1]+self.dic[0]['time'])
                    
                #xnumerical.append(self.dic[j]['time'])
                
            plt.plot(xnumerical,ynumerical,'ro',color=color_level[stage[i]],linestyle="dashed", marker="o",)
            
            for h in range(len(self.Texperiment)-3):
                yexperiment.append(self.Texperiment[h][i-1])
                xexperiment.append(self.Texperiment[h]['time'])
            plt.plot(xexperiment,yexperiment,'ro',color=color_level[stage[i]],linestyle="solid", marker="x",)
            
            xlabel='Time (s)'
            ylabel= 'Temperature (C)'
            name='tc'+str(i)+'.png'
            title='Thermocouple '+str(i)
            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel=(ylabel)
            import os
            
            plt.savefig(name,facecolor='w', edgecolor='w')
            plt.close()
#END OF EXPERIMENT CLASS
#############################################################################




#############################################################################
#MAIN CODE            

#Inputs:
Tinjection=10                          #Temperature of the CO2 input
Pressure=200                           #Pressure in the vessel. It is possible to use pressure(vertical_position,pressure_top) if thermocouples[i][2] - temperature of thermocouples - are already set
nombre_step=10                         #Number of step in time discretization
timetotal=300                          #Total time of the simulation
timestep=float(timetotal)/nombre_step  #Time for each simulation step 
experimentfile='6.18.2013_11.07 AM.txt'#Name for the experiment file
namedir='test1'                        #Name for the directory created (with figures). If the directory already exists, the name is incremented 
max_rate=1.6e-3                        #flow rate aimed during the experiment



for compteur in range(nombre_step):              #Loop: number of steps
    name='step_'+str(compteur)                   #Main name for the step
    outfile='step_'+str(compteur)+'.listing'     #Name for the TOUGH2 output file
    outfilecsv='step_'+str(compteur)+'.csv'      #Name for the .csv file
    
    if compteur==0:                              #Is it the first step?
        previousstepcsv=None                     #You do not use a previous step for initial data, but the experiment file
        time_inj=[]                              #Time dependent injection: linear for the first step from 0 to max rate
        rate_inj=[]                              #idem
        for i in range(4):                       #idem
            rate=max_rate*float(i)/4             #idem
            time_inj.append(i*float(timestep)/4) #idem
            rate_inj.append(rate)                #idem
        infileexperiment=experimentfile          #Define the experimentfile in order to set initial data    
        
    else:                                        #If it is not the first step
        previousstepcsv='step_'+str(compteur-1)+'.csv' #Define the name of .csv file previous step
        infileexperiment=None                    #Do not use the experiment file
        time_inj=[]                              #Time dependent injection: constant rate, equal to max_rate
        rate_inj=[]                              #idem
        for i in range(4):                       #idem   
            time_inj.append(i*float(timestep)/4) #idem
            rate_inj.append(max_rate)            #idem
    a=step(name,previousstepcsv,infileexperiment,timestep,outfile,outfilecsv)#Create a step object with all the previous datas

    a.vessel_definition()                       #Define t2data object (pytough)
    a.mesh_definition()                         #Define the mesh
    
    if compteur==0:                             #Is it the first step?
        a.setinitialT()                         #
        a.setinitialTw()                        #
    else: a.importpreviousnumericaldatas()      #
    
    a.rocktype_definition()                     #Defines rocktype
    a.create_sink('PIPEO')                      #Defines a sink
    a.parameters_definition()                   #Defines parameters
    a.generation_definition(rate_inj,time_inj)  #Defines a generation object
    a.incon_definition()                        #Define incon    
    #a.write_vtk('vtk')                         DO NOT WORK!!!   
    a.write()                                   #Write .core file
    a.run()                                     #Run TOUGH2
    a.exportnumericaldatas()                    #Export all the datas in a .csv file
    print name+' COMPLETED'
    
b=experiment(name=name,experimentfile=experimentfile,namedir=namedir)#Create ane experiment object
b.gatherdatas(filename='step_',steps=nombre_step)                    #Gathers all the step in a .csv file: out.csv
b.plot()                                                             #Figure out figures,experimentfile and out.csv in a new directory   
    
        
        
    
    
            



        
                    
            
            
        
        
            
            