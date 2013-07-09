#!/usr/bin/python

'''
Created on 18 juin 2013

@author: Pierre
'''

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
theta_step=4
layer_step_middle=5
layer_step_top=3
layer_step_bottom=3
middle_radii_step=5
top_radii_step=int(middle_radii_step/rwall*rcap)
bottom_radii_step=int(middle_radii_step/rwall*rcap)


#######################################################################################
#Initialisation Thermocouples


number_thermocouples=22
#thermocouples[i][0]=radial position
#thermocouples[i][1]=theta position
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






number_thermocouplesw=9;
#thermocouples[i][0]=vertical position
#thermocouples[i][1]=temperature
#thermocouples[i][2]=indice block
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
#    Ecrire une fonction qui prend en argument thermocouples[i][2]
#    Qui effectue une interpolation logarithmique radiale
#    Qui effectue une interpolation verticale lineaire
#    Qui retourne la temperarure interpolee
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
    #find the nearest block
    #return indice in blocklist
    distance_min_square=hsand**2
    indice=None
    for i in range(len(gridblock.blocklist)):
        distance_square=(position[0]-gridblock.blocklist[i].centre[0])**2+(position[1]-gridblock.blocklist[i].centre[1])**2+(position[2]-gridblock.blocklist[i].centre[2])**2
        if distance_square<distance_min_square:
            indice=i
            distance_min_square=distance_square
    return indice
#######################################################################################
import math
from scipy.interpolate import interp1d
import scipy.io


from scipy.interpolate import interp2d

mat = scipy.io.loadmat('co2data3.mat')


list1=[]
for i in range(600):
    list1.append(mat['P'][i][0])

list2=mat['T'][0]
f=interp2d(list2,list1,mat['tc'],)


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
    
phi1_interpolate=interp1d(phi1_x,phi1_y)
print phi1_interpolate

phi2_x=[]
phi2_y=[]
for compteur in range(len(phi2_X)):
    phi2_x.append(Xtox(phi2_X[compteur]))
    phi2_y.append(Ytoy(phi2_Y[compteur]))
    
phi2_interpolate=interp1d(phi1_x,phi1_y)


epsilon=0.41
k_solid=1 #W/(m.K)
Pressure=200



def thermal_conductivity(Temperature,Pressure):
    ratio_k=ratio(Temperature,Pressure)
    phi1=phi1_interpolate(ratio_k)
     
    phi2=phi2_interpolate(ratio_k)
     
    phi=phi2+(phi1-phi2)*(epsilon-0.260)/0.216
     
    return float(phi)
    

def ratio(Temperature,Pressure):
    ratio=k_fluid(Temperature,Pressure)/k_solid
    return ratio

def k_fluid(Temperature,Pressure):
    k_fluid=f(Temperature,Pressure)
    return k_fluid








#######################################################################################

    
    
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
            
            conductivity=thermal_conductivity(thermocouples[i+1][3],Pressure)
            self.vessel.grid.add_rocktype(rocktype(name=name_sand,nad=2,permeability=[9.3e-12]*3,density=2600.0,porosity=0.3942,conductivity=conductivity,specific_heat=specific_heat,compressibility=4.5e-10))
            self.vessel.grid.rocktype[name_sand].relative_permeability={'type':7, 'parameters':[.457,.3,1,.05]}
            self.vessel.grid.rocktype[name_sand].capillarity={'type':8, 'parameters':[.457,0,5.1e-5,1e7,.999]}
        
            
        for blk in self.vessel.grid.blocklist[:]:
            blk_radii=math.sqrt(blk.centre[0]**2+blk.centre[1]**2)
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
                            '_option_str':'','option':(1,1,0,0,0,0,0,0,0,0,0,2,0,0,0,0,4,0,0,0,0,4,0,0,0), 'diff0':None, 'texp':None, 'tstart':0.0, 'tstop':20000,
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
        injection=t2generator(name='inj 1',block='at 0',type='COM3',gx=None,itab='',ltab=4,rate=rate_inj,time=time_inj,ex=None,fg=None,hg=None)
        self.vessel.generator=injection
        self.vessel.generatorlist=[injection]
    
    
    def incon_definition(self):
        incon=t2incon()
        name=map(chr, range(97, 102))+map(chr,range(103,107))
        print self.thermocouples
        for blk in self.vessel.grid.blocklist[:]:
            blk_radii=math.sqrt(blk.centre[0]**2+blk.centre[1]**2)
            if str(blk.rocktype)[0:4]=='SAND':
                
                T=initial_interp(blk.centre,self.thermocouples)
                print T
                t=t2blockincon(block=blk.name, porosity=0,variable=[.1379e8,0,1,T])
                incon.add_incon(t)
            
            elif blk.centre[2]<hbottom+hbottom3:
                if blk_radii>rpipe:
                    T=thermocouplesw[name[0]][1]
                    t=t2blockincon(block=blk.name, porosity=0,variable=[.1379e8,0,1,T])
                    incon.add_incon(t)
                else:
                    T=Tinjection
                    t=t2blockincon(block=blk.name, porosity=0,variable=[.1379e8,0,1,T])
                    incon.add_incon(t)
            elif blk.centre[2]>hbottom+hbottom3 +hsand:  
                T=thermocouplesw[name[number_thermocouplesw-1]][1]
                t=t2blockincon(block=blk.name, porosity=0,variable=[.1379e8,0,1,T])
                incon.add_incon(t) 
            else:    
                T=initial_interpw(blk.centre,thermocouplesw)
                t=t2blockincon(block=blk.name, porosity=0,variable=[.1379e8,0,1,T])
                incon.add_incon(t)
        self.vessel.incon=incon
   
    def write_vtk(self,vtkfile):
        self.vessel.grid.write_vtk(geo=self.mulgrid,savefile=vtkfile,wells=False)
    
    
    def write(self):
        print self.vessel.grid.blocklist[0]
        self.vessel.grid.delete_block('0b001')
        
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
        print resultat['time']
         
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
        for line in open(input1):
            if line[0]=='e':
                elements=int(line[9:])
            elif line[0]=='t':
                time=int(line[5:])
            else:
                fout.write(line)
        # now the rest:
        digit=len(str(time))
        import csv
        
         
        for num in range(2,steps):
            f = open(filename+str(num)+".csv")
            
            
            for line in f:
                number=time*(num-1)
                if line[0]!='e' and line[0]!='t':
                    number_line=number+int(line[0:digit-1])
                    line=line[digit:]
                    
                    line=str(number_line)+','+line
                    fout.write(line)
                    
            f.close() # not really needed
        time_str='time'+','+str(time)
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
                xnumerical.append(self.dic[j]['time'])
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
            
            
Tinjection=40
Pressure=200
nombre_step=2
timetotal=200
timestep=float(timetotal)/nombre_step
experimentfile='6.18.2013_11.07 AM.txt'
namedir='test1'
for compteur in range(nombre_step):
    name='step_'+str(compteur)
    outfile='step_'+str(compteur)+'.listing'
    outfilecsv='step_'+str(compteur)+'.csv'
    
    if compteur==0:
        previousstepcsv=None
        time_inj=[]
        rate_inj=[]
        for i in range(10):
            rate=1.6e-3*float(i)/10
            time_inj.append(float(timestep)/10)
            rate_inj.append(rate)
        infileexperiment=experimentfile
        
    else: 
        previousstepcsv='step_'+str(compteur-1)+'.csv'
        infileexperiment=None
        time_inj=[]
        rate_inj=[]
        for i in range(10):
            time_inj.append(float(timestep)/10)
            rate_inj.append(1.6e-3)
    a=step(name,previousstepcsv,infileexperiment,timestep,outfile,outfilecsv)
    print 'compteur'+str(compteur)
    print a.infileexperiment
    
    a.setinitialT()
    a.setinitialTw()
    a.vessel_definition()
    a.mesh_definition()
    if compteur==0:
        a.setinitialT()
        a.setinitialTw()
    else: a.importpreviousnumericaldatas()
    
    a.rocktype_definition()
    a.parameters_definition()
    a.generation_definition(rate_inj,time_inj)
    a.incon_definition()
    #a.write_vtk('vtk')
    a.write()
    a.run()
    a.exportnumericaldatas()
    
b=experiment(name=name,experimentfile=experimentfile,namedir=namedir)
b.gatherdatas(filename='step_',steps=nombre_step)
b.plot()
    
        
        
    
    
            



        
                    
            
            
        
        
            
            