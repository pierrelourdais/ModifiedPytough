'''
Created on 24 mai 2013

@author: Pierre
Tu use with Pytough
'''


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

f=open('vessel.listing', 'r')
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
resultat={}
resultat['time']=len(separateur)
resultat['elements']=elements
# for i in range(len(separateur)):
for i in range(len(separateur)):
    difference=0
    while (content[separateur[i]+difference]+'None').split()[0]!='at':
        difference+=1
    resultat[i]={}
    resultat[i]['time']=time[i]
    for j in range(elements):
        blank=0
        test=content[separateur[i]+difference+j]
        
        test1=test.split()
        
        while (len(test1)==0 or test1[-1]!='0.00'):
            blank+=1
            test=content[separateur[i]+difference+j+blank]
            test1=test.split()
            

            
            
        resultat[i][j]={}
        a=content[separateur[i]+j+blank+difference].split()
        resultat[i][j]['T']=float(a[4])
        resultat[i][j]['P']=float(a[3])    
        resultat[i][j]['name']=a[0]+a[1] 
        resultat[i][j]['DL']=float(a[12]) 
        resultat[i][j]['DG']=float(a[11])
        resultat[i][j]['xco2aq']=float(a[8])

return resultat

