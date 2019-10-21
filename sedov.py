# -*- coding: utf-8 -*-
"""
Created on Wed Oct 01 10:21:01 2014

@author: QuantumMonkey
"""


import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
#from sympy import *
import ps2 


pi = 3.141592653
small = 1e-5 #define a really small number

# Put gamma to be global variable instead. Change the value of gamma number here
gamma = 5.0/3.0
#gamma = 1.0001
 
# sets the minimum and maximum possibility for u (otherwise get imaginary numbers)
umin = (gamma+1.0)/(2.0*gamma) 
umax1 = 5.0*(gamma+1.0) *0.5 / (3.0*gamma-1.0) 
umax2 = (gamma+1.0) / 2.0 
umax = min(umax1,umax2)
#print umin,umax


#Expression that gives the value of \ita a.k.a. Eqn (12)
def func(u1) : 
    v1 = (13.0*(gamma**2.0) - 7.0*gamma + 12.0)/(2.0*gamma + 1.0)/( 3.0*gamma - 1.0)
    v2 = -5.0*( gamma -1.0) /(2.0*gamma +1.0)

    chk1 = (u1**2.0 )
    chk2 = ((5.0*gamma + 5.0- 2.0* u1* (3.0*gamma-1.0))/(7.0-gamma))**v1
    chk3  = ((2.0*gamma*u1- gamma - 1.0)/(gamma-1.0))**v2 
#    print chk1,chk2,chk3
    rhs = chk1 * chk2 *chk3
    return rhs
#    return rhs - ita**(-5.0)   


    
def findu(func,c) : # this is the equation we want to solve : func(u') = c
    flist = []
    urange = np.linspace(umin,umax,100)
    for u in urange:
        flist.append(math.log(func(u)))
#    print flist
        
#    checking, plotting        
#    plt.plot(urange,flist)
#    plt.xlabel('u')
#    plt.ylabel ('log RHS')
#    plt.xlim(umin,umax)
#    plt.show()

    for i in range(99) :
#        print flist[i]
        if (flist[i] > c) & (flist[i+1] < c) :
            return urange[i]
     
    return urange[0]
            

# get the expression of u',\rho',p'
def expressions(ita) :
    num = ita**(-5.0)
#    u1 = ps2.bisec(func,umin+small,umax-small,num) #giving ita, solve for u'
    u1 = findu(func,num)
#    print "u' : ", u1
    v1 = (13.0*gamma**2.0 - 7.0*gamma + 12.0)/(2.0*gamma + 1.0)/( 3.0*gamma - 1.0)
    v3 = 3.0 /(2.0* gamma +1.0 )
    v4 = v1/(2.0 - gamma)
    v5 = 2.0/(gamma - 2.0)
    v6 = (13.0*gamma**2.0 - 7.0*gamma + 12.0)/5.0/(2.0 -gamma)/( 3.0*gamma - 1.0)
    v7 = gamma/(gamma - 2.0)

    rho1 = ((2.0*gamma*u1-gamma - 1.0)/(gamma-1.0))**v3 * \
        ((5.0*(gamma+1.0)-2.0*u1*(3.0*gamma - 1.0))/(7.0-gamma))**v4 * \
        ((gamma+1.0-2.0*u1)/(gamma-1.0))**v5  # given the u' solved from before, get rho' p'
    p1 = ita **(-2.0) * ((5.0*(gamma+1.0)-2.0*u1*(3.0*gamma-1.0))/(7.0-gamma))**v6 * \
        ((gamma+1.0-2.0*u1)/(gamma-1.0))**v7 * u1**1.2
  
    return u1,rho1,p1
    
    
# find the value of \xi_0
def findxi0(ita) :
    u,rho,p = expressions(ita)
    integ = (rho* u**2.0 + p )* ita**4.0         
    xi0 = (integ*32.0*pi/25.0/(gamma**2.0 -1))**(-0.2)
    
    return xi0

# Eqn (11). Here output the u,\rho,p that we want    
def compute(ita):
    u1,rho1,p1 = expressions(ita)
        
    ut = ita*u1
    rhot = rho1
    pt = ita**2.0*p1
    tt= pt/rhot  
    
    return ut,rhot,pt,tt
    
def mass(ita) : # integrate the non-dimensional mass
    u,rho,p = expressions(ita)
    ex = 3.0 * ita**2.0 *rho
    return ex

    
    
# main function: Plotting the variables with \ita and m(\ita)
ulist = []
rholist = []
plist = []
tlist = []
mlist = []
(xi,err) = quad(findxi0,small,1.0)
print "the xi0 value for this case is : ", xi
italist = np.linspace(small, 1., 1000) # cannot start from zero!
#print italist
for ita in italist :
#    print 'ita :', ita
    ui,rhoi,pi,ti = compute(ita)    
#    (mita,err) = quad(mass,small,ita)
#    print mita
#    mlist.append(mita)
#    ulist.append(ui)
#    rholist.append(rhoi)
#    plist.append(pi)
#    tlist.append(ti)    
#    print ui,rhoi,pi,ti
 
# Plotting with function of ita   
plt.subplot(1,2,1)     
plt.plot(italist,ulist,label = 'Velocity',linewidth = 2.0)
plt.plot(italist,rholist,label = 'Density',linestyle ='--',linewidth = 2.0)
plt.plot(italist,plist, label = 'Pressure',linestyle = ':',linewidth = 2.0)
plt.plot(italist,tlist,label = 'Temperature',linestyle = '-.',linewidth = 2.0)
plt.legend(loc=2)
plt.xlim(0,1)
plt.ylim(-0.1,2)
plt.xlabel('${\eta}$')
plt.title('${\gamma}$ = 1.0001')

 # with function of mass
plt.subplot(1,2,2)
plt.plot(mlist,ulist,label = 'Velocity',linewidth = 2.0)
plt.plot(mlist,rholist,label = 'Density',linestyle = '--',linewidth = 2.0)
plt.plot(mlist,plist, label = 'Pressure',linestyle = ':',linewidth = 2.0)
plt.plot(mlist,tlist,label = 'Temperature',linestyle = '-.',linewidth = 2.0)
plt.legend(loc=2)
plt.xlabel('M(${\eta}$)')
plt.ylim(0,2)
plt.show()  

    


        



    
    