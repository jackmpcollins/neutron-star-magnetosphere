## Jack Collins, Final Year Project, National University of Ireland Galway, 2017
## written for Python2.7

## This file defines the class 'Star'

from visual import *
from math import pi

global c
c = 3*(10**5)

class Star(object):
    def __init__(self, radius, period, chi):
        self.radius = radius #radius of star in km
        self.period = period #period of rotation in seconds
        self.chi = chi*pi/180 #inclination of magnetic axis to rotational in radians

        self.B_0 = 1 #magnetic field strength at magnetic pole
        self.xi = 0.9 #the dimensionless magnetic colatitude of open field lines

        self.lc_radius = c*self.period/(2*pi)
        self.lc_height = 4*self.lc_radius
        self.eta_lc = self.lc_radius/self.radius

        self.omega = 2*pi/self.period #stellar angular velocity
        self.theta_0 = (self.omega*self.radius/c)**0.5 #canonial polar cap half-angle

        self.k = 0.836 #parameter measuring general-relativistic effect of frame dragging at the stellar surface in units of stellar angular velocity
        self.phi_0 = 25*pi/180
        self.lambda_0 = (1-self.k)*cos(self.chi) + (3/2)*self.theta_0*self.xi*sin(self.chi)*cos(self.phi_0)

    def draw(self):
        #create scene and set as main
        scene = display(title="Neutron Star Magnetosphere") #, stereo="redcyan", stereodepth=0)
        scene.select()

        #draw star and axes
        sphere(pos=vector(0,0,0), radius=self.radius, color=color.white)
        arrow(pos=vector(0,0,self.radius), axis=vector(0,0,self.radius), color=color.red)
        arrow(pos=vector(self.radius,0,0), axis=vector(self.radius,0,0), color=color.green)

        #draw rotation axis
        schi = -sin(self.chi)
        cchi = cos(self.chi) #may be other side of magnetic axis - investigate
        rotationAxis = vector(schi,0,cchi)
        for i in range(-20, 21):
            arrow(pos=i*self.radius*rotationAxis, axis=self.radius*rotationAxis, color=color.blue)

        #draw light cylinder
        cylinder(pos=-self.lc_height/2*rotationAxis, axis=self.lc_height*rotationAxis, radius=self.lc_radius, opacity=0.15, color=color.yellow)

# k (estimate), the parameter measuring the general-relativistic effect of frame dragging at the stellar surface in units of stellar angular velocity
'''
M = (4/3)*pi*((1000*R)**3)*(10**17) # NS mass (volume by density), kg
I = (2/5)*M*R**2 #moment of inertia of NS with radius R
I_45 = I/(10**45) #g cm^2
R_6 = R/(10**6) #cm
k = 0.15*I_45/(R_6**3)
'''