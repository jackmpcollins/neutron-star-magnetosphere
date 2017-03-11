## Jack Collins, Final Year Project, National University of Ireland Galway, 2017
## written for Python2.7

## This file defines the class 'Star'

from visual import *
from classFieldline import *

global c
c = 3.0*(10**5) #km/s

class Star(object):
    def __init__(self, radius, period, chi):
        self.radius = radius #radius of star in km
        self.period = period #period of rotation in seconds
        self.chi = chi*pi/180 #inclination of magnetic axis to rotational in radians

        self.B_0 = 1 #magnetic field strength at magnetic pole
        self.xi = 0.9 #the dimensionless magnetic colatitude of open field lines, xi=1 corresponds to the last open field lines, xi=0 corresponds to the magnetic axis

        self.lc_radius = c*self.period/(2*pi)
        self.lc_height = 4*self.lc_radius
        self.eta_lc = self.lc_radius/self.radius

        self.omega = 2*pi/self.period #stellar angular velocity
        self.theta_0 = (self.omega*self.radius/c)**0.5 #canonical polar cap half-angle

        self.k = 0.836 #parameter measuring general-relativistic effect of frame dragging at the stellar surface in units of stellar angular velocity
        
    def draw(self, *initialAngles):
        #create scene and set as main
        self.scene = display(title="Neutron Star Magnetosphere") #, stereo="redcyan", stereodepth=0)
        self.scene.select()

        #draw star and axes
        sphere(pos=vector(0,0,0), radius=self.radius, color=color.white)
        arrow(pos=vector(0,0,self.radius), axis=vector(0,0,self.radius), color=color.red)
        arrow(pos=vector(self.radius,0,0), axis=vector(self.radius,0,0), color=color.green)

        #draw rotation axis
        schi = sin(self.chi)
        cchi = cos(self.chi) #Muslimov and Harding (2005) says phi measured counterclockwise from meridian passing through rotation axis
        self.rotationAxis = vector(schi,0,cchi)
        for i in range(-20, 21):
            arrow(pos=i*self.radius*self.rotationAxis, axis=self.radius*self.rotationAxis, color=color.blue)
        self.rotAxPerp = self.radius*vector(cchi,0,-schi)
        arrow(pos=self.rotAxPerp, axis=self.rotAxPerp, color=color.yellow)
        #draw light cylinder
        cylinder(pos=-self.lc_height/2*self.rotationAxis, axis=self.lc_height*self.rotationAxis, radius=self.lc_radius, opacity=0.15, color=color.yellow)

        #draw magnetosphere
        for theta in initialAngles:
            for phi in range(0,360,30):
                Fieldline(self, phi, theta).draw()

    def setOtherk(self):
        # k (estimate), the parameter measuring the general-relativistic effect of frame dragging at the stellar surface in units of stellar angular velocity
        '''
        M = (4.0/3)*pi*((1000*self.radius)**3)*(10**17) # NS mass (volume by density), kg
        I = (2.0/5)*M*self.radius**2 #moment of inertia of NS with radius R
        I_45 = I/(10**45) #g cm^2
        R_6 = self.radius/(10**6) #cm
        k = 0.15*I_45/(R_6**3)
        '''
        print("old k =", self.k)
        self.k = (1.0/125)*pi*self.radius**2 #above simplified down CHECK!
        print("new k =", self.k)

    def animate(self):
        self.scene.userspin = False
        self.scene.forward = self.rotAxPerp
        self.scene.up = self.rotationAxis
        while True:
            rate(30)
            currentCamera = vector(self.scene.forward)
            self.scene.forward = currentCamera.rotate(pi/90, self.rotationAxis)