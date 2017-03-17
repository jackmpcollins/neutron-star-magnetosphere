## Jack Collins, Final Year Project, National University of Ireland Galway, 2017
## written for Python2.7

## This file defines the class 'Fieldline'

from visual import *

#converts a point from 'magnetic coordinates' to coordinates with z along rotational axis
def mag2rot(pos, chi):
    m = matrix([[cos(chi),0,-sin(chi)], [0,1,0], [sin(chi),0,cos(chi)]])
    new = dot(m,pos)
    return vector(new.tolist()[0])

#changes emission from magnetic coords to rotational coords
def calculateEmissionDirection(B, chi):
    B = mag2rot(B, chi)
    theta = acos(B.z/mag(B))
    phi = atan2(B.y, B.x) + (B.y < 0)*2*pi
    return (rad2deg(phi), rad2deg(theta))

class Fieldline(object):
    def __init__(self, star, phi_deg, theta_deg, stopAtEmission=False, emissionRadPercentLC=[0.5], fromSouthPole=False, rep=6000):
        eta_lc = star.eta_lc
        k = star.k
        B_0 = star.B_0
        theta_0 = star.theta_0
        xi = star.xi

        phi = phi_deg*pi/180
        phi_0 = phi
        #print("emissionRadPercentLC: ", emissionRadPercentLC)
        emissionRadii = [item*star.lc_radius for item in emissionRadPercentLC+[100]] #+[100] gives a ridiculous outermost emission distance which will never be reached
        #print("emissionRadii: ", emissionRadii)
        lambda_0 = (1-star.k)*cos(star.chi) + 1.5*star.theta_0*star.xi*sin(star.chi)*cos(phi_0)
        #x = raw_input("continue?")
        #initial position
        theta = pi*theta_deg/180
        position = star.radius*vector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
        self.initialPosition = position

        self.isOpen = False #assume fieldline is closed until proven open

        self.emissionDirections = [] #star list of emissions

        #if following field lines backwards from south pole, need to change how to check if returned to star, and how position changes due to field
        if fromSouthPole:
            def backAtStar(theta):
                return (theta < 0.05*pi)
            def movePos(position, B):
                return (position - B)
        else:
            def backAtStar(theta):
                return (theta > 0.95*pi)
            def movePos(position, B):
                return (position + B)

        #start list of points
        self.points = []
        
        cleanExit = False #this gets set to true if for loop exits because of hitting light cylinder or returning to star. False indicates 6000 repetitions
        for i in range(rep):
            r = mag(position)
            if backAtStar(theta):
            #if r + 0.5 < star.radius: #better condition to use if visualising
                cleanExit = True
                print("Back at star")
                break

            theta = acos(position.z/r)
            phi = atan2(position.y, position.x) + (position.y < 0)*2*pi #gives positive angle in radians
            eta = (r/star.radius)
            
            dist = sqrt( position.y**2 + (r*cos(atan(position.z/position.x) + star.chi))**2)
            #print("dist: ", dist)
            if dist > star.lc_radius:
                cleanExit = True
                #print("Hit the light cylinder")
                self.isOpen = True
                break
            
            #calculate trig and common terms now so they aren't repeatedly done later
            #this also makes formulae more readable
            sin_chi = sin(star.chi)
            cos_chi = cos(star.chi)
            sin_phi = sin(phi)
            cos_phi = cos(phi)
            sin_theta = sin(theta)
            cos_theta = cos(theta)
            tan_theta = tan(theta)
            sin_theta2 = sin_theta**2
            
            #spherical unit vectors
            e_r = vector(sin_theta*cos_phi, sin_theta*sin_phi, cos_theta)
            e_theta = vector(cos_theta*cos_phi, cos_theta*sin_phi, -sin_theta)
            e_phi = vector(-sin_phi, cos_phi, 0)

            if star.chi == 0:
                #simplified case for aligned rotator
                B_r_s = cos_theta*( 1 + ((eta/eta_lc)**2)*( sin_theta2 - 2*k*(1-k) ) )
                B_theta_s = 0.5*sin_theta*( 1 - 0.5*((eta/eta_lc)**2)*( sin_theta2 - 4*k*(1-k) ) )
                B_phi_s = -(1-k)*(eta/eta_lc)*sin_theta
                B_s = (B_0/eta**3)*(B_r_s*e_r + B_theta_s*e_theta + B_phi_s*e_phi)
                B = 10*B_s/mag(B_s)

            else:
                #pure dipole magnetic field
                B_d = (B_0/eta**3)*( cos_theta*e_r + 0.5*sin_theta*e_theta )
                
                #first correction due to charge flow along open field lines
                B_r_1 = (1.5)*theta_0*xi*(1-theta/tan_theta)*sin_chi*sin_phi
                B_theta_1 = -(1.5)*theta_0*xi*((theta - sin_theta*cos_theta)/sin_theta2)*sin_chi*sin_phi
                B_phi_1 = -(1-k)*cos_chi*sin_theta + (1.5)*theta_0*xi*( (sin_theta**3 + sin_theta - theta*cos_theta)/sin_theta2 )*sin_chi*cos_phi
                B_1 = (eta/eta_lc)*(B_0/eta**3)*(B_r_1*e_r + B_theta_1*e_theta + B_phi_1*e_phi)
                
                #second correction due to rotation
                B_r_2 = cos_theta*( cos_chi*( cos_chi*sin_theta2 + 2*lambda_0*(1-k) ) + (sin_chi**2)*(1 - 0.5*sin_theta2) - 2*( cos_chi*sin_theta*cos_theta + 1.5*lambda_0*theta_0*xi*(theta/(sin_theta*cos_theta) - 1) )*sin_chi*cos_phi - 0.5*(sin_chi**2)*sin_theta2*cos(2*phi) )
                B_theta_2 = -sin_theta*( ( 0.25*cos_chi*sin_theta2 + lambda_0*(1-k) )*cos_chi + 0.5*(sin_chi**2)*(1 - 0.25*sin_theta2) + ( ( (0.25/tan_theta)*cos(2*theta) + (1 - cos_theta)/(sin_theta**3) )*cos_chi - 1.5*lambda_0*theta_0*xi*(3*(1 - theta/tan_theta)/sin_theta2 - 1) )*sin_chi*cos_phi - (5.0/8)*( (1/sin_theta2)*( (1-cos_theta)/sin_theta2 - 0.5 ) + 0.2*sin_theta2 )*(sin_chi**2)*cos(2*phi) )
                B_phi_2 = ( ( (0.25 + (1-cos_theta)/sin_theta2 )*cos_chi + 1.5*lambda_0*theta_0*xi*( (3*sin_theta*cos_theta - theta*(3-2*sin_theta2))/sin_theta2 ) )*sin_chi*sin_phi - (5.0/8)*( (1-cos_theta)/(sin_theta**3) - (1/(2*sin_theta)) )*(sin_chi**2)*sin(2*phi) )
                B_2 = ((eta/eta_lc)**2)*(B_0/eta**3)*(B_r_2*e_r + B_theta_2*e_theta + B_phi_2*e_phi)
                
                #Effect of ExB drift
                B_r_3 = -2*(cos_chi*cos_theta + sin_chi*sin_theta*cos_phi)
                B_theta_3 = cos_chi*sin_theta - sin_chi*cos_theta*cos_phi
                B_phi_3 = sin_chi*sin_phi
                B_3 = ((eta/eta_lc)**2)*(B_0/eta**3)*lambda_0*(B_r_3*e_r + B_theta_3*e_theta + B_phi_3*e_phi)
                
                #calculate total magnetic field
                B = 10*(B_d + B_1 + B_2 + B_3)/mag(B_d + B_1 + B_2 + B_3)
            #print("emissioRadii[0]: ", emissionRadii[0])

            #check if at an emission distance, and if there then remove that distance from check
            if r > emissionRadii[0]:
                #print("Emission!")
                self.emissionDirections.append(calculateEmissionDirection(B, star.chi))
                del emissionRadii[0]
                if stopAtEmission:
                    cleanExit = True
                    break
            
            self.points.append(position)

            #current position
            position = movePos(position, B)

        if not cleanExit:
            print("Hit "+str(rep)+" rep limit")

    #draw the fieldline
    def draw(self):
        curve(pos=self.points, radius=5)