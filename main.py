## Jack Collins, Final Year Project, National University of Ireland Galway, 2017
## written for Python2.7

from visual import *
from time import time

'''
## PARAMETERS
R = 10 #radius of star in km
P = 0.3 #period of rotation in seconds
chi = 0 #inclination of magnetic axis to rotational degrees
B_0 = 1 #magnetic field strength at magnetic pole
xi = 0.9 #the dimensionless magnetic colatitude of open field lines

#create scene and set as main
scene = display(title="Neutron Star Magnetosphere") #, stereo="redcyan", stereodepth=0)
scene.select()

## CALCULATE
# stellar angular velocity
omega = 2*pi/P

# theta_0, the canonical PC half-angle
c = 3*10**5 # speed of light, km/s
theta_0 = (omega*R/c)**0.5
print("theta_0:", theta_0)

# chi in radians
chi = chi*pi/180

# k (estimate), the parameter measuring the general-relativistic effect of frame dragging at the stellar surface in units of stellar angular velocity

M = (4/3)*pi*((1000*R)**3)*(10**17) # NS mass (volume by density), kg
I = (2/5)*M*R**2 #moment of inertia of NS with radius R
I_45 = I/(10**45) #g cm^2
R_6 = R/(10**6) #cm
k = 0.15*I_45/(R_6**3)

#above simplified down by hand
k = 0.008*pi*(R**2)

#
#phi_0 = ? #HOW TO CALCULATE THIS?
phi_0 = 25*pi/180
lambda_0 = (1-k)*cos(chi) + (3/2)*theta_0*xi*sin(chi)*cos(phi_0)
print("lambda_0:", lambda_0)
'''
def drawBasics():
    # star and rotational axis
    sphere(pos=vector(0,0,0), radius=R, color=color.white)
    arrow(pos=vector(0,0,R), axis=vector(0,0,R), color=color.red)
    arrow(pos=vector(R,0,0), axis=vector(R,0,0), color=color.green)

    # rotation axis
    schi = -sin(chi)
    cchi = cos(chi) #may be other side of magnetic axis - investigate
    rotationAxis = vector(schi,0,cchi)
    for i in range(-20, 21):
        arrow(pos=i*R*rotationAxis, axis=R*rotationAxis, color=color.blue)

    # light cylinder
    global lightCylinderRadius
    lightCylinderRadius = c*P/(2*pi)
    print("lc radius:", lightCylinderRadius)
    lightCylinderHeight = 4*lightCylinderRadius
    global eta_lc
    eta_lc = lightCylinderRadius/R
    print("eta_lc:", eta_lc)
    cylinder(pos=-lightCylinderHeight/2*rotationAxis, axis=lightCylinderHeight*rotationAxis, radius=lightCylinderRadius, opacity=0.15, color=color.yellow)

#samples the calculated points to reduce amount to plot, and plots curve
def samplePointsPlotCurve(points, initialAngle):
    usedPoints = []
    for i in range(len(points)):
        if i % 30 == 0:
            usedPoints.append(points[i])
    curve(pos=points, radius=25/initialAngle)


#draws a single magnetic field line given it's starting point on the surface of the star
def drawMagLine(phi, theta_init_deg, rep):
    phi_0 = phi
    ## PARAMETERS

    #cartesian unit vectors
    x = vector(1,0,0)
    y = vector(0,1,0)
    z = vector(0,0,1)
    
    #initial position
    theta_init = pi*theta_init_deg/180
    x_pos = sin(theta_init)*cos(phi)*x
    y_pos = sin(theta_init)*sin(phi)*y
    z_pos = cos(theta_init)*z
    position = R*(x_pos + y_pos + z_pos)

    #start a curve
    #fieldCurve = curve(pos=[position], radius=25/theta_init_deg)

    #start list of points
    points = [position]
    
    for i in range(rep):
        #current values
        r = mag(position)
        #check if loop complete
        if r + 0.5 < R:
            print("r:",r)
            print("back at star")
            break
        theta = acos(position.z/r)
        phi = atan(position.y/position.x)
        if position.x < 0:
                    phi += pi
        eta = (r/R)
        if r*sin(theta) > lightCylinderRadius:
        	print("hit light cylinder")
        	break
        
        #calculate trig and common terms now so they aren't repeatedly done later
        #this also makes formulae more readable
        sin_chi = sin(chi)
        cos_chi = cos(chi)
        sin_phi = sin(phi)
        cos_phi = cos(phi)
        sin_theta = sin(theta)
        cos_theta = cos(theta)
        tan_theta = tan(theta)
        sin_theta2 = sin_theta**2
        
        ##PRECALCULATED TRIG
        #spherical vectors
        e_r = sin_theta*cos_phi*x + sin_theta*sin_phi*y + cos_theta*z
        e_theta = cos_theta*cos_phi*x + cos_theta*sin_phi*y - sin_theta*z
        e_phi = -sin_phi*x + cos_phi*y
        
        if chi == 0:
            #simplified case for aligned rotator
            B_r_s = cos_theta*( 1 + ((eta/eta_lc)**2)*( sin_theta2 - 2*k*(1-k) ) )
            B_theta_s = (0.5)*sin_theta*( 1 - (0.5)*((eta/eta_lc)**2)*( sin_theta2 - 4*k*(1-k) ) )
            B_phi_s = -(1-k)*(eta/eta_lc)*sin_theta
            B_s = (B_0/eta**3)*(B_r_s*e_r + B_theta_s*e_theta + B_phi_s*e_phi)
            B = 10*B_s/mag(B_s)

        else:
            #pure dipole magnetic field
            B_d = (B_0/eta**3)*( cos_theta*e_r + 0.5*sin_theta*e_theta )
            
            #first correction due to charge flow along open field lines
            B_r_1 = (3/2)*theta_0*xi*(1-theta/tan_theta)*sin_chi*sin_phi
            B_theta_1 = -(3/2)*theta_0*xi*(theta - sin_theta*cos_theta/sin_theta2)*sin_chi*sin_phi
            B_phi_1 = -( (1-k)*cos_chi*sin_theta + (3/2)*theta_0*xi*( (sin_theta**3 + sin_theta - theta*cos_theta)/sin_theta2 )*sin_chi*cos_phi )
            B_1 = (eta/eta_lc)*(B_0/eta**3)*(B_r_1*e_r + B_theta_1*e_theta + B_phi_1*e_phi)
            
            #second correction due to rotation
            B_r_2 = cos_theta*( cos_chi*( cos_chi*sin_theta2 + 2*lambda_0*(1-k) ) + (sin_chi**2)*(1 - 0.5*sin_theta2) - 2*( cos_chi*sin_theta*cos_theta + (3/2)*lambda_0*theta_0*xi*(theta/(sin_theta*cos_theta) - 1) )*sin_chi*cos_phi - 0.5*(sin_chi**2)*(sin_theta2)*cos(2*phi) )
            B_theta_2 = -sin_theta*( ( 0.25*cos_chi*sin_theta2 + lambda_0*(1-k) )*cos_chi + 0.5*(sin_chi**2)*(1 - 0.25*sin_theta2) + ( ( (0.25/tan_theta)*cos(2*theta) + (1 - cos_theta)/(sin_theta**3) )*cos_chi - (3/2)*lambda_0*theta_0*xi*(3*(1 - theta/tan_theta)/sin_theta2 - 1) )*sin_chi*cos_phi - (5/8)*( (1/sin_theta2)*( (1-cos_theta)/sin_theta2 - 0.5 ) + 0.2*sin_theta2 )*sin_theta2*cos(2*theta) )
            B_phi_2 = ( ( (0.25 + (1-cos_theta)/sin_theta2 )*cos_chi + (3/2)*lambda_0*theta_0*xi*((3*sin_theta*cos_theta - theta*(3-2*sin_theta2))/sin_theta2) )*sin_chi*sin_phi - (5/8)*( (1-cos_theta)/(sin_theta**3) - (2*sin_theta)**(-1) )*sin_theta2*sin(2*phi) )
            B_2 = ((eta/eta_lc)**2)*(B_0/eta**3)*(B_r_2*e_r + B_theta_2*e_theta + B_phi_2*e_phi)
            
            #Effect of ExB drift
            B_r_3 = -2*(cos_chi*cos_theta + sin_chi*sin_theta*cos_phi)
            B_theta_3 = cos_chi*sin_theta - sin_chi*cos_theta*cos_phi
            B_phi_3 = sin_chi*sin_phi
            B_3 = ((eta/eta_lc)**2)*(B_0/eta**3)*lambda_0*(B_r_3*e_r + B_theta_3*e_theta + B_phi_3*e_phi)
            
            #calculate total magnetic field
            B = 10*(B_d + B_1 + B_2 + B_3)/mag(B_d + B_1 + B_2 + B_3)
        
        #fieldCurve.append(position)
        points.append(position)

        #current position
        position = position + B

        #print feedback
        '''
        print("\n x, z, theta:", position.x, position.z, theta*180/pi)
        print("position:", position)
        print("mag vector:", B)
        print("mag field strength:", mag(B))
        '''
    samplePointsPlotCurve(points, theta_init_deg)

def allPhi(phiRange):
    for j in range(phiRange):
        phi = pi*j/10
        #drawMagLine(phi, 25, 1000)
        #drawMagLine(phi, 20, 42)
        #drawMagLine(phi, 15, 75)
        drawMagLine(phi, 10, 2000)
        drawMagLine(phi, 5, 2000)
        #drawMagLine(phi, 4, 3000)


def main():
    global R, P, chi, c, B_0, xi, theta_0, k, lambda_0, eta_lc
    #load list of already saved pulsars
    starInfo = {}
    try:
        f = open("savedStars.txt", 'r')
        print("Saved stars are ")
        for line in f.readlines():
            info = line.split(",")
            starName = info.pop(0)
            print(starName)
            starInfo[starName] = info
        f.close()
    except:
        print("No saved stars")

    #ask user input
    starName = raw_input("Type name of saved pulsar or \'new\': ")
    if starName == "new":
        R = float(raw_input("Radius of star (km): "))
        P = float(raw_input("Period of rotation (s): "))
        chi = float(raw_input("Inclination of magnetic axis to rotational axis (degrees): "))
        chi = chi*pi/180 #convert to radians
        name = raw_input("If you would like to save this star enter a name (to skip hit Enter): ")
        if name != "":
            f = open("savedStars.txt", 'w')
            f.write(",".join([name,R,P,chi]))
    else:
        (R, P, chi) = [float(item) for item in starInfo[starName]]

    B_0 = 1 #magnetic field strength at magnetic pole
    xi = 0.9 #the dimensionless magnetic colatitude of open field lines

    #create scene and set as main
    scene = display(title="Neutron Star Magnetosphere") #, stereo="redcyan", stereodepth=0)
    scene.select()

    # stellar angular velocity
    omega = 2*pi/P

    # theta_0, the canonical PC half-angle
    c = 3*10**5 # speed of light, km/s
    theta_0 = (omega*R/c)**0.5

    # k (estimate), the parameter measuring the general-relativistic effect of frame dragging at the stellar surface in units of stellar angular velocity
    '''
    M = (4/3)*pi*((1000*R)**3)*(10**17) # NS mass (volume by density), kg
    I = (2/5)*M*R**2 #moment of inertia of NS with radius R
    I_45 = I/(10**45) #g cm^2
    R_6 = R/(10**6) #cm
    k = 0.15*I_45/(R_6**3)
    '''
    #above simplified down by hand
    k = 0.008*pi*(R**2)
    k = 0.836
    print("k:", k)

    #phi_0 = ?
    phi_0 = 25*pi/180
    lambda_0 = (1-k)*cos(chi) + (3/2)*theta_0*xi*sin(chi)*cos(phi_0)
    print("lambda_0:", lambda_0)

    startTime = time()
    #draw star, rotational axis, light cylinder
    drawBasics()
    #calculate points using values given
    allPhi(21)
    #plot field
    #output instructions for manipulation
    print("To zoom in/out, hold down the scroll wheel and move the mouse forward/backward.")
    print("To move star, hold down the right mouse button to grab, and move the mouse.")
    #ask user if want to do a different pulsar
    print("Time elapsed:", time() - startTime)


main()