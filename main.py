
from visual import *

## PARAMETERS
R = 10 #radius of star in km
P = 0.3 #period of rotation in seconds
chi = 30 #inclination of magnetic axis to rotational degrees
B_0 = 1 #magnetic field strength at magnetic pole
xi = 0.9 #the dimensionless magnetic colatitude of open field lines


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
'''
M = (4/3)*pi*((1000*R)**3)*(10**17) # NS mass (volume by density), kg
I = (2/5)*M*R**2 #moment of inertia of NS with radius R
I_45 = I/(10**45) #g cm^2
R_6 = R/(10**6) #cm
k = 0.15*I_45/(R_6**3)
'''
#above simplified down by hand
k = 0.008*pi*(R**2)

#
#phi_0 = ? #HOW TO CALCULATE THIS?
phi_0 = 25*pi/180
lambda_0 = (1-k)*cos(chi) + (3/2)*theta_0*xi*sin(chi)*cos(phi_0)
print("lambda_0:", lambda_0)


## DRAW
# star and rotational axis
star = sphere(pos=vector(0,0,0), radius=R, color=color.white)
arrow(pos=vector(0,0,R), axis=vector(0,0,R), color=color.red)
arrow(pos=vector(R,0,0), axis=vector(R,0,0), color=color.green)

# rotation axis
schi = -sin(chi)
cchi = cos(chi) #may be other side of magnetic axis - investigate
rotationAxis = vector(schi,0,cchi)
for i in range(-20, 21):
    arrow(pos=i*R*rotationAxis, axis=R*rotationAxis, color=color.blue)

# light cylinder
lightCylinderRadius = c*P/(2*pi)
print("lc radius:", lightCylinderRadius)
lightCylinderHeight = 4*lightCylinderRadius
eta_lc = lightCylinderRadius/R
print("eta_lc:", eta_lc)
cylinder(pos=-lightCylinderHeight/2*rotationAxis, axis=lightCylinderHeight*rotationAxis, radius=lightCylinderRadius, opacity=0.15, color=color.yellow)

def drawMagLine(phi, theta_init_deg, rep):
    phi_0 = phi
    
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
        
        #calculate trig and common terms now so they aren't repeatedly done later
        #this also makes formulae more readable
        #commented number is number of times term appears in formulae
        sin_chi = sin(chi) #111111111
        cos_chi = cos(chi) #111111
        sin_phi = sin(phi) #111111
        cos_phi = cos(phi) #111111
        sin_theta = sin(theta) #11111111111111111
        cos_theta = cos(theta) #111111111111111
        tan_theta = tan(theta) #111
        sin_theta2 = sin_theta**2 #11111111111111
        
        ##PRECALCULATED TRIG
        #spherical vectors
        e_r = sin_theta*cos_phi*x + sin_theta*sin_phi*y + cos_theta*z
        e_theta = cos_theta*cos_phi*x + cos_theta*sin_phi*y - sin_theta*z
        e_phi = -sin_phi*x + cos_phi*y
        
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
        
        '''
        ##ORIGINAL
        #spherical vectors
        e_r = sin(theta)*cos(phi)*x + sin(theta)*sin(phi)*y + cos(theta)*z
        e_theta = cos(theta)*cos(phi)*x + cos(theta)*sin(phi)*y - sin(theta)*z
        e_phi = -sin(phi)*x + cos(phi)*y
        
        #pure dipole magnetic field
        B_d = (B_0/eta**3)*( cos(theta)*e_r + 0.5*sin(theta)*e_theta )
        
        #first correction due to charge flow along open field lines
        B_r_1 = (3/2)*(B_0/eta**3)*theta_0*xi*(1-theta/tan(theta))*sin(chi)*sin(phi)
        B_theta_1 = -(3/2)*(B_0/eta**3)*theta_0*xi*(theta - sin(theta)*cos(theta)/(sin(theta))**2)*sin(chi)*sin(phi)
        B_phi_1 = -(B_0/eta**3)*( (1-k)*cos(chi)*sin(theta) + (3/2)*theta_0*xi*( (sin(theta)**3 + sin(theta) - theta*cos(theta))/sin(theta)**2 )*sin(chi)*cos(phi) )
        B_1 = (eta/eta_lc)*(B_r_1*e_r + B_theta_1*e_theta + B_phi_1*e_phi)
        
        #second correction due to rotation
        B_r_2 = ((eta/eta_lc)**2)*(B_0/eta**3)cos(theta)*( cos(chi)*( cos(chi)*sin(theta)**2 + 2*lambda_0*(1-k) ) + (sin(chi)**2)*(1 - 0.5*sin(theta)**2) - 2*( cos(chi)*sin(theta)*cos(theta) + (3/2)*lambda_0*theta_0*xi*(theta/(sin(theta)*cos(theta)) - 1) )*sin(chi)*cos(phi) - 0.5*(sin(chi)**2)*(sin(theta)**2)*cos(2*phi) )
        B_theta_2 = -((eta/eta_lc)**2)*(B_0/eta**3)*sin(theta)*( ( 0.25*cos(chi)*(sin(theta)**2) + lambda_0*(1-k) )*cos(chi) + 0.5*(sin(chi)**2)*(1 - 0.25*(sin(theta)**2)) + ( ( (0.25/tan(theta))*cos(2*theta) + (1 - cos(theta))/(sin(theta)**3) )*cos(chi) - (3/2)*lambda_0*theta_0*xi*(3*(1 - theta/tan(theta))/(sin(theta)**2) - 1) )*sin(chi)*cos(phi) - (5/8)*( (1/(sin(theta)**2))*( (1-cos(theta))/(sin(theta)**2) - 0.5 ) + 0.2*(sin(theta)**2) )*(sin(theta)**2)*cos(2*theta) )
        B_phi_2 = ((eta/eta_lc)**2)*(B_0/(eta**3))*( ( (0.25 + (1-cos(theta))/(sin(theta)**2) )*cos(chi) + (3/2)*lambda_0*theta_0*xi*((3*sin(theta)*cos(theta) - theta*(3-2*(sin(theta)**2)))/(sin(theta)**2)) )*sin(chi)*sin(phi) - (5/8)*( (1-cos(theta))/(sin(theta)**3) - (2*sin(theta))**(-1) )*(sin(theta)**2)*sin(2*phi) )
        B_2 = B_r_2*e_r + B_theta_2*e_theta + B_phi_2*e_phi
        
        #simplified case for aligned rotator
        B_r_s = (B_0/eta**3)*cos(theta)*( 1 + ((eta/eta_lc)**2)*( sin(theta)**2 - 2*k*(1-k) ) )
        B_theta_s = (0.5)*(B_0/eta**3)*sin(theta)*( 1 - (0.5)*((eta/eta_lc)**2)*( sin(theta)*2 - 4*k*(1-k) ) )
        B_phi_s = (-B_0/eta**3)*(1-k)*(eta/eta_lc)*sin(theta)
        B_s = B_r_s*e_r + B_theta_s*e_theta + B_phi_s*e_phi
        '''
        
        #plot magnetic field
        B = 10*(B_d + B_1 + B_2)/mag(B_d + B_1 + B_2)
        #B = 10*B_s/mag(B_s)
        arrow(pos=position, axis=B, color=color.white) #magnetosphere
                
        #current position
        position = position + B

        #print feedback
        '''
        print("\n x, z, theta:", position.x, position.z, theta*180/pi)
        print("position:", position)
        print("mag vector:", B)
        print("mag field strength:", mag(B))
        '''

for j in range(1):
    phi = pi*j/10
    #drawMagLine(phi, 25, 1000)
    drawMagLine(phi, 20, 42)
    #drawMagLine(phi, 15, 75)
    drawMagLine(phi, 10, 165)
    #drawMagLine(phi, 5, 2000)
