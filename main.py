## Jack Collins, Final Year Project, National University of Ireland Galway, 2017
## written for Python2.7

## This file defines calculates the emission an observer at various angles to the
## rotation axis would see

from classStar import *
from classFieldline import *
import matplotlib.pyplot as plt
from visual import *
from time import time
import numpy as np
from scipy.ndimage.filters import gaussian_filter

#asks user which star to model
def specifyStarParameters():
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
        name = raw_input("If you would like to save this star enter a name (to skip hit Enter): ")
        if name != "":
            f = open("savedStars.txt", 'w')
            f.write(",".join([name,R,P,chi]))
    else:
        (R, P, chi) = [float(item) for item in starInfo[starName]]
    star = Star(R, P, chi)
    return star


#plots graph of emissions
def plotEmissions(emissions, chi=-1, rad=-1):
    em = gaussian_filter(np.asarray(emissions), sigma=3)
    plt.figure(figsize=(10,5))
    plt.pcolor(em, vmin=0, vmax=1.8)
    plt.axis([360,0,180,0])
    #show color scale
    plt.xlabel(r'Phase $\phi$ (degrees)', {'size':'15'})
    plt.ylabel("Observer Viewing Angle (degrees)", {'size':'15'})
    if chi >= 0:
        plt.text(360, 175, r'$\chi=$' + str(chi)+r'$^{\circ}$', size=20, color='white')
    if rad >= 0:
        plt.text(360, 175, str(rad)+r'$R_{LC}$', size=20, color='white')
    plt.colorbar()
    plt.show()
    return 0


#saves matrix of emissions WITHOUT BLUR to file
def saveEmissions(emissions):
    fname = raw_input("Enter a name for the emissions file: ")
    if fname == "":
        print("Not saved")
        return 0
    f = open("savedEmissions/"+fname+".txt", 'w')
    for line in emissions:
        f.write(str(line)[1:-1] + "\n")
    f.close()
    print("File saved")
    return 0


#saves polar cap angles
def savePC(angles, filename):
    f = open("savedPolarCapPoints/"+filename+".txt", 'w')
    f.write(str(angles)[1:-1])
    f.close()

def readPC(filename):
    f = open("savedPolarCapPoints/"+filename+".txt", "r")
    points = f.read()[1:-1].split("), (")
    f.close()

    thetas = []
    for point in points:
        thetas.append(float(point.split(", ")[1]))
    return(thetas)


#plots the shape of the polar cap for a given magnetic inclination angle
def plotPCPoints(angle):
    file = "fifthBandAngle" + str(angle) + "PolarCapPoints"
    thetas = readPC(file)

    #create x-y points to plot birds eye view
    x = []
    y = []
    for i in range(360):
        phi = deg2rad(i - 90) # rotate plot by 90 so phi = 0 is in negative y direction - easier to see change
        x.append(thetas[i]*cos(phi))
        y.append(thetas[i]*sin(phi))

    plt.figure(figsize=(16,16))
    plt.scatter(x,y)
    plt.scatter(0,0, marker='x', color='black', s=1000)
    plt.arrow(0, 0, 0, -3, head_width=0.5, head_length=1, fc='k', ec='k')
    plt.text(1, 0.5, "velocity", size=40)

    plt.arrow(0, 0, 3, 0, head_width=0.5, head_length=1, fc='k', ec='k')
    plt.text(0, -2.5, r'$\phi = 0^{\circ}$', size=40)

    plt.text(-10, 10, r'$\chi = $'+str(angle)+r'$^{\circ}$', size=80)
    plt.grid(True)
    plt.axis([-12,12,-12,12])
    plt.show()
    return 0


#finds the staring position of the last open field line for a particular value of phi
#probably should start with a smaller theta, or take a minimum theta guess as a parameter
#probably won't work for crab
def findPC(star, phi, decimalPlaces):
    theta = 0.1
    line = Fieldline(star, phi, theta)
    for i in range(decimalPlaces+1):
        stepSize = 10**(-i)
        while line.isOpen:
            theta += stepSize
            line = Fieldline(star, phi, theta)
        theta -= stepSize
        line = Fieldline(star, phi, theta)
    return theta #degrees


#plots theta vs phi for the polar cap
def plotPC():
    ns = specifyStarParameters()
    #ns.setOtherk()
    startTime = time()

    PCAngles = []
    phi_range = np.arange(0,360,30)
    for phi in phi_range:
        print("phi: ", phi)
        PCAngles.append(findPC(ns, phi, 1))
    
    print("Time elapsed:", time() - startTime)

    print(PCAngles)
    plt.plot(PCAngles)
    plt.show()
    return PCAngles


#outputs the emissions view of a particular observer angle
def observerView(emissions, obsAng):
    em = gaussian_filter(np.asarray(emissions), sigma=3)
    observedEmissions = np.fliplr(em)[obsAng] #flipping because star rotates oposite direction to increasing phi
    plt.plot(observedEmissions)
    plt.xlabel("phase (degrees)")
    plt.ylabel("relative emission")
    plt.text(5, 1.55, str(obsAng)+r'$^{\circ}$', size=20)
    plt.axis([0,360,0,1.7])
    plt.show()
    return 0


#plots emissions from a band with specified inner radius to the edge of the polar cap
def innerCircleToPolarCap():
    ns = specifyStarParameters()
    #ns.setOtherk()
    startTime = time()

    emissions = [[0.0]*360 for i in range(180)]
    for phi in np.arange(0,360,1):
        print("phi: ", phi)
        theta = 2.0 #inner radius
        line = Fieldline(ns, phi, theta, False)
        #for theta in np.arange(1,theta_max, 0.1):
        while line.isOpen:
            (phi_emis, theta_emis) = line.emissionDirections[0]
            emissions[int(round(theta_emis))][int(round(phi_emis)) % 360] += 1.0
            theta += 0.5
            line = Fieldline(ns, phi, theta, False)
        print("theta PC (degrees): ", theta)

    plotEmissions(emissions)

    print("Time elapsed:", time() - startTime)

    obsAng = int(raw_input("Enter the observers angle (or 0 to quit): "))
    while obsAng > 0:
        plt.plot(emissions[obsAng])
        plt.xlabel("phi (degrees)")
        plt.ylabel("relative emission")
        plt.show()
        obsAng = int(raw_input("Enter the observers angle (or 0 to quit): "))

    return 0


#draws a star with field lines at theta = 5, 10, 15 degrees
def visualise3D():
    ns = specifyStarParameters()
    #ns.setOtherk()
    startTime = time()

    #draw star, rotational axis, light cylinder, magnetosphere
    ns.draw(3, 5, 10)
    ns.animate()
    
    #output instructions for manipulation
    print("To zoom in/out, hold down the scroll wheel and move the mouse forward/backward.")
    print("To move star, hold down the right mouse button to grab, and move the mouse.")

    print("Time elapsed:", time() - startTime)
    return 0


#this plots emissions using square angular coordinates to get a more equal distribution
def plotUsingNewCoords():
    ns = specifyStarParameters()
    #ns.setOtherk()
    startTime = time()

    theta_max = max(findPC(ns, 0, 2), findPC(ns, 90, 2))
    increment_size= 0.2
    limit = deg2rad(0.01)

    emissions = [[0.0]*360 for i in range(180)]
    for A in np.arange(-theta_max,theta_max,increment_size):
        print("A ", A)
        Adeg = A
        A = deg2rad(A)
        x = sin(A)
        for B in np.arange(-theta_max,theta_max,increment_size):
            #print("A, B: ", Adeg, B)
            B = deg2rad(B)
            y = sin(B)
            z = cos(A)*cos(B)
            position = vector(x,y,z)
            #print("A,B ", A,B)

            if abs(A) < limit and abs(B) < limit:
                #don't try to create a field line too close to the magnetic pole
                print("here!")
                break

            (x, y, z) = position/mag(position)
            theta = rad2deg(acos(z))
            phi = rad2deg(atan2(y, x))
            #print(position)
            #print(ns, phi, theta, False, 100)
            #rate(5)
            
            line = Fieldline(ns, phi, theta, False)
            if line.isOpen:
                print(phi, theta)
                (phi_emis, theta_emis) = line.emissionDirections[0]
                emissions[int(round(theta_emis))][int(round(phi_emis)) % 360] += 1.0

    print("Time elapsed:", time() - startTime)

    plotEmissions(emissions)
    saveEmissions(emissions)

    obsAng = int(raw_input("Enter the observers angle (or 0 to quit): "))
    while obsAng > 0:
        plt.plot(emissions[obsAng])
        plt.xlabel("phi (degrees)")
        plt.ylabel("relative emission")
        plt.show()
        obsAng = int(raw_input("Enter the observers angle (or 0 to quit): "))

    return 0


#One function to rule them all, one function to find them, One function to bring them all and in the darkness bind them.
def emissionsFromBand(chi, emissionRadii):
    ns = Star(10, 0.03, chi)
    #ns.setOtherk()
    startTime = time()

    emissionBandThickness = 0.2 #degrees
    stepSize = 0.02

    #find (roughly) biggest incircle of polar cap and come in slightly more than required for calculations
    polarCapDistances = []
    for angle in [0, 45, 90, 135, 180, 200, 225, 270, 315, 360]:
        polarCapDistances.append(findPC(ns, angle, 2))
    theta_min = min(polarCapDistances) - emissionBandThickness - 5*stepSize

    emissions = [ [[0.0]*360 for i in range(180)] for j in range(len(emissionRadii))]
    polarCapPoints = []

    for phi in np.arange(0,360,1):
        print("phi: ", phi)
        theta = theta_min
        line = Fieldline(ns, phi, theta, False, emissionRadii)

        listOfEmittingFieldlinesFromThisPhi = []

        while line.isOpen:
            listOfEmittingFieldlinesFromThisPhi.append(line.emissionDirections)
            theta += stepSize
            line = Fieldline(ns, phi, theta, False, emissionRadii)

        print("This should be the same number always: ", len(listOfEmittingFieldlinesFromThisPhi[-int(round(emissionBandThickness/stepSize)):]))
        for individualFieldline in listOfEmittingFieldlinesFromThisPhi[-int(round(emissionBandThickness/stepSize)):]: #we only want emissions from within 'emissionBandThickness' of last open field lines
            for i in range(len(individualFieldline)):
                (phi_emis, theta_emis) = individualFieldline[i]
                emissions[i][int(round(theta_emis))][int(round(phi_emis)) % 360] += 1.0


        #print("theta PC (degrees): ", theta)
        polarCapPoints.append((phi,theta))

    #savePC(polarCapPoints, "fifthBandAngle"+str(chi)+"PolarCapPoints")

    for i in range(len(emissionRadii)):
        fname = "distChangeRadius"+str(emissionRadii[i])+"Angle"+str(chi)
        f = open("savedEmissions/"+fname+".txt", 'w')
        for line in emissions[i]:
            f.write(str(line)[1:-1] + "\n")
        f.close()
        print("File saved")


    print("Time elapsed:", time() - startTime)

    return 0


#calculates and saves the emissions from stars with multiple of 10 degree magnetic inclination
#can take long to run
def batchCalculate():
    emissionRadius = int(raw_input("Enter emission distance to light cylinder (% LC): "))
    startTime = time()

    for chi in range(0, 100, 10):
        print("chi: ", chi)
        try:
            emissionsFromBand(chi, [emissionRadius])
        except:
            print("It failed oh no .. KEEP GOING")

    print("TOTAL TIME TAKEN = ", time() - startTime)


#Reads from emissions file and returns array of emissions
def readFromFile(filename):
    f = open("savedEmissions/"+filename, 'r')
    emissions = []
    for line in f:
        emissions.append([float(item) for item in line[:-1].split(", ")])
    return emissions

'''
#sample code for plotting from saved files
name1 = "fifthBandRadius0.5Angle"
name2 = "distChangeRadius0.5Angle60.txt"
for i in range(1, 10):
    emissions = readFromFile("distChangeRadius"+str(i/10.0)+"Angle45.txt")
    plotEmissions(emissions, rad=i/10.0)
'''

'''
#sample code for calculating emission from range of distances for specid chi
chi = int(raw_input("Enter chi: "))
emRads = [item/10.0 for item in [1,2,3,4,5,6,7,8,9]]
print(emRads)
emissionsFromBand(chi, emRads)
'''

if __name__ == "__main__":
    # your code goes here
    print("edit main function to your own specifications")



'''NOTES
SAVE THESE USING 3D ARRAY - done
WHEN SAVING EMISSIONS, SPLIT INTO SEVERAL FILES SO THAT PLOTTING STILL WORKS - done
MAKE FUNCTION THAT RUNS FOR DIFFERENT CHI AND RUN THESE IN PARALLEL - done

PLOT ALL FIELD LINES IN POLAR CAP - avoided? Avoided woo!
CHECK MORE POINTS FOR MIN POLAR CAP THETA - done


'''