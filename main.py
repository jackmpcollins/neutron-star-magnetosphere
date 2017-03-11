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
def plotEmissions(emissions):
    em = gaussian_filter(np.asarray(emissions), sigma=3)
    plt.pcolor(em)
    plt.xlabel("phi (degrees)")
    plt.ylabel("theta (degrees)")
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


#allows for different blurring/scattering
#Does this just do the same as blurring the image after making it the usual way?
#this is redundant - Used gaussian_filter from scipy instead
def addWithBlur(direc, data):
    (phi_emis, theta_emis) = direc
    blurring = [[1,2,1],[2,4,2],[1,2,1]] #Approximate Gaussian blur
    row = int(round(theta_emis))
    print("row: ", row)
    col = int(round(phi_emis) % 360)
    print("col: ", col)
    for i in range(3):
        for j in range(3):
            print("i,j: ", i , j)
            print("blurring i j: ", blurring[i][j])
            data[row-1+i][col-1+j] = blurring[i][j]
    return(data)


#finds the staring position of the last open field line for a particular value of phi
def findPC(star, phi, decimalPlaces):
    theta = 2.0
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

'''
def findPC(star, phi, decimalPlaces, theta_guess):
    while stepSize > 0.01:
        line = Fieldline(star, phi, theta_guess)
        if line.isOpen:
            theta_guess += stepSize
'''

def main():
    ns = specifyStarParameters()
    #ns.setOtherk()
    startTime = time()

    emissions = [[0.0]*360 for i in range(180)]
    for phi in np.arange(0,360,1):
        print("phi: ", phi)
        theta = 2.0
        line = Fieldline(ns, phi, theta, False)
        #for theta in np.arange(1,theta_max, 0.1):
        while line.isOpen:
            if line.alreadyEmitted:
                (phi_emis, theta_emis) = line.emissionDirection
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
def main2():
    ns = specifyStarParameters()
    #ns.setOtherk()
    startTime = time()

    #draw star, rotational axis, light cylinder, magnetosphere
    ns.draw(5, 10, 15)
    ns.animate()
    
    #output instructions for manipulation
    print("To zoom in/out, hold down the scroll wheel and move the mouse forward/backward.")
    print("To move star, hold down the right mouse button to grab, and move the mouse.")

    print("Time elapsed:", time() - startTime)
    return 0


#this plots emissions using square angular coordinates to get a more equal distribution
def main3():
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
                (phi_emis, theta_emis) = line.emissionDirection
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


#This plots emissions from the last open field lines and lines x degrees towards the magnetic pole from these
#TOO SLOW
def main4():
    ns = specifyStarParameters()
    #ns.setOtherk()
    startTime = time()

    print("ns.theta_0: ", ns.theta_0)

    emissions = [[0.0]*360 for i in range(180)]
    for phi in np.arange(0,360,1):
        print("phi: ", phi)

        theta_max = findPC(ns, phi, 2)
        print("theta_max: ", theta_max)

        for theta in np.arange(theta_max-2,theta_max, 0.1):
            line = Fieldline(ns, phi, theta, False)
            (phi_emis, theta_emis) = line.emissionDirection
            emissions[int(round(theta_emis))][int(round(phi_emis)) % 360] += 1.0

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


#Plots emissions only from within a band of given thickness at the end of the polar cap
def main5():
    ns = specifyStarParameters()
    #ns.setOtherk()
    startTime = time()

    emissionBandThickness = 1 #degrees
    theta_min = min(findPC(ns, 0, 2), findPC(ns, 45, 2), findPC(ns, 90, 2)) - emissionBandThickness - 0.1

    stepSize = 0.05

    emissions = [[0.0]*360 for i in range(180)]
    for phi in np.arange(0,360,1):
        print("phi: ", phi)
        theta = theta_min
        line = Fieldline(ns, phi, theta, False)

        listOfEmissionsFromThisPhi = []

        while line.isOpen:
            if line.alreadyEmitted:
                #(phi_emis, theta_emis) = line.emissionDirection
                listOfEmissionsFromThisPhi.append(line.emissionDirection)
            theta += stepSize
            line = Fieldline(ns, phi, theta, False)


        for coords in listOfEmissionsFromThisPhi[-int(emissionBandThickness/stepSize):]: #we only want emissions from within 'emissionBandThickness' of last open field lines
            (phi_emis, theta_emis) = coords
            emissions[int(round(theta_emis))][int(round(phi_emis)) % 360] += 1.0


        print("theta PC (degrees): ", theta)

    plotEmissions(emissions)
    saveEmissions(emissions)

    print("Time elapsed:", time() - startTime)

    obsAng = int(raw_input("Enter the observers angle (or 0 to quit): "))
    while obsAng > 0:
        plt.plot(emissions[obsAng])
        plt.xlabel("phi (degrees)")
        plt.ylabel("relative emission")
        plt.show()
        obsAng = int(raw_input("Enter the observers angle (or 0 to quit): "))

    return 0

main5()