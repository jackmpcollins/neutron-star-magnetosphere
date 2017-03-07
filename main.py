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
    plt.pcolor(np.asarray(emissions))
    plt.xlabel("phi (degrees)")
    plt.ylabel("theta (degrees)")
    plt.show()
    return 0

#allows for different blurring/scattering
#Does this just do the same as blurring the image after making it the usual way?
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

def testBlur():
    a = [[0]*10 for i in range(5)]
    direc = (4,3)
    print(a)
    print(addWithBlur(direc, a))
    return 0

def main():
    ns = specifyStarParameters()
    #ns.setOtherk()
    startTime = time()

    theta_max = 1.1 #remove
    emissions = [[0]*360 for i in range(180)]
    for phi in np.arange(0,360,1):
        print("phi: ", phi)
        theta = 2.0
        line = Fieldline(ns, phi, theta, False)
        #for theta in np.arange(1,theta_max, 0.1):
        while line.isOpen:
            if line.alreadyEmitted:
                (phi_emis, theta_emis) = line.emissionDirection
                emissions[int(round(theta_emis))][int(round(phi_emis)) % 360] += 1
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


main()