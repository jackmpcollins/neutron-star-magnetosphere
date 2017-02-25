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


def main():
    ns = specifyStarParameters()
    #ns.setOtherk()
    startTime = time()

    #draw star, rotational axis, light cylinder, magnetosphere
    #ns.draw(5, 10, 15)
    
    #need to determine polar cap so know what values to use
    theta_max = 5
    emissions = [[0]*360 for i in range(180)]
    for phi in range(0,360,1):
        print("phi: ", phi)
        for theta in np.arange(1,theta_max, 0.1):
            (phi_emis, theta_emis) = Fieldline(ns, phi, theta).emissionDirection
            emissions[int(round(theta_emis))][int(round(phi_emis))] += 1

    plotEmissions(emissions)
    #plotEmissions()
    #output instructions for manipulation
    print("To zoom in/out, hold down the scroll wheel and move the mouse forward/backward.")
    print("To move star, hold down the right mouse button to grab, and move the mouse.")
    #ask user if want to do a different pulsar
    print("Time elapsed:", time() - startTime)
    return 0

main()