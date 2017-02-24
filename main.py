## Jack Collins, Final Year Project, National University of Ireland Galway, 2017
## written for Python2.7

from visual import *
from time import time

from classStar import Star
from functionsMagneticField import drawMagnetosphere, plotEmissions


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

def main():
    ns = specifyStarParameters()
    #ns.setOtherk()
    startTime = time()
    #draw star, rotational axis, light cylinder
    ns.draw()
    #calculate points using values given
    #plot field
    phiAngle = float(raw_input("phi angle: "))
    drawMagnetosphere(ns, phiAngle)
    #output instructions for manipulation
    print("To zoom in/out, hold down the scroll wheel and move the mouse forward/backward.")
    print("To move star, hold down the right mouse button to grab, and move the mouse.")
    #ask user if want to do a different pulsar
    plotEmissions()
    print("Time elapsed:", time() - startTime)
    return 0

main()
