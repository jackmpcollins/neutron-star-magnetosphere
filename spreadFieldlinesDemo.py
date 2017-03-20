## This program (hopefully) plots an equally distributed set of points on the surface of a sphere

from visual import *

sphere(radius = 1)
sphere(pos=vector(0,0,1), radius=0.006, color=color.red)

def new():
	for i in range(50):
		A = i*pi/100
		for j in range(-50,50):
			B = j*pi/100
			x = sin(A)
			y = sin(B)
			z = cos(A)*cos(B)
			position = vector(x,y,z)
			sphere(radius=0.005, pos=position/mag(position), color=color.blue)

def old():
	for i in range(50):
		theta = i*pi/100
		for j in range(-50,50):
			phi = j*pi/100
			x = sin(theta)*cos(phi)
			y = sin(theta)*sin(phi)
			z = cos(theta)
			position = vector(-x,-y,z)
			sphere(radius=0.005, pos=position/mag(position), color=color.black)

new()
old()
			
