"""Programmed by Binay Limbu as part of the third year project
University of Southampton, 2013
"""

try:
    vel = float(raw_input("please enter the velocity of the object in m/s >>>"))
except ValueError:
    pass
try:
    reentry_angle = float(raw_input("Please enter the reentry angle in degrees >>>"))
except ValueError:
    pass
try:
    drop_height = float(raw_input("Please enter the height in meters >>>"))
except ValueError:
    pass
try:
    rad = float(raw_input("please enter the radius in meters >>>"))
except ValueError:
    pass
try:
    mat_density = float(raw_input("Please enter the density of material in kg/m^3 >>>"))
except ValueError:
    pass
try:
    c_d = float(raw_input("Please enter the drag coefficient >>>"))
except ValueError:
    pass
try:
    e_h_o_a = float(raw_input("please enter the effective heat of ablation in J/kg >>>"))
except ValueError:
    pass
try:
    h_l_t_b_a = float(raw_input("please enter the heat load to begin ablation in W/m^2>>>"))
except ValueError:
    pass
try:
    t_final = float(raw_input("please enter the number of seconds you want it to run >>>"))
except ValueError:
    pass
try:
    n = int(raw_input("Please enter the number of iterations(integer). >>>"))
except ValueError:
    pass

import math
import numpy
import pylab

def solver( t_final,vel ,drop_height,n):
    """This function will solve a set of first order ordinary differential equations."""
    #global gamma
    time_step = (t_final - 0)/float(n)  #time steps
    v = vel
    h = drop_height
    d = 0
    t = 0
    r = rad
    angle = math.radians(reentry_angle)
    y_speed = [v*math.sin(angle)]
    x_speed = [v*math.cos(angle)]
    velocity = [v]
    height = [h]
    distance = [d]
    theta = [angle]
    time = [0]
    rad_map = [r]
    heat_rate_map = [0]
    recession_rate_map = [0]
    bcplot = [0]
    for i in range(n):
        print i
        dense = density(h)
        heat_rate = convective_heat(dense,v,r)

        print 'heat rate',heat_rate
        recession_rate = 0
        if heat_rate > h_l_t_b_a:
            recession_rate = ablate(heat_rate, mat_density,e_h_o_a)
            print 'recession rate', recession_rate
            r -= (recession_rate)*time_step
            print r
            if r <= 0:
                break
        area = area_sphere(r)
        print 'area' ,area
        kg = mass(r, mat_density)
        print 'mass' ,kg
        b_c = ballistic(area, kg, c_d)

        print 'ballist', b_c
        print 'density = %f'%dense
        accel = vel_eqn(dense,v,angle,b_c)
        v += time_step*accel
        angle += time_step*angle_eqn(v,angle)
        vertical_speed = height_eqn(v,angle)
        h += time_step*vertical_speed
        horizontal_speed = distance_eqn(v,angle)
        d += time_step*horizontal_speed
        t += time_step
        print 'angle = %f,velocity = %f , distance = %f, acceleration = %f, horizontal speed = %f'%(angle, v, d,accel, horizontal_speed)
        print 'height %f'%h
        print 'Hello'
        bcplot.append(b_c)
        heat_rate_map.append(heat_rate)
        recession_rate_map.append(recession_rate)
        rad_map.append(r)
        time.append(t)
        height.append(h)
        distance.append(d)
        velocity.append(v)
        theta.append(angle)
        y_speed.append(vertical_speed)
        x_speed.append(horizontal_speed)
        if h <= 0:
            break
        #else:
            #continue
    a = numpy.asarray([height,distance,velocity,y_speed,x_speed,heat_rate_map,recession_rate_map, rad_map, bcplot,time])
    numpy.savetxt("Best_case.csv",a, delimiter = ",")
    bc_graph(height,bcplot)
    graph_q_dot(height, heat_rate_map)
    graph_recession_rate(height, recession_rate_map)
    height_distance_graph(height, distance)
    x_graphs(time,distance,x_speed)
    y_graphs(time,height, y_speed)
    vel_graph(time,velocity)
    theta_graph(time,theta)
    radius_graph(time,rad_map)
    print 'finish'

def bc_graph(height, bcplot):
    """plots the graph of bccoeff against height"""
    pylab.plot(bcplot,height,label = 'ballistic coefficient')
    pylab.grid()
    pylab.xlabel('ballistic coeff')
    pylab.ylabel('height')
    pylab.show()

def graph_q_dot(he, heats):
    """maps heat rate against height"""
    pylab.plot(heats, he, label = 'heat rate')
    pylab.grid()
    pylab.xlabel('heat rate')
    pylab.ylabel('height')
    pylab.show()

def graph_recession_rate(height, recession_rate_map):
    """maps recession rate against height."""
    pylab.plot(recession_rate_map, height , label = 'recession rate')
    pylab.grid()
    pylab.xlabel('recession rate')
    pylab.ylabel('height')
    pylab.show()

def radius_graph(time,rad_map):
    """prints the radius of the sphere with time."""
    pylab.plot(time, rad_map, label = 'radius size')
    pylab.grid()
    pylab.xlabel('time')
    pylab.ylabel('radius')
    pylab.show()

def height_distance_graph(height, distance):
    """This prints the height against distance graph."""
    pylab.plot(distance, height, label = 'reentry trajectory')
    pylab.grid()
    pylab.xlabel('distance')
    pylab.ylabel('height')
    pylab.show()


def theta_graph(time, theta):
    """This prints the graph for theta against time"""
    pylab.plot(time, theta, label = 'reentry angle')
    pylab.grid()
    pylab.xlabel('time')
    pylab.ylabel('gamma')
    pylab.show()


def vel_graph(time,speed):
    """This prints the graph for velocity against time."""
    pylab.plot(time, speed, label = 'velocity')
    pylab.grid()
    pylab.xlabel('time')
    pylab.ylabel('velocity')
    pylab.show()


def x_graphs(time,distance,speed_hor):
    """Prints graph relating to the x components."""
    pylab.plot(time,distance,label = 'distance')
    #pylab.legend(('distance','speed'))
    pylab.grid()
    pylab.xlabel('time')
    pylab.ylabel('x - distance')
    pylab.show()
    pylab.plot(time,speed_hor)
    pylab.grid()
    pylab.xlabel('time')
    pylab.ylabel('x-speed')
    pylab.show()


def y_graphs(time,height,speed_ver):
    """Prints graph relating to y components."""
    pylab.plot(time,height,label = 'height')
    #pylab.legend(('height','speed'))
    pylab.grid()
    pylab.xlabel('time')
    pylab.ylabel('height')
    pylab.show()
    pylab.plot(time, speed_ver, label = 'speed')
    pylab.grid()
    pylab.xlabel('time')
    pylab.ylabel('y-speed')
    pylab.show()



def vel_eqn(dense, v, angle,b_c):
    """This is used to calculate the slope or the acceleration of the equation."""
    #global b_coeff
    return -9.8*math.sin(angle) - (0.5*dense*v**2)/b_c

def angle_eqn(v,angle):
    """:must return the angle in radians,
:the angle must be negative."""
    return (-9.8*math.cos(angle)/v)


def height_eqn(v,angle):
    """returns the vertical velocity."""
    return v*math.sin(angle)

def distance_eqn(v,angle):
    """returns the horizontal speed at the time."""
    return v*math.cos(angle)



def density(height):
    """This is the function that will be called by the other equations to compute the
density."""
    if height <=  11000:
        T = 15.04 - (0.00649*height)
        P = 101.29 *((T + 273.1)/288.08)**(5.256)
        rho = P/(.2869*(T + 273.1))
    elif height >11000 and height <= 25000:
        T = -56.46
        P = 22.65 * math.exp(1.73 - 0.000157*height)
        rho = P/(.2869*(T + 273.1))
    elif height>25000:
        T = - 131.21 + (.00299*height)
        P = 2.488 *((T+273.1)/216.6)**(-11.388)
        rho = P/(.2869 * (T +273.1))
    return rho

def area_sphere(r):
    """returns the area!"""
    return math.pi*(r**2)

def convective_heat(atmos_dense, vel_body, rad_body):
    """needs density, velocity and radius to find the convective heating."""
    return 0.000183*(vel_body**3)*math.sqrt(atmos_dense/rad_body)

def ballistic(a,kilo,drag_c):
    """calculates new B_C."""
    #a = area(rad)
    return kilo/(drag_c*a)
         
def mass(rad_body,mat_dens):
    """finds the new mass"""
    return (4/3)*math.pi*(rad_body**3)*mat_dens

def ablate(q, mat_dens, h_o_a):
    """This caluculates the recession rate."""
    return q/(mat_dens*h_o_a)

if raw_input("Do this ? (y/n) >>> ") != 'n':
    solver(t_final, vel, drop_height,n)
