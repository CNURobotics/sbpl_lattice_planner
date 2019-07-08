"""
This code generates the motion primitive generators.
Select the relevant subset for generating the full motion primitive file using
genmprim_extended.py.

"""

import sys#, types, time, random, os
import numpy as np
import math
#from mprim_generator_utils import *
#from mprim_generator_generators import *
from mprim_generator_primitives import *
#from mprim_generator_io import *

# if available import pylab (from matlibplot)
visualize_plt = False
try:
    import matplotlib.pylab as plt
    visualize_plt = True
except ImportError:
    pass

n_angles = 16
#resolution = 0.100
resolution = 0.050
use_uniform_angles = False

min_radius = 0.230*0.25 # Will allow some negative motion of inside wheel of Turtlebot

min_turning_radius_m = min_radius

grid_max = 8
gridx = np.array(list(range(-grid_max,grid_max+1,1)))
gridy = np.array(list(range(-grid_max,grid_max+1,1)))
gridx,gridy = np.meshgrid(gridx,gridy)

grid = np.vstack((gridx.flatten(),gridy.flatten()))
print(grid)


if use_uniform_angles:
    print "Using uniform angle spacing (but straight motions don't align with grid)"
    # Use uniform angles instead of
    ang_inc = 2*math.pi/n_angles
    angles = [ang_inc*i for i in range(n_angles//4)]
else:
    # Use non-uniform angles tp match the grid
    print "Using non-uniform angles that match the grid"
    angles = [math.atan2(y,x) for x,y in [(1,0), (2,1), (2,2), (1,2)] ]

print angles
degrees=[angle*180.0/math.pi for angle in angles]
print "angle (degrees)=", degrees

angids = list(range(4)) # Angle ID for generator
deltas = list(range(1,5)) # Angles we will move to from starting angle

def validate_grid(args, angleind, primind,
                  ex, ey, delta, resolution,
                  uv, nv, costmult=1):

    basemprimendpts_c = np.array((ex, ey, delta, costmult),dtype=int)
    baseendpose_c = np.array((ex, ey, delta, costmult))
    endtheta_c    = (angid+basemprimendpts_c[2])%n_angles
    endpose_c     = np.array((ex, ey, endtheta_c))

    startpt = np.array((0., 0., angles[angid]))
    endpt   = np.array((endpose_c[0]*resolution,
                        endpose_c[1]*resolution,
                        args['angles'][endtheta_c] ) )
    currentangle = args['angles'][angleind]
    endangle     = args['angles'][endtheta_c]

    prim_def_primitive={}

    #print 40*"="
    #print "base=",basemprimendpts_c
    if not generate_first_quadrant_primitive(angleind, primind,
                                             currentangle, basemprimendpts_c,
                                             args,
                                             prim_def_primitive ):
        #print "Failed to define valid motion primitive for first quadrant", basemprimendpts_c
        # print "Invalid:  {} {} {} {} #{} {} {}".format(
        #     basemprimendpts_c[0],basemprimendpts_c[1],
        #     basemprimendpts_c[2],basemprimendpts_c[3],
        #     prim_def_primitive['style'],
        #     prim_def_primitive['line length'],
        #     prim_def_primitive['turning_radius'])
        return False
    else:

        line_seg = prim_def_primitive['line length']
        radius   = prim_def_primitive['turning_radius']
        del_angle = shortest_angle_dist(endangle,currentangle)
        cost      = prim_def_primitive['actioncost']

        #print line_seg, radius
        if math.isinf(radius) or line_seg < 2.25*resolution:
            # Otherwise we can definitely get there by straight move then turn
            print "{:3d},{:3d},{:3d},{:3d},{:2d}, #{} {:.3f} {:.3f} {}".format(
                angleind,
                basemprimendpts_c[0],basemprimendpts_c[1],
                basemprimendpts_c[2],basemprimendpts_c[3],
                prim_def_primitive['style'],
                line_seg,
                radius,
                cost
                )
        #else:
        #    print "too long line!"
        return True

full_angles =angles[:]
for quad in [math.pi/2, math.pi, 3*math.pi/2]:
    for ang in angles:
        full_angles.append(ang+quad)

args = dict()
args['resolution'] = resolution
args['numberofanglesperquadrant'] = len(angles)
args['numberofprimsperangle'] = 99
args['numberOfSamples'] = 10
args['minTurnRadius'] = min_turning_radius_m
args['check_arc_line_arc_generator'] = True
args['wheelbase'] = 0.230 # Turtlebot
args['cost_conversion_factor']  = 1000

args['angles'] = full_angles

for key, value in args.iteritems():
    print key,value


print "current angle, dx, dy, angle increment, cost mult, # style, line segment length, radius"

for angid in angids:

    primind = 0
    print 20*"=", " ", angid, " (", degrees[angid], ") ", 20*"="
    # Unit and Normal vector for the current direction
    uv = np.array(( np.cos(angles[angid]),np.sin(angles[angid])))
    nv = np.array((-np.sin(angles[angid]),np.cos(angles[angid])))
    #print "uv=",uv
    #print "nv=",nv

    # Handle the straight motions
    # Select grid points in the proper direction
    gd = np.dot(nv,grid)
    #print "gd=",gd
    gl = np.logical_and(np.logical_and(gd > -0.2, gd < 0.2), np.dot(uv,grid) > 0.2)
    #print "gd[logical]=",gd[gl]
    gs = grid[:, gl ]
    print "# Straight motion"
    #print "gs=",gs.shape, gs

    for ndx in range(gs.shape[1]):
        primind+= 1
        validate_grid(args, angid, primind, gs[0,ndx], gs[1,ndx], 0, resolution, uv, nv, 1)

    print "# Reverse motion"
    gl = np.logical_and(np.logical_and(gd > -0.2, gd < 0.2), np.dot(uv,grid)< - 0.2, )
    #print "gd[logical]=",gd[gl]
    gs = grid[:, gl ]
    gm = np.dot(uv,gs)
    ndx = np.argmax(gm)
    primind+= 1
    validate_grid(args, angid, primind, gs[0,ndx], gs[1,ndx], 0, resolution, uv, nv, 5)

    print "# ZRT motion"
    primind+= 1
    validate_grid(args, angid, primind,0, 0, 1, resolution, None, None, 5)
    primind+= 1
    validate_grid(args, angid, primind,0, 0,-1, resolution, None, None, 5)


    # Handle turns relative to current direction
    gd = np.dot(nv,grid)
    for delta in deltas:
        print "# ",delta,"/",n_angles," motion"
        for direction in (1,-1):
            #print angid, delta, direction

            # Select grid points in the proper direction
            gl = np.logical_and(direction*gd > 0.2, np.dot(uv,grid) > 0.2)
            gs = grid[:, gl]
            for ndx in range(gs.shape[1]):
                primind+= 1
                #print 30*"=",primind
                validate_grid(args, angid, primind, gs[0,ndx], gs[1,ndx], delta*direction, resolution, uv, nv, 1)

    for direction in (1,-1):
        gl = np.logical_and(direction*gd > 0.2, np.dot(uv,grid) > 0.2)
        gs = grid[:, gl]
        print "# Handle swerves to ",direction," where we don't change direction"
        for ndx in range(gs.shape[1]):
            primind+= 1
            #print 30*"=",primind
            validate_grid(args, angid, primind, gs[0,ndx], gs[1,ndx], 0, resolution, uv, nv, 4)
