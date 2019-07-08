#  /*
#   * Copyright (c) 2016-2019
#   *    Capable Humanitarian Robotics and Intelligent Systems Lab (CHRISLab)
#   *    Christopher Newport University
#   *    chrislab@cnu.edu
#   *
#   * All rights reserved.
#   *
#   * Based on genmprim_unicycle.m  Copyright (c) 2008, Maxim Likhachev
#   *
#   * Redistribution and use in source and binary forms, with or without
#   * modification, are permitted provided that the following conditions are met:
#   *
#   *     * Redistributions of source code must retain the above copyright
#   *       notice, this list of conditions and the following disclaimer.
#   *     * Redistributions in binary form must reproduce the above copyright
#   *       notice, this list of conditions and the following disclaimer in the
#   *       documentation and/or other materials provided with the distribution.
#   *     * Neither the name of the Carnegie Mellon University or
#   *       Christopher Newport University nor the names of its
#   *       contributors may be used to endorse or promote products derived from
#   *       this software without specific prior written permission.
#   *
#   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#   * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#   * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#   * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#   * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#   * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#   * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#   * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#   * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#   * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#   * POSSIBILITY OF SUCH DAMAGE.
#   */

import sys, math
import numpy as np

# Our custom utility methods
from mprim_generator_utils import *
from mprim_generator_generators import *


def get_motion_primitive_endpoint(angleind, primind, currentangle,
                                  primitive_generator_definitions,args):
    """
        Get the primitive generator and the angle of rotation
        based on the designated quadrant
        angleind     - index of angle in total angles (all quadrants)
        primind      - which generator to use for given angle
        currentangle - Starting angle
        primitive_generator_definitions - all the generator definitions
        args         - system options
    """
    resolution                 = args['resolution']
    numberofanglesperquadrant  = args['numberofanglesperquadrant']
    numberofprimsperangle      = args['numberofprimsperangle']

    quadrant       = int(angleind/numberofanglesperquadrant)
    generatorIndex = angleind%numberofanglesperquadrant

    quadrant_angle = quadrant*np.pi/2.0
    assert 0 <= quadrant <= 3, str("Invalid quadrant =  %d"%quadrant)

    #print "angle,quadrant, qa",angleind, quadrant, quadrant_angle
    return (quadrant_angle, primitive_generator_definitions[generatorIndex][primind][:])

    # # Return the rotation angle based on quadrant and the specific primitive generator
    # if (quadrant == 0):
    #     # First quadrant 0-90
    #     return (0., primitive_generator_definitions[generatorIndex][primind][:])
    # elif (quadrant == 1):
    #     # Second quadrant 90-180
    #     return (np.pi/2.0, primitive_generator_definitions[generatorIndex][primind][:])
    # elif (quadrant == 2):
    #     # Third quadrant 180-270
    #     return (np.pi, primitive_generator_definitions[generatorIndex][primind][:])
    # elif (quadrant == 3):
    #     # Fourth quadrant 270-360
    #     return (3.0*np.pi/2.0, primitive_generator_definitions[generatorIndex][primind][:])
    #
    # # Should not get here
    # raise Exception(str("Invalid quadrant =  %d"%quadrant))
    # return None


def generate_arc_line_arc_primitive(prim_def_primitive,args):
    """
        Connect start and end point using an arc, straightline, arc motion pattern
        CSC - Circle-straight line - Circle style
    """
    angle               = prim_def_primitive['angle']
    basemprimendpts_c   = prim_def_primitive['basemprim_end_pt']
    endpose_c           = prim_def_primitive['endpose']
    intermcells_m       = prim_def_primitive['intermcells']
    startpt             = prim_def_primitive['startpoint']
    endpt               = prim_def_primitive['endpoint']
    baseendpose_c       = basemprimendpts_c[0:3]
    #print "start pt=", startpt
    #print " End pt =", endpt

    # Internal values
    resolution              = args['resolution']
    numberofangles          = args['numberofanglesperquadrant']*4
    numberofprimsperangle   = args['numberofprimsperangle']
    numofsamples            = args['numberOfSamples']
    min_turning_radius_m    = args['minTurnRadius']
    wheelbase               = args['wheelbase']
    cost_conversion_factor   = args['cost_conversion_factor']


    # Standardize the geometry, then rotate solution
    if (startpt[0] != 0.0 or startpt[1] != 0.0):
        raise Exception(" Invalid starting point!")

    #print "Start pt=",startpt
    #print " End pt =",endpt
    local_startpt = np.array((0., 0., 0.))

    # Rotation matrix based on starting point heading
    rot = np.array(np.vstack((np.hstack(( np.cos(startpt[2]), np.sin(startpt[2]), 0.0)),
                              np.hstack((-np.sin(startpt[2]), np.cos(startpt[2]), 0.0)),
                              np.hstack((        0.0        ,          0.0      , 1.0)))))

    local_endpt = np.dot(rot, endpt)
    local_endpt[2] = endpt[2] - startpt[2]
    #print "Local start pt=",local_startpt
    #print "Local  End pt =",local_endpt

    v0 = np.array((1.0, 0.0))
    v1 = np.array((np.cos(local_endpt[2]), np.sin(local_endpt[2])))

    # Connecting vector between the start and end point positions
    vl = local_endpt[:2] - local_startpt[:2]
    shortest = np.linalg.norm(vl)
    #print "  vl=",vl," shortest=",shortest

    # Determine the directions orthogonal to the connecting vector
    (vc0,vc1,DT,DTO,DTx,DTy,DTquad) = get_normal_from_heading_vectors(v0,v1,vl)
    # print " v0=",v0," vc0=",vc0, " |vc0|=",np.linalg.norm(vc0)
    # print " v1=",v1," vc1=",vc1, " |vc1|=",np.linalg.norm(vc1)
    # print "  DT  =", DT, " DTO=",DTO," x,y=",DTx,", ",DTy," quad=",DTquad

    # Find the largest radius circle we can use to transition between points
    # Using two arcs  (solve r for |vl|^2 - (2r)^2 == 0)
    radius = min_turning_radius_m
    try:
        radius1 = get_max_radius_circles(local_startpt, vc0,local_endpt, vc1,min_turning_radius_m)
        radius2 = get_max_radius_circles(local_startpt, vc0,local_endpt,-vc1,min_turning_radius_m)
        radius3 = get_max_radius_circles(local_startpt,-vc0,local_endpt, vc1,min_turning_radius_m)
        radius4 = get_max_radius_circles(local_startpt,-vc0,local_endpt,-vc1,min_turning_radius_m)
        # Use the minimum radius, but should be > min_turning
        radius = np.min((radius1, radius2, radius3, radius4))
    except OSError as err:
        print "Failed to get maximum radius!"
        print("OS error: {0}".format(err))
    except ValueError as err:
        print "Failed to get maximum radius!"
        print("Value error: {0}".format(err))
    except TypeError as err:
        print "Failed to get maximum radius!"
        print("Type error: {0}".format(err))
    except KeyError as err:
        print "Failed to get maximum radius!"
        print("Key error: {0}".format(err))
    except NameError as err:
        print "Failed to get maximum radius!"
        print("Name error: {0}".format(err))
    except:
        print "Failed to get maximum radius!"
        print("Unexpected error:", sys.exc_info()[0])

    # Get a list of all possible tangents between all circles
    #    tangent to start/end points and headings
    potential_tangents = []

    potential_tangents.append(get_potential_internal_tangents(local_startpt,  vc0, local_endpt, v1,  vc1,radius))
    potential_tangents.append(get_potential_internal_tangents(local_startpt,  vc0, local_endpt, v1, -vc1,radius))
    potential_tangents.append(get_potential_internal_tangents(local_startpt, -vc0, local_endpt, v1,  vc1,radius))
    potential_tangents.append(get_potential_internal_tangents(local_startpt, -vc0, local_endpt, v1, -vc1,radius))
    potential_tangents.append(get_potential_external_tangents(local_startpt,  vc0, local_endpt, v1,  vc1,radius))
    potential_tangents.append(get_potential_external_tangents(local_startpt,  vc0, local_endpt, v1, -vc1,radius))
    potential_tangents.append(get_potential_external_tangents(local_startpt, -vc0, local_endpt, v1,  vc1,radius))
    potential_tangents.append(get_potential_external_tangents(local_startpt, -vc0, local_endpt, v1, -vc1,radius))
    sys.stdout.flush()

    for possibility in potential_tangents:
        if (possibility is not None):
            # Workable tangent (ptp,atp,vtn,pi,ptp1,atp1)
            #print possibility
            (valid,vc0,pt0,at0,rad0,vtn,pi,vc1,pt1,at1,rad1) = possibility
            if (valid):
                # print "Ready to rock! radius=",radius
                # print "  startpt=",startpt  , "  vc0=",vc0
                # print "     pt0 =",pt0      , "  at0=",at0, " rad0=",rad0
                # print "      pi =",pi       , "  vtn=",vtn
                # print "     pt1 =",pt1      , "  at1=",at1, " rad1=",rad1
                # print "   endpt =",endpt    , "  vc1=",vc1

                linear = np.linalg.norm(pt1-pt0)
                arcA = np.fabs(at0)*radius
                arcB = np.fabs(at1)*radius
                total_length = arcA + linear + arcB

                vel = total_length
                omega0 = vel/rad0
                omega1 = vel/rad1

                # print "vel=",vel," omega0=",omega0,"  at0=",at0," omega1=",omega1,"  at1=",at1," linear=",linear," at1=",at1
                # print "total length=",total_length

                # Define the intermediate points
                # Scaled time for linear section (Eqn (12) in www.sbpl.net/node/53)
                tl0  = arcA/vel
                tl1  = (arcA + linear)/vel
                # print "tl0=",tl0,"  tl1=",tl1

                # Calculate from standardized start, then rotate the points
                # generate samples assuming unit time interval
                for iind in range(0, numofsamples):
                    time = float(iind+1)/(numofsamples)
                    if (time < tl0):
                        # Eqn (9) in www.sbpl.net/node/53
                        dtheta = omega0*time
                        intermcells_m[iind,:] = [local_startpt[0] + rad0*np.sin(dtheta),
                                                 local_startpt[1] + rad0 - rad0*np.cos(dtheta),
                                                 dtheta ]
                        # print "arc A : ",time," ",time*vel," ",intermcells_m[iind,:]," dTheta=",dtheta," of ",at0
                    elif (time < tl1):
                        # Eqn (2) in www.sbpl.net/node/53
                        intermcells_m[iind,:] = [pt0[0]+(time-tl0)*vel*np.cos(at0),
                                                 pt0[1]+(time-tl0)*vel*np.sin(at0),
                                                 at0 ]
                        # print "line1: ",time," ",time*vel," ",intermcells_m[iind,:]

                    else:
                        # Eqn (9) in www.sbpl.net/node/53
                        omega1 = np.sign(at1)*vel/radius
                        dtheta = omega1*(time-tl1)
                        intermcells_m[iind,:] = [pt1[0] - rad1*np.sin(at0) + rad1*np.sin(dtheta)*np.cos(at0) + rad1*np.cos(dtheta)*np.sin(at0),
                                                 pt1[1] + rad1*np.cos(at0) + rad1*np.sin(dtheta)*np.sin(at0) - rad1*np.cos(dtheta)*np.cos(at0),
                                                 dtheta + at0 ]
                        # print "arc B : ",time," ",time*vel," ",intermcells_m[iind,:]," dTheta=",dtheta," of ",at1

                # Rotate relative to actual starting angle
                # Rotation matrix based on starting point heading
                #print "---- Before rotating Arc-Line-Arc positioning ----"
                #print "intermcells_m=",intermcells_m
                #print "---- Rotate to final Arc-Line-Arc positioning ----"
                for iind in range(0, numofsamples):
                    intermcells_m[iind,:3] = np.dot(intermcells_m[iind,:3],rot)
                    intermcells_m[iind,2] += startpt[2]

                # print "intermcells_m=",intermcells_m

                # correct
                errorxy = np.array((endpt[0] -intermcells_m[numofsamples-1,0],
                                    endpt[1] -intermcells_m[numofsamples-1,1]))

                error_theta = math.fabs(shortest_angle_dist(endpt[2], intermcells_m[int(numofsamples)-1,2]))

                if np.linalg.norm(errorxy) > 1e-5 or error_theta > 1.e-3:
                    # We have an error, so correct it, but this should not be frequent
                    print "---- Correct Arc-Line-Arc positioning ----"
                    print "post rotated intermcells_m=",intermcells_m
                    print "endpt=",endpt
                    print('errx=%f erry=%f\n'%(errorxy[0], errorxy[1]))
                    interpfactor = \
                         np.array((np.arange(0.,
                                             1.+(1./(numofsamples)),
                                             1./(numofsamples-1))))

                    print "interp'=",interpfactor.conj().T

                    intermcells_m[:,0] = intermcells_m[:,0]+errorxy[0]*interpfactor.conj().T
                    intermcells_m[:,1] = intermcells_m[:,1]+errorxy[1]*interpfactor.conj().T

                    print " corrected:"
                    print " startpt=",startpt
                    print "intermcells_m=",intermcells_m
                    print "  endpt =",endpt

                errorxy = np.array((endpt[0] -intermcells_m[numofsamples-1,0],
                                    endpt[1] -intermcells_m[numofsamples-1,1]))

                error_theta = math.fabs(shortest_angle_dist(endpt[2], intermcells_m[int(numofsamples)-1,2]))

                if np.linalg.norm(errorxy) > 1e-5 or error_theta > 1.e-3:
                    print "---- Error in Arc-Line-Arc positioning ----"
                    print "post rotated intermcells_m=",intermcells_m
                    print "endpt=",endpt
                    print('errx=%f erry=%f errtheta=%f\n'%(errorxy[0], errorxy[1],error_theta))
                    return False

                #plt.figure()
                #plt.plot(intermcells_m[:,0],intermcells_m[:,1],marker='.')
                #plt.axis('equal');
                #print "----------------"

                #print "Wait for button press for next primitive generation"
                #plt.waitforbuttonpress()  # uncomment to plot each primitive set one at a time

                # Yeah! We found a solution using arc-line-arc
                prim_def_primitive['style']            = "Arc-Line-Arc"
                prim_def_primitive['turning_radius']   = radius
                prim_def_primitive['line length']      = vel*(tl1-tl0)

                # Radius to outside wheel * rotation angle is max arc length + linear distance
                # Add resolution to encourage straight motions
                travel_dist = (math.fabs(at0) + math.fabs(at1))*(math.fabs(radius) + 0.5*wheelbase) + linear + resolution
                prim_def_primitive['actioncost'] = int(cost_conversion_factor*travel_dist + 0.5)
                # print "   ac line arc=", prim_def_primitive['actioncost'] ,travel_dist, args['wheelbase'], radius, cost_conversion_factor
                return True

            else:
                pass
                # print "Invalid vc0=",vc0," vc1=",vc1
                #print possibility


    #print "No valid tangents for arc-line-arc ..."
    if False:
        # Visualize the definitions during debug
        get_potential_internal_tangents(local_startpt,  vc0, local_endpt, v1,  vc1,radius,True)
        get_potential_internal_tangents(local_startpt,  vc0, local_endpt, v1, -vc1,radius,True)
        get_potential_internal_tangents(local_startpt, -vc0, local_endpt, v1,  vc1,radius,True)
        get_potential_internal_tangents(local_startpt, -vc0, local_endpt, v1, -vc1,radius,True)
        get_potential_external_tangents(local_startpt,  vc0, local_endpt, v1,  vc1,radius,True)
        get_potential_external_tangents(local_startpt,  vc0, local_endpt, v1, -vc1,radius,True)
        get_potential_external_tangents(local_startpt, -vc0, local_endpt, v1,  vc1,radius,True)
        get_potential_external_tangents(local_startpt, -vc0, local_endpt, v1, -vc1,radius,True)
    #print "                  - return failure!"

    return False

def generate_line_arc_primitive(prim_def_primitive,args):
    """
        Original approach used in SBPL (http://sbpl.net/node/53)
        Move straight from current heading, then arc to final heading
    """
    # unicycle-based move forward or backward  (http://sbpl.net/node/53)
    angle               = prim_def_primitive['angle']
    basemprimendpts_c   = prim_def_primitive['basemprim_end_pt']
    endpose_c           = prim_def_primitive['endpose']
    intermcells_m       = prim_def_primitive['intermcells']
    startpt             = prim_def_primitive['startpoint']
    endpt               = prim_def_primitive['endpoint']
    baseendpose_c       = basemprimendpts_c[0:3]

    # Internal values
    resolution               = args['resolution']
    numberofangles           = args['numberofanglesperquadrant']*4
    numberofprimsperangle    = args['numberofprimsperangle']
    numofsamples             = args['numberOfSamples']
    minTurnRadius            = args['minTurnRadius']
    cost_conversion_factor   = args['cost_conversion_factor']
    wheelbase                = args['wheelbase']

    # unicycle-based move forward or backward  (http://sbpl.net/node/53)
    R = np.array(np.vstack((np.hstack((np.cos(startpt[2]),
                                       np.sin(endpt[2])-np.sin(startpt[2]))),
                            np.hstack((np.sin(startpt[2]),
                                      -np.cos(endpt[2])+np.cos(startpt[2]))) )))

    S = np.dot(np.linalg.pinv(R), np.array(np.vstack((endpt[0]-startpt[0],
                                                      endpt[1]-startpt[1])) ))
    l = S[0][0]       # l in www.sbpl.net/node/53  eqn (10)
    if l < -1.e-3 or math.fabs(S[1][0]) < minTurnRadius:
        # Requires a starting point before current point
        # if a very small number, we will just presume that we can make it  and clamp to zero

        # print ' The line-arc formulation from www.sbpl.net/node/53 is NOT allowed l=%.6f'%(l)
        return False

    else:

        if l < 0.0:
            #print "Clamping line to zero for pure turn! l=%.6f"%(l)
            l = 0.0

        # Can do a line-arc (as original www.sbpl.net/node/53) or line-arc-line
        # radius is the maximum radius that we can use with line-arc-0
        #print('The line-arc formulation from www.sbpl.net/node/53 is allowed l=%f'%(l))
        radius = S[1][0]  # r in eqn (10)

        rotation_angle = wrap_angle(endpt[2]-startpt[2])
        omega = rotation_angle+ l/radius # w in eqn (12) www.sbpl.net/node/53 (unit time scaling)
        vel   = radius*omega # v in eqn (12) www.sbpl.net/node/53  (in terms of unit time scaling)


        # Scaled time for linear section (Eqn (12) in www.sbpl.net/node/53)
        tl  = l/vel

        #print "omega=",omega
        #print "vel=",vel
        #print "scaled tl=",tl

        # compute rv
        # rv = baseendpose_c(3)*2*pi/numberofangles;
        # compute tv
        # tvx = (endpt(1) - startpt(1))*rv/(sin(endpt(3)) - sin(startpt(3)))
        # tvy = -(endpt(2) - startpt(2))*rv/(cos(endpt(3)) - cos(startpt(3)))
        # tv = (tvx + tvy)/2.0;
        if math.isnan(radius) or math.isnan(vel) or math.fabs(radius) < 1e-8:
            print "Invalid radius or velocity in line-arc generator", radius, vel
            print "R=\n",R
            print "Rpi=\n",np.linalg.pinv(R)
            print "S=\n",S
            print "l=",l
            print "radius=",radius
            return False

        # generate samples
        # print "angle=",angle, " startpt= ",startpt

        for iind in range(0, numofsamples):
            dt = float(iind+1)/(numofsamples)
            # dtheta = rv*dt + startpt(3);
            # intermcells_m(iind,:) = [startpt(1) + tv/rv*(sin(dtheta) - sin(startpt(3))) ...
            #                         startpt(2) - tv/rv*(cos(dtheta) - cos(startpt(3))) ...
            #                         dtheta];
            if math.isnan(dt) or dt < 0.0:
                print "Invalid time step in line-arc generator", dt, iind
                return False

            if (dt*vel)<l:
                # Eqn (2) in www.sbpl.net/node/53
                intermcells_m[iind,:] = [startpt[0]+dt*vel*np.cos(startpt[2]),
                                         startpt[1]+dt*vel*np.sin(startpt[2]),
                                         startpt[2] ]

                first_arc = 1.*intermcells_m[iind,:]
            else:
                # Eqn (9) in www.sbpl.net/node/53
                dtheta = omega*(dt-tl) + startpt[2]
                intermcells_m[iind,:] = [startpt[0] + l*np.cos(startpt[2]) + radius*(np.sin(dtheta)-np.sin(startpt[2])),
                                         startpt[1] + l*np.sin(startpt[2]) - radius*(np.cos(dtheta)-np.cos(startpt[2])),
                                         dtheta ]

        # correct
        errorxy = np.array((endpt[0] -intermcells_m[int(numofsamples)-1,0],
                            endpt[1] -intermcells_m[int(numofsamples)-1,1]))

        if np.linalg.norm(errorxy) > 1e-3:
            print('l=%f errx=%f erry=%f\n'%(l, errorxy[0], errorxy[1]))
            interpfactor = np.array((np.arange(0., 1.+(1./(numofsamples)), 1./(numofsamples-1))))

            #print "intermcells_m=",intermcells_m
            print "interp'=",interpfactor.conj().T

            intermcells_m[:,0] = intermcells_m[:,0]+errorxy[0]*interpfactor.conj().T
            intermcells_m[:,1] = intermcells_m[:,1]+errorxy[1]*interpfactor.conj().T

        errorxy = np.array((endpt[0] -intermcells_m[int(numofsamples)-1,0],
                            endpt[1] -intermcells_m[int(numofsamples)-1,1]))

        error_theta = math.fabs(shortest_angle_dist(endpt[2], intermcells_m[int(numofsamples)-1,2]))

        prim_def_primitive['turning_radius']   = radius
        prim_def_primitive['style']            = "Line-Arc"
        prim_def_primitive['line length']      = l

        # Radius to outside wheel * rotation angle is max arc length + linear distance
        # Add resolution to encourage straight motions
        travel_dist = math.fabs(rotation_angle)*(math.fabs(radius) + 0.5*wheelbase) + l + resolution
        prim_def_primitive['actioncost'] = int(cost_conversion_factor*travel_dist + 0.5)
        #print "   ac arc=", prim_def_primitive['actioncost'] ,rotation_angle, args['wheelbase'], radius, cost_conversion_factor

        if np.linalg.norm(errorxy) > 1e-3 or error_theta > 1e-3:
            print "Invalid error for line-arc!",errorxy, error_theta
            return False
        #else:
        #    print "Final line-arc error=",errorxy, error_theta

        return True


    raise Exception(" Invalid return ")
    return False


def generate_arc_based_primitives(prim_def_primitive,args):
    """
        Determine which type of generator we intend to use
    """
    if (args['check_arc_line_arc_generator']):
        # Prefer the original genmprim approach
        if not generate_line_arc_primitive(prim_def_primitive,args):
            #print "line-arc mode failed - try arc-line-arc shortest path"
            return generate_arc_line_arc_primitive(prim_def_primitive,args)
        else:
            #print "use valid line-arc construct instead of shortest arc-line-arc path"
            return True
    else:
        # Use the original approach of line followed by single arc
        return generate_line_arc_primitive(prim_def_primitive,args)

    raise Exception(" Invalid return")
    return

def generate_first_quadrant_primitive(angleind,     primind,
                                      currentangle, basemprimendpts_c,
                                      args, prim_def_primitive ):
    """
    Define primitive generator for first quadrant angle
    prim_def_primitive dictionary to hold output
    """

    resolution               = args['resolution']
    numberofangles           = args['numberofanglesperquadrant']*4
    numberofprimsperangle    = args['numberofprimsperangle']
    numofsamples             = args['numberOfSamples']
    minTurnRadius            = args['minTurnRadius']
    angles                   = args['angles']
    wheelbase                = args['wheelbase']
    cost_conversion_factor   = args['cost_conversion_factor']

    # Define action to terminal end point
    baseendpose_c = basemprimendpts_c[0:3]
    #print "bep=",baseendpose_c

    endx_c        = np.round(baseendpose_c[0]) # First quadrant angle is 0.0, not need to rotate
    endy_c        = np.round(baseendpose_c[1])
    endtheta_c    = (angleind+baseendpose_c[2])%numberofangles
    endpose_c     = np.array((endx_c, endy_c, endtheta_c))
    #print "epc=",endpose_c

    startpt = np.array((0., 0., currentangle))
    endpt   = np.array(((endpose_c[0]*resolution),
                        (endpose_c[1]*resolution),
                         angles[endtheta_c] ) )

    # print "\n------ Base Generator Definition -----------"
    # print "     angleind=",angleind,"  primind=",primind
    # print " basemprimendpts_c=",basemprimendpts_c
    # print " baseendpose_c =",baseendpose_c
    # print "     endpose_c =",endpose_c
    # print "      startpt  =",startpt
    # print "      endpt    =",endpt

    # Generate intermediate poses (remember they are w.r.t 0,0 (and not centers of the cells)
    intermcells_m = np.zeros((numofsamples, 3))

    # Store the specific motion primitive data in dictionary
    # Fill in the primitive definition data
    prim_def_primitive['angle']            = 0.0
    prim_def_primitive['basemprim_end_pt'] = basemprimendpts_c
    prim_def_primitive['endpose']          = endpose_c
    prim_def_primitive['intermcells']      = intermcells_m
    prim_def_primitive['startpoint']       = startpt
    prim_def_primitive['endpoint']         = endpt
    prim_def_primitive['style']            = "None"
    prim_def_primitive['turning_radius']   = np.Inf
    prim_def_primitive['line length']      = 0.0
    # Retrive the cost multiplier (-99 for invalid primitive)
    prim_def_primitive['actioncostmult']   = basemprimendpts_c[3]
    prim_def_primitive['actioncost']       = 10**9

    line_angle = np.math.atan2(endpt[1]-startpt[1], endpt[0]-startpt[0])
    line_angle = wrap_angle(line_angle-startpt[2])
    check_line_angle = np.math.fabs(np.math.sin(line_angle)) < 0.1

    # print " line_angle=",line_angle," check=",check_line_angle, " end x=",endx_c," y=",endy_c," endpose change=",baseendpose_c[2]
    # print "    zrt=",np.logical_and(endx_c == 0., endy_c == 0.),"  same heading=",np.logical_and(check_line_angle,baseendpose_c[2] == 0.)
    if np.logical_and(endx_c == 0., endy_c == 0.):             # turn in place or move forward along current heading has simple interpolation
        prim_def_primitive['style']            = "Zero-radius turn"
        prim_def_primitive['turning_radius']   = 0.0
        rotation_angle = wrap_angle(endpt[2] - startpt[2])
        for iind in range(0, numofsamples):
            fraction = float(iind+1)/(numofsamples)
            intermcells_m[iind,:] = np.array((startpt[0]+(endpt[0]-startpt[0])*fraction,
                                              startpt[1]+(endpt[1]-startpt[1])*fraction, 0))

            intermcells_m[iind,2] = (startpt[2]+rotation_angle*fraction)%(2.*np.pi)
            #print " ",iind,"  of ",numofsamples," fraction=",fraction," rotation=",rotation_angle

        # Define action cost as the total distance each wheel moves during motion
        # Deliberately double counting wheel motion for zero-radius-turns to discourage slightly
        # Add the minimum resolution to the arc length to prefer straighter motions with turn than full stop
        prim_def_primitive['actioncost']       = int((math.fabs(rotation_angle)*wheelbase + resolution)*cost_conversion_factor + 0.5)
        # print "   ac zrt=", prim_def_primitive['actioncost'] ,rotation_angle, wheelbase, cost_conversion_factor

    elif np.logical_and(check_line_angle,baseendpose_c[2] == 0.):
        prim_def_primitive['style']            = "Straight"
        prim_def_primitive['line length']      = np.linalg.norm(np.array([startpt[0]-endpt[0], startpt[1]-endpt[1]]))


        for iind in range(0, numofsamples):
            fraction = float(iind+1)/(numofsamples)
            intermcells_m[iind,:] = np.array((startpt[0]+(endpt[0]-startpt[0])*fraction,
                                              startpt[1]+(endpt[1]-startpt[1])*fraction, 0))
            rotation_angle = wrap_angle(endpt[2] - startpt[2])

            intermcells_m[iind,2] = (startpt[2]+rotation_angle*fraction)%(2.*np.pi)
            # print " ",iind,"  of ",numofsamples," fraction=",fraction," rotation=",rotation_angle

        # Define action cost as the maximum distance a wheel moves during motion
        prim_def_primitive['actioncost']       = int(prim_def_primitive['line length']*cost_conversion_factor + 0.5)
        # print "   ac lin=", prim_def_primitive['actioncost'] , prim_def_primitive['line length'], wheelbase, cost_conversion_factor
    else:
        # Use several geometric techniques to determine the motion primitives
        if (not generate_arc_based_primitives(prim_def_primitive, args)):
            basemprimendpts_c[3] = -99 # Flag invalid data (negative cost multiplier)
            #print "Failed to generate a valid motion primitive"
        #else:
        #    print "Successful generation"

    #print "   action cost=", prim_def_primitive['actioncost']

    return basemprimendpts_c[3] > 0


def get_primitive_definition(angleind, primind, currentangle, primitive_definitions, primitive_generator_definitions,args):
    """
        Define the motion primitives for this grid point
    """
    # Get parameters used in generation
    resolution                = args['resolution']
    numberofanglesperquadrant = args['numberofanglesperquadrant']
    numberofprimsperangle     = args['numberofprimsperangle']
    numofsamples              = args['numberOfSamples']
    numberofangles            = numberofanglesperquadrant*4

    # Get the references and add as necessary
    while (angleind >= len(primitive_definitions)):
        # print " Adding primitive definition for ",len(primitive_definitions), " of ",angleind," angles"
        primitive_definitions.append([])

    if (len(primitive_definitions) != angleind+1):
        print "len(primitive_definitions) = ",len(primitive_definitions)
        print "   primind = ",angleind
        sys.stdout().flush()
        raise Exception(" Invalid handling of primitive_definitions!")

    prim_def_angle_list = primitive_definitions[angleind] # reference to angle list

    while (primind >= len(prim_def_angle_list)):
        # print " Adding primitive definition for ",len(prim_def_angle_list), " of ",primind," for angle ",angleind
        prim_def_angle_list.append(dict())

    # Store the specific motion primitive data in dictionary
    if (len(prim_def_angle_list) != primind+1):
        print "len(prim_def_angle_list) = ",len(prim_def_angle_list)
        print "   primind = ",primind
        sys.stdout().flush()
        raise Exception(" Invalid handling of prim_def_angle_list!")

    # Reference to the dictionary for primitives we are definining here
    prim_def_primitive  = prim_def_angle_list[primind]

    # Get standardized end point data
    quadrant_angle, basemprimendpts_c = get_motion_primitive_endpoint(angleind,
                                                primind, currentangle,
                                                primitive_generator_definitions,args)

    if (quadrant_angle == 0.):
        #print "Define base primitive from generator ..."
        if not generate_first_quadrant_primitive(angleind, primind,
                                                 currentangle, basemprimendpts_c,
                                                 args,
                                                 prim_def_primitive ):
            #print "Failed to define valid motion primitive for first quadrant", basemprimendpts_c
            pass


    else:
        #print " Rotate existing primitive by "+ str(angle)+" to get next quadrant (" + str(primind)+" " + str(angleind) + ") ..."

        # Get corresponding primitive in the first quadrant
        corresponding_angleind = angleind%numberofanglesperquadrant
        corresponding_prim_def_angle_list = primitive_definitions[corresponding_angleind] # reference to angle list
        corresponding_angleind_primitive_def = corresponding_prim_def_angle_list[primind]

        # Define action to terminal end point for this particular angle
        baseendpose_c = basemprimendpts_c[0:3]
        endx_c        = np.round((baseendpose_c[0]*np.cos(quadrant_angle))-(baseendpose_c[1]*np.sin(quadrant_angle)))
        endy_c        = np.round((baseendpose_c[0]*np.sin(quadrant_angle))+(baseendpose_c[1]*np.cos(quadrant_angle)))
        endtheta_c    = (angleind+baseendpose_c[2])%numberofangles
        endpose_c     = np.array((endx_c, endy_c, endtheta_c))


        startpt = np.array((0., 0., currentangle))
        endpt   = np.array(((endpose_c[0]*resolution),
                            (endpose_c[1]*resolution),
                             args['angles'][endtheta_c] ) )


        #print "\n----- Rotated Primitive from Base Generator -----------"
        #print "    angleind=",angleind,"  primind=",primind
        #print "    endpose_c=",endpose_c
        #print "    endtheta_c=",endtheta_c, " baseendpose_c =",baseendpose_c
        #print("    quadrant rotation angle=%f\n"% (angle*180./np.pi))
        #print "    startpt =",startpt
        #print "    endpt   =",endpt

        # Generate intermediate poses (remember they are w.r.t 0,0 (and not centers of the cells)
        intermcells_m = np.zeros((numofsamples, 3))

        # Fill in the primitive definition data
        prim_def_primitive['angle']            = quadrant_angle
        prim_def_primitive['basemprim_end_pt'] = basemprimendpts_c
        prim_def_primitive['endpose']          = endpose_c
        prim_def_primitive['intermcells']      = intermcells_m
        prim_def_primitive['startpoint']       = startpt
        prim_def_primitive['endpoint']         = endpt
        prim_def_primitive['actioncost']       = corresponding_angleind_primitive_def['actioncost']
        prim_def_primitive['actioncostmult']   = corresponding_angleind_primitive_def['actioncostmult']
        prim_def_primitive['style']            = corresponding_angleind_primitive_def['style']
        prim_def_primitive['turning_radius']   = corresponding_angleind_primitive_def['turning_radius']

        # Rotate the corresponding primitive into this quadrant
        # Use transpose and post multiply row vector of positions
        rot = np.array(np.vstack((np.hstack(( np.cos(quadrant_angle), np.sin(quadrant_angle), 0.0)),
                                  np.hstack((-np.sin(quadrant_angle), np.cos(quadrant_angle), 0.0)),
                                  np.hstack((          0            ,             0.        , 1.0)))))
        for iind in range(0, numofsamples):
            intermcells_m[iind,:] = np.dot(corresponding_angleind_primitive_def['intermcells'][iind,:],rot)
            intermcells_m[iind,2] += quadrant_angle


        # End of rotating exiting primitive block

def generate_motion_primitives(primitive_generator_definitions,args):
    """
        Generate all motion primitives for this setup
    """
    print "Generate motion primitive definitions ..."
    primitive_definitions = [] # List to hold data for each angle index

    resolution                = args['resolution']
    numberofanglesperquadrant = args['numberofanglesperquadrant']
    numberofprimsperangle     = args['numberofprimsperangle']
    numofsamples              = args['numberOfSamples']
    numberofangles            = numberofanglesperquadrant*4

     # iterate over angles
    for angleind in range(0, numberofangles):
        primitive_definitions.append([]) # list to hold data to hold data for each primitive within angle

        prim_def_angle_list = primitive_definitions[angleind] # reference to angle list

        currentangle = args['angles'][angleind]

        # iterate over primitives
        for primind in range(0, numberofprimsperangle):
            # Add the newest generator
            prim_def_angle_list.append(dict())                # Add this primitive to angle list

            # Get the specific primitive definition and add
            #     to primitive_definitions list
            get_primitive_definition(angleind, primind, currentangle, primitive_definitions, primitive_generator_definitions,args)

    print "Finished primitive definitions!"
    return primitive_definitions
