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
#from genmprim_extended_utils import *

def get_default_primitive_generators(args):
    """
        Get the primitive generators for the first quadrant
            (assume symmetric for other for quadrants)
        based on the default implementation of SBPL
    """
    base_mprim_end_points = dict()  # Create a dictionary to hold the generators


    if (args['numberofprimsperangle'] != 10):
        print "default primitives use 10 primitives per angle!"
        args['numberofprimsperangle'] = 10

    if (args['numberofanglesperquadrant'] != 4):
        print "default primitives use 4 angles per quadrant!"
        args['numberofanglesperquadrant'] = 4


    numberofanglesperquadrant = args['numberofanglesperquadrant']
    numberofangles = numberofanglesperquadrant*4

    numberofprimsperangle     = args['numberofprimsperangle']

    if (numberofanglesperquadrant != 4):
        raise Exception(" This file assumes 4 angles per quadrant!")

    # Cost multipliers (multiplier is used as costmult*cost)
    forwardcostmult         = args['forwardcostmult']
    forwardandturncostmult  = args['forwardandturncostmult']
    backwardcostmult        = args['backwardcostmult']
    turninplacecostmult     = args['turninplacecostmult']
    sidestepcostmult        = args['sidestepcostmult']

    # Define the list of angles for the first quadrant (assuming 4 angles)
    args['angles'] = np.zeros(numberofangles)
    args['angles'][0] = 0.0                     #  0.0 degrees
    args['angles'][1] = math.atan2(2.,2.)*0.5  # 22.5 degrees
    args['angles'][2] = math.atan2(2.,2.)      # 45.0 degrees
    args['angles'][3] = math.atan2(2.,2.)*1.5  # 67.5 deg = 90-22.5
    #args['angles'][1] = math.atan2(1.,3.)       # 18.4 degrees
    #args['angles'][2] = math.atan2(2.,2.)       # 45.0 degrees
    #args['angles'][3] = math.atan2(3.,1.)       # 71.6 degrees = 90 - 18.4

    # Initialize the other 3 quadrants
    for ndx in range(1,4):
        args['angles'][ndx*4 + 0] = float(ndx)*np.pi/2.0 + args['angles'][0]
        args['angles'][ndx*4 + 1] = float(ndx)*np.pi/2.0 + args['angles'][1]
        args['angles'][ndx*4 + 2] = float(ndx)*np.pi/2.0 + args['angles'][2]
        args['angles'][ndx*4 + 3] = float(ndx)*np.pi/2.0 + args['angles'][3]


    # Initialize list to hold primitive data for each angle
    for ang in range(0,numberofanglesperquadrant):
        base_mprim_end_points[ang] = [] # x,y,delta theta,costmult

    # Define the end points for each angle in the first quadrant
    # x aligned with the heading of the robot, angles are positive counter clockwise
    # note, what is shown x,y,theta changes (not absolute numbers)

    # Aligned 0 degrees
    # 0 theta change
    base_mprim_end_points[0].append(np.array(( 1.,  0.,  0., forwardcostmult)))
    base_mprim_end_points[0].append(np.array((10.,  0.,  0., forwardcostmult)))
    base_mprim_end_points[0].append(np.array(( 5.,  0.,  0., forwardcostmult)))
    base_mprim_end_points[0].append(np.array((-1.,  0.,  0., backwardcostmult)))
    # 1/16 theta change
    base_mprim_end_points[0].append(np.array(( 6.,  1.,  1., forwardandturncostmult)))
    base_mprim_end_points[0].append(np.array(( 6., -1., -1., forwardandturncostmult)))
    base_mprim_end_points[0].append(np.array(( 9.,  1.,  1., forwardandturncostmult)))
    base_mprim_end_points[0].append(np.array(( 9., -1., -1., forwardandturncostmult)))
    # turn in place
    base_mprim_end_points[0].append(np.array(( 0.,  0.,  1., turninplacecostmult)))
    base_mprim_end_points[0].append(np.array(( 0.,  0., -1., turninplacecostmult)))

    #Aligned to 22.5 degrees on grid
    # 0 theta change
    base_mprim_end_points[1].append(np.array(( 2.,  1.,  0., forwardcostmult)))
    base_mprim_end_points[1].append(np.array(( 8.,  4.,  0., forwardcostmult)))
    base_mprim_end_points[1].append(np.array(( 4.,  2.,  0., forwardcostmult)))
    base_mprim_end_points[1].append(np.array((-2., -1.,  0., backwardcostmult)))
    # 1/16 theta change
    base_mprim_end_points[1].append(np.array(( 5.,  4.,  1., forwardandturncostmult)))
    base_mprim_end_points[1].append(np.array(( 6.,  1., -1., forwardandturncostmult)))
    base_mprim_end_points[1].append(np.array(( 8.,  5.,  1., forwardandturncostmult)))
    base_mprim_end_points[1].append(np.array(( 8.,  3., -1., forwardandturncostmult)))
    # turn in place
    base_mprim_end_points[1].append(np.array((0., 0.,  1., turninplacecostmult)))
    base_mprim_end_points[1].append(np.array((0., 0., -1., turninplacecostmult)))

    # Aligned to 45 degrees on grid
    # 0 theta change
    base_mprim_end_points[2].append(np.array(( 1.,  1.,  0., forwardcostmult)))
    base_mprim_end_points[2].append(np.array(( 7.,  7.,  0., forwardcostmult)))
    base_mprim_end_points[2].append(np.array(( 4.,  4.,  0., forwardcostmult)))
    base_mprim_end_points[2].append(np.array((-1., -1.,  0., backwardcostmult)))
    # 1/16 theta change
    base_mprim_end_points[2].append(np.array(( 4.,  6.,  1., forwardandturncostmult)))
    base_mprim_end_points[2].append(np.array(( 6.,  4., -1., forwardandturncostmult)))
    base_mprim_end_points[2].append(np.array(( 6.,  8.,  1., forwardandturncostmult)))
    base_mprim_end_points[2].append(np.array(( 8.,  6., -1., forwardandturncostmult)))
    # turn in place
    base_mprim_end_points[2].append(np.array(( 0.,  0.,  1., turninplacecostmult)))
    base_mprim_end_points[2].append(np.array(( 0.,  0., -1., turninplacecostmult)))

    # Aligned to 67.5 degrees on grid  (map the 22.5 values)
    for primind in range(0,numberofprimsperangle):
        base_mprim_end_points[3].append(np.array(( base_mprim_end_points[1][primind][1],
                                                   base_mprim_end_points[1][primind][0],
                                               -1.*base_mprim_end_points[1][primind][2],
                                                   base_mprim_end_points[1][primind][3] )))

    return base_mprim_end_points

def read_primitive_generator_definitions(args):
    """
        Read in the primitive generator definitions from a CSV text file
    """
    import csv
    base_mprim_end_points = dict()  # define dictionary to hold generator data

    with open(args['input'], 'r') as csvfile:
        primreader = csv.reader(csvfile)
        initialized = False
        getAngles   = False
        getPrims    = False
        angle_cnt   = 0
        for i, row in enumerate(primreader):
            #print "row:",row
            if (row == None or len(row) < 1 or row[0][0] == '#'):
                    # This line is a comment, so ignore
                    #print "Skip: ",row
                    continue

            if (not initialized):
                if (len(row) == 2 or len(row) == 3):
                    if (len(row) == 3 and row[2][0] != '#'):
                        print "Invalid data before initialization!"
                        print row
                        raise Exception(str("Invalid data in %s"%args['input']))

                    print "Initialize primitive setup ..."
                    initialized = True

                    args['numberofanglesperquadrant'] = int(row[0])
                    args['numberofprimsperangle'] = int(row[1])


                    numberofangles = args['numberofanglesperquadrant']*4
                    args['angles'] = np.zeros(numberofangles) # initialize

                    if (args['numberofanglesperquadrant']< 1):
                        print "Invalid number of angles=",args['numberofangles']," - must be divisible by 4 quadrants!"
                        raise Exception(str("Invalid number of angles=%d - must be divisible by 4 quadrants!"%args['numberofangles']))

                    if (args['numberofprimsperangle'] < 1):
                        print "Invalid number of primitives=",args['numberofprimsperangle']
                        raise Exception(str("Invalid number of primitives=%d !"%args['numberofprimsperangle']))

                    # Only need to initialize first quadrant angles
                    for ang in range(0,args['numberofanglesperquadrant']):
                        base_mprim_end_points[ang] = [] # initialize elements of dictionary into a list

                else:
                    print "Invalid data before initialization!"
                    print row
                    raise Exception(str("Invalid data in %s"%args['input']))
            else:
                # start angle, end_x, end_y, end_angle, costmult
                #print "getAngles=",getAngles,"  getPrims=",getPrims

                if (getAngles):
                    # Process angle data until we see primitive flag
                    if (len(row) == 1):
                        if (row[0] != "Primitives"):
                            # Expect primitives next!
                            print "Invalid field!"
                            raise Exception(str("Invalid field during get angles : %s"%row))
                        else:
                            # Switch to getting primitives
                            print "Flag ready to load primitives"
                            getPrims = True
                            getAngles = False

                            # Validate the angles
                            if (angle_cnt != args['numberofanglesperquadrant']):
                                print "Invalid count of initial angles !"
                                raise Exception(str("Invalid count of initial angles ! : %s"%row))

                            # Validate the angles
                            start = -0.0001
                            valid = True
                            for ang in args['angles']:
                                if ang < start:
                                    valid = False
                                    print "ERROR: ",ang," is decreasing"
                            if not valid:
                                raise Exception(str("Invalid angle - must be monotonic increasing : %s"%args['angles']))

                            # Initialize the remaining angle list
                            for ndx1 in range(1,4):
                                for ndx2 in range(0,args['numberofanglesperquadrant']):
                                    args['angles'][ndx1*4 + ndx2] = float(ndx1)*np.pi/2.0 + args['angles'][ndx2]

                            print "Angle List: radians (degrees)"
                            for ang in args['angles']:
                                print ang, " (",ang*180/np.pi,")"


                    else:
                      # Getting the angle data now
                      try:
                        if len(row) == 1 or \
                            (len(row) == 2 and row[1].strip()[0] == '#'):
                            args['angles'][angle_cnt] = float(row[0]) # Just radians

                        elif len(row) == 2 or \
                            (len(row) == 3 and row[2].strip()[0] == '#'):
                            # get dx,dy
                            dx = float(row[0])
                            dy = float(row[1])

                            args['angles'][angle_cnt] = math.atan2(dy,dx)
                            #print "  angle ",angle_cnt,"  = ",args['angles'][angle_cnt]
                            angle_cnt = angle_cnt + 1

                        else:
                            print "Invalid field!"
                            raise Exception(str("Invalid field during get angles : %s"%row))

                        if (angle_cnt > args['numberofanglesperquadrant']):
                            print " Only set the first quadrant angles!"
                            raise Exception(str("Invalid Only set the first quadrant angles! : %s"%row))

                      except Exception as e:
                          print e
                          raise Exception(str("Invalid field during get angles : %s"%row))

                elif (getPrims):
                    # Now get the primitive defintions
                    if (len(row) < 5 or len(row) > 6):
                        print "Invalid data in row ",i," <",row,">"
                        raise Exception(str("Invalid data in %s"%args['input']))

                    if (len(row) ==6):
                        comment = row[5].strip()
                        if (len(comment)  and comment[0] != '#'):
                            print "Invalid data in row ",i," <",row,"> - sixth field must be comment with #"
                            raise Exception(str("Invalid data in %s"%args['input']))

                    angleind        = int(row[0])
                    end_x           = int(row[1])
                    end_y           = int(row[2])
                    angle_increment = int(row[3])
                    costmult        = int(row[4])

                    if (angleind >= args['numberofanglesperquadrant']):
                        print "Invalid angle ID=",angleind," >= ",args['numberofanglesperquadrant']," first quadrant angles!"
                        raise Exception(" Invalid number of angles defined in this file!")

                    base_mprim_end_points[angleind].append(np.array((end_x, end_y, angle_increment, costmult)))
                    print "  base_mprim_end_points[",angleind,"]","[",len(base_mprim_end_points[angleind])-1,"]",base_mprim_end_points[angleind][-1]
                elif (not getAngles) and (not getPrims):
                    if (row[0] == "Angles"):
                        getAngles = True
                    elif (row[0] == "Resolution"):
                        print "Reset resolution from {:.3f} to {:.3f}".format(args['resolution'],float(row[1]))
                        args['resolution'] = float(row[1])
                    else:
                        print "Invalid row searching for Angles keyword"
                        print row

                else:
                        print "Invalid processing for row=",row,">"
                        raise Exception(" Invalid processing for row=%s"%(row))

    print "Finished reading input file!"

    consistent = True

    # Check that the number of motion primitives is consistent with parameters
    mprims_per_angle = 0
    for angle, angle_def in base_mprim_end_points.iteritems():
        mprim_size = len(angle_def)
        if (mprims_per_angle > 0 and mprim_size != mprims_per_angle):
            print "Inconsistent number of motion primitives for angle=",angle," mprims=",mprim_size
            consistent = False

        if (mprim_size > mprims_per_angle):
            mprims_per_angle = mprim_size

    if (args['numberofprimsperangle'] != mprims_per_angle):
        print "Number of loaded motion primitives is not consistent with parameters ",mprims_per_angle," != ",args['numberofprimsperangle']
        args['numberofprimsperangle'] = mprims_per_angle
        consitent = False

    if not consistent:
        raise Exception(" Inconsistent primitives!")


    return base_mprim_end_points

def get_primitive_generators(args):
    """
        Get the primitive generators from defaults or input file
    """
    primitive_generator_definitions = None
    if (args['input'] != ""):
        print "Load primitive generator definitions from file <",args['input'],"> ..."
        primitive_generator_definitions = read_primitive_generator_definitions(args)
    else:
        print "Use default primitive generators ..."
        primitive_generator_definitions = get_default_primitive_generators(args)

    #print "primitive_generator_definitions=",primitive_generator_definitions
    return primitive_generator_definitions
