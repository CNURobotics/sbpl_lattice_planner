
#% /*
#%  * Copyright (c) 2016, David Conner (Christopher Newport University)
#%  *
#%  * Based on genmprim_unicycle.m  Copyright (c) 2008, Maxim Likhachev
#%  * All rights reserved.
#%  * 
#%  * Redistribution and use in source and binary forms, with or without
#%  * modification, are permitted provided that the following conditions are met:
#%  * 
#%  *     * Redistributions of source code must retain the above copyright
#%  *       notice, this list of conditions and the following disclaimer.
#%  *     * Redistributions in binary form must reproduce the above copyright
#%  *       notice, this list of conditions and the following disclaimer in the
#%  *       documentation and/or other materials provided with the distribution.
#%  *     * Neither the name of the Carnegie Mellon University or
#%  *       Christopher Newport University nor the names of its
#%  *       contributors may be used to endorse or promote products derived from
#%  *       this software without specific prior written permission.
#%  * 
#%  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#%  * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#%  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#%  * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#%  * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#%  * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#%  * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#%  * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#%  * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#%  * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#%  * POSSIBILITY OF SUCH DAMAGE.
#%  */

import sys#, types, time, random, os
import numpy as np

# if available import pylab (from matlibplot)
visualize_plt = False
try:
    import matplotlib.pylab as plt
    visualize_plt = True
except ImportError:
    pass

def matrixSize(mat, elem=None):
    if elem is None:
        return mat.shape
    else:
        return mat.shape[elem]


# Get the primitive generators for the first quadrant (assume symmetric for other for quadrants)
def getDefaultPrimitiveGenerators(args):
    base_mprim_end_points = dict()


    if (args['numberofprimsperangle'] != 10):
        print "default primitives use 10 primitives per angle!"
        args['numberofprimsperangle'] = 10


    numberofangles          = args['numberofangles']
    numberofprimsperangle   = args['numberofprimsperangle']
    
    if (numberofangles%4 != 0):
        raise Exception(" Invalid number of angles - not even by quadrants!")

    increments = numberofangles / 4

    if (increments != 4):
        raise Exception(" This file assumes 4 increments for first quadrant!")
    
    # Cost multipliers (multiplier is used as costmult*cost)
    forwardcostmult         = args['forwardcostmult']
    forwardandturncostmult  = args['forwardandturncostmult']
    backwardcostmult        = args['backwardcostmult']
    turninplacecostmult     = args['turninplacecostmult']
    sidestepcostmult        = args['sidestepcostmult']
    


    # Define matrices to hold primitive data
    for ang in range(0,increments):                             #%x,y,theta,costmult 
        base_mprim_end_points[ang] = np.zeros((numberofprimsperangle, 4))

    # Define the end points for each angle in the first quadrant
    #%x aligned with the heading of the robot, angles are positive counter clockwise
    #%note, what is shown x,y,theta changes (not absolute numbers)


    # Aligned 0 degrees 
    # 0 theta change
    base_mprim_end_points[0][0,:] = np.array(np.hstack(( 1.,  0.,  0., forwardcostmult)))
    base_mprim_end_points[0][1,:] = np.array(np.hstack((10.,  0.,  0., forwardcostmult)))
    base_mprim_end_points[0][2,:] = np.array(np.hstack(( 5.,  0.,  0., forwardcostmult)))
    base_mprim_end_points[0][3,:] = np.array(np.hstack((-1.,  0.,  0., backwardcostmult)))
    #%1/16 theta change
    base_mprim_end_points[0][4,:] = np.array(np.hstack(( 6.,  1.,  1., forwardandturncostmult)))
    base_mprim_end_points[0][5,:] = np.array(np.hstack(( 6., -1., -1., forwardandturncostmult)))
    base_mprim_end_points[0][6,:] = np.array(np.hstack(( 9.,  1.,  1., forwardandturncostmult)))
    base_mprim_end_points[0][7,:] = np.array(np.hstack(( 9., -1., -1., forwardandturncostmult)))
    #%turn in place
    base_mprim_end_points[0][8,:] = np.array(np.hstack(( 0.,  0.,  1., turninplacecostmult)))
    base_mprim_end_points[0][9,:] = np.array(np.hstack(( 0.,  0., -1., turninplacecostmult)))

    #Aligned to 22.5 degrees on grid
    #%0 theta change     
    base_mprim_end_points[1][0,:] = np.array(np.hstack(( 2.,  1.,  0., forwardcostmult)))
    base_mprim_end_points[1][1,:] = np.array(np.hstack(( 8.,  4.,  0., forwardcostmult)))
    base_mprim_end_points[1][2,:] = np.array(np.hstack(( 4.,  2.,  0., forwardcostmult)))
    base_mprim_end_points[1][3,:] = np.array(np.hstack((-2., -1.,  0., backwardcostmult)))
    #%1/16 theta change
    base_mprim_end_points[1][4,:] = np.array(np.hstack(( 5.,  4.,  1., forwardandturncostmult)))
    base_mprim_end_points[1][5,:] = np.array(np.hstack(( 6.,  1., -1., forwardandturncostmult)))
    base_mprim_end_points[1][6,:] = np.array(np.hstack(( 8.,  5.,  1., forwardandturncostmult)))
    base_mprim_end_points[1][7,:] = np.array(np.hstack(( 8.,  3., -1., forwardandturncostmult)))
    #%turn in place
    base_mprim_end_points[1][8,:] = np.array(np.hstack((0., 0.,  1., turninplacecostmult)))
    base_mprim_end_points[1][9,:] = np.array(np.hstack((0., 0., -1., turninplacecostmult)))

    # Aligned to 45 degrees on grid
    #%0 theta change 
    base_mprim_end_points[2][0,:] = np.array(np.hstack(( 1.,  1.,  0., forwardcostmult)))
    base_mprim_end_points[2][1,:] = np.array(np.hstack(( 7.,  7.,  0., forwardcostmult)))
    base_mprim_end_points[2][2,:] = np.array(np.hstack(( 4.,  4.,  0., forwardcostmult)))
    base_mprim_end_points[2][3,:] = np.array(np.hstack((-1., -1.,  0., backwardcostmult)))
    #%1/16 theta change
    base_mprim_end_points[2][4,:] = np.array(np.hstack(( 4.,  6.,  1., forwardandturncostmult)))
    base_mprim_end_points[2][5,:] = np.array(np.hstack(( 6.,  4., -1., forwardandturncostmult)))
    base_mprim_end_points[2][6,:] = np.array(np.hstack(( 6.,  8.,  1., forwardandturncostmult)))
    base_mprim_end_points[2][7,:] = np.array(np.hstack(( 8.,  6., -1., forwardandturncostmult)))
    #%turn in place
    base_mprim_end_points[2][8,:] = np.array(np.hstack(( 0.,  0.,  1., turninplacecostmult)))
    base_mprim_end_points[2][9,:] = np.array(np.hstack(( 0.,  0., -1., turninplacecostmult)))

    # Aligned to 67.5 degrees on grid 
    for primind in range(0,numberofprimsperangle):
        # reverse x and y, and negate angle
        base_mprim_end_points[3][primind,0] =  base_mprim_end_points[1][primind,1]
        base_mprim_end_points[3][primind,1] =  base_mprim_end_points[1][primind,0]
        base_mprim_end_points[3][primind,2] = -base_mprim_end_points[1][primind,2]
        base_mprim_end_points[3][primind,3] =  base_mprim_end_points[1][primind,3]

    return base_mprim_end_points

# Read in the primitive generator definitions from a CSV text file
def readPrimitiveGeneratorDefinitions(args):
    raise Exception(" Not defined yet!")
    base_mprim_end_points = dict()

    return base_mprim_end_points


# Get the primitive generator and the angle of rotation based on the designated quadrant
def getMotionPrimitiveEndpoint(angleind, primind, currentangle, primitive_generator_definitions,args):
    resolution              = args['resolution']
    numberofangles          = args['numberofangles']
    numberofprimsperangle   = args['numberofprimsperangle']

    quadrant = angleind/4
    generatorIndex = angleind%4

    # Return the rotation angle based on quadrant and the specific primitive generator
    if (quadrant == 0):
        # First quadrant 0-90
        return (0., primitive_generator_definitions[generatorIndex][primind,:])

    elif (quadrant == 1):
        # Second quadrant 90-180
        return (np.pi/2.0, primitive_generator_definitions[generatorIndex][primind,:])

    elif (quadrant == 2):
        # Third quadrant 180-270
        return (np.pi, primitive_generator_definitions[generatorIndex][primind,:])
    elif (quadrant == 3):
        # Fourth quadrant 270-360
        return (3.0*np.pi/2.0, primitive_generator_definitions[generatorIndex][primind,:])
    else:
        raise Exception(str("Invalid quadrant =  %d"%quadrant))

    return None

def writeMotionPrimitiveFile(primitive_definitions,args):
    resolution              = args['resolution']
    numberofangles          = args['numberofangles']
    numberofprimsperangle   = args['numberofprimsperangle']
    
    print " Write primitives to file <",args['output'],"> ..."
    fout = open(args['output'], 'w')

    #%write the header
    fout.write('resolution_m: %f\n'%(resolution))
    fout.write('numberofangles: %d\n'%(numberofangles))
    fout.write('totalnumberofprimitives: %d\n'%(numberofprimsperangle*numberofangles))

    #%iterate over angles
    for angleind in range(0, numberofangles):
        currentangle = (angleind*2.*np.pi)/numberofangles
        currentangle_36000int = np.round(angleind*36000./numberofangles)

        #%iterate over primitives    
        for primind in range(0, numberofprimsperangle):
            fout.write('primID: %d\n'%primind)
            fout.write('startangle_c: %d\n'%angleind)

            prim_def_primitive = primitive_definitions[angleind][primind]

            endpose_c      = prim_def_primitive['endpose']
            intermcells_m  = prim_def_primitive['intermcells']
            actioncostmult = prim_def_primitive['actioncostmult']

            #%write out
            fout.write('endpose_c: %d %d %d\n'% ( endpose_c[0], endpose_c[1], endpose_c[2]))
            fout.write('additionalactioncostmult: %d\n'%( actioncostmult))
            fout.write('intermediateposes: %d\n'%( matrixSize(intermcells_m, 0)))
            for interind in range(0, matrixSize(intermcells_m, 0)):
                fout.write('%.4f %.4f %.4f\n'%( intermcells_m[interind,0], 
                                                intermcells_m[interind,1], 
                                                intermcells_m[interind,2]))

    fout.close()
    print "     Finished writing the motion primitive file!"
    return

def visualizeMotionPrimitives(primitive_definitions,args):
    resolution              = args['resolution']
    numberofangles          = args['numberofangles']
    numberofprimsperangle   = args['numberofprimsperangle']
    
    #%iterate over angles
    for angleind in range(0, numberofangles):
        currentangle = (angleind*2.*np.pi)/numberofangles
        currentangle_36000int = np.round(angleind*36000./numberofangles)

        plt.figure(angleind)  # use angleind to plot each primitive set on different window
        plt.hold(True)
        plt.grid(True)
        plt.axis('equal')
        plt.axis([-10*resolution,10*resolution,-10*resolution,10*resolution])
        plt.title(str(angleind)+" "+str(currentangle_36000int/100.))

        #%iterate over primitives    
        for primind in range(0, numberofprimsperangle):
            prim_def_primitive = primitive_definitions[angleind][primind]
            endpose_c      = prim_def_primitive['endpose']
            intermcells_m  = prim_def_primitive['intermcells']
            actioncostmult = prim_def_primitive['actioncostmult']
            endpt          = prim_def_primitive['endpoint']
            plt.plot(intermcells_m[:,0], intermcells_m[:,1],linestyle="-",marker="o")
            plt.text(endpt[0], endpt[1], str(primind)+" ("+str(endpose_c[2])+")")

        plt.waitforbuttonpress()  # uncomment to plot each primitive set one at a time
        
    #    plt.waitforbuttonpress()  # Uncomment to hold until buttom pressed
    print "Hold windows open until all closed"
    plt.show()  # Keep windows open until the program is terminated
    return

def generateMotionPrimitives(primitive_generator_definitions,args):

    print "Generate motion primitive definitions ..."
    primitive_definitions = [] # List to hold data for each angle index

    resolution              = args['resolution']
    numberofangles          = args['numberofangles']
    numberofprimsperangle   = args['numberofprimsperangle']
    numofsamples            = args['numberOfSamples']

     #%iterate over angles
    for angleind in range(0, numberofangles):
        primitive_definitions.append([]) # list to hold data to hold data for each primitive within angle

        prim_def_angle_list = primitive_definitions[angleind] # reference to angle list

        currentangle = (angleind*2.*np.pi)/numberofangles
        currentangle_36000int = np.round(angleind*36000./numberofangles)

        #%iterate over primitives    
        for primind in range(0, numberofprimsperangle):

            angle, basemprimendpts_c = getMotionPrimitiveEndpoint(angleind, primind, currentangle, primitive_generator_definitions,args)

            # Define action to terminal end point     
            baseendpose_c = basemprimendpts_c[0:3]
            endx_c = np.round((baseendpose_c[0]*np.cos(angle))-(baseendpose_c[1]*np.sin(angle)))
            endy_c = np.round((baseendpose_c[0]*np.sin(angle))+(baseendpose_c[1]*np.cos(angle)))
            endtheta_c = (angleind+baseendpose_c[2])%numberofangles
            endpose_c = np.array(np.hstack((endx_c, endy_c, endtheta_c)))


            startpt = np.array(np.hstack((0., 0., currentangle)))
            endpt   = np.array(np.hstack(((endpose_c[0]*resolution), 
                                          (endpose_c[1]*resolution), 
                                          (( (angleind+baseendpose_c[2])%numberofangles*2.*np.pi)/numberofangles) ) ) )
            

            #print "\n-----------------\nangleind=",angleind,"  primind=",primind
            #print "endpose_c=",endpose_c            
            #print( 'rotation angle=%f\n'% (angle*180./np.pi))
            #print "startpt =",startpt
            #print "endpt   =",endpt


            # Generate intermediate poses (remember they are w.r.t 0,0 (and not centers of the cells)
            intermcells_m = np.zeros((numofsamples, 3))

            # Store the specific motion primitive data in dictionary
            prim_def_angle_list.append([])                    # Add this primitive to angle list
            prim_def_angle_list[primind] = dict()             # Define a dictionary to hold the data
            prim_def_primitive = prim_def_angle_list[primind] # reference to specific primitive dictionary
            prim_def_primitive['angle']            = angle
            prim_def_primitive['basemprim_end_pt'] = basemprimendpts_c
            prim_def_primitive['endpose']          = endpose_c
            prim_def_primitive['intermcells']      = intermcells_m
            prim_def_primitive['startpoint']       = startpt
            prim_def_primitive['endpoint']         = endpt 
            prim_def_primitive['actioncostmult']   = basemprimendpts_c[3]

            if np.logical_or(np.logical_and(endx_c == 0., endy_c == 0.), baseendpose_c[2] == 0.):
                #%turn in place or move forward            
                for iind in range(0, numofsamples):
                    fraction = float(iind+1)/(numofsamples)                 
                    intermcells_m[iind,:] = np.array((startpt[0]+(endpt[0]-startpt[0])*fraction, 
                                                      startpt[1]+(endpt[1]-startpt[1])*fraction, 0))
                    rotation_angle = baseendpose_c[2]*(2.*np.pi/numberofangles)

                    intermcells_m[iind,2] = (startpt[2]+rotation_angle*fraction)%(2.*np.pi)
                    #print " ",iind,"  of ",numofsamples," fraction=",fraction," rotation=",rotation_angle
                
            else:
                #%unicycle-based move forward or backward  (http://sbpl.net/node/53) 
                R = np.array(np.vstack((np.hstack((np.cos(startpt[2]),  np.sin(endpt[2])-np.sin(startpt[2]))), 
                                        np.hstack((np.sin(startpt[2]), -np.cos(endpt[2])+np.cos(startpt[2]))) )))

                S = np.dot(np.linalg.pinv(R), np.array(np.vstack((endpt[0]-startpt[0], endpt[1]-startpt[1])) ))
                l = S[0]       # l in www.sbpl.net/node/53  eqn (10)
                radius = S[1]  # r in eqn (10)

                omega = (baseendpose_c[2]*2.*np.pi/numberofangles)+ l/radius # w in eqn (12) www.sbpl.net/node/53 (unit time scaling)
                vel   = radius*omega # v in eqn (12) www.sbpl.net/node/53  (in terms of unit time scaling)

                #print "R=\n",R
                #print "Rpi=\n",np.linalg.pinv(R)
                #print "S=\n",S
                #print "l=",l
                #print "radius=",radius

                if l<0.:
                    print('WARNING: l = %f < 0 -> bad action start/end points - no straight section\n'%(l))
                    l = 0.
                    chord2 = np.array((endpt[0]-startpt[0], endpt[1]-startpt[1]))
                    chord  = np.sqrt(np.dot(chord2,chord2)) # Chord length
                    dangle = baseendpose_c[2]*2.*np.pi/numberofangles
                    radius = (chord/2.0)/np.sin(dangle/2.0);
                    #print "     new radius = ",radius
                    omega  = dangle
                    vel    = omega*radius

                
                # Scaled time for linear section (Eqn (12) in www.sbpl.net/node/53)
                tl  = l/vel

                #print "omega=",omega
                #print "vel=",vel
                #print "scaled tl=",tl
                
                #%compute rv
                #%rv = baseendpose_c(3)*2*pi/numberofangles;
                #%compute tv
                #%tvx = (endpt(1) - startpt(1))*rv/(sin(endpt(3)) - sin(startpt(3)))
                #%tvy = -(endpt(2) - startpt(2))*rv/(cos(endpt(3)) - cos(startpt(3)))
                #%tv = (tvx + tvy)/2.0;              
                #%generate samples
                for iind in range(0, numofsamples):
                    dt = float(iind+1)/(numofsamples)
                    #%dtheta = rv*dt + startpt(3);
                    #%intermcells_m(iind,:) = [startpt(1) + tv/rv*(sin(dtheta) - sin(startpt(3))) ...
                    #%                        startpt(2) - tv/rv*(cos(dtheta) - cos(startpt(3))) ...
                    #%                        dtheta];
                    if (dt*vel)<l:
                        # Eqn (2) in www.sbpl.net/node/53 
                        intermcells_m[iind,:] = np.array(np.hstack((startpt[0]+dt*vel*np.cos(startpt[2]), 
                                                                    startpt[1]+dt*vel*np.sin(startpt[2]), 
                                                                    startpt[2])))

                        first_arc = 1.*intermcells_m[iind,:]
                    else:
                        # Eqn (9) in www.sbpl.net/node/53 
                        dtheta = omega*(dt-tl) + startpt[2]
                        intermcells_m[iind,:] = np.array(np.hstack((startpt[0] + l*np.cos(startpt[2]) + radius*(np.sin(dtheta)-np.sin(startpt[2])), 
                                                                    startpt[1] + l*np.sin(startpt[2]) - radius*(np.cos(dtheta)-np.cos(startpt[2])), 
                                                                    dtheta)))
                        
                #%correct
                errorxy = np.array(np.hstack((endpt[0] -intermcells_m[int(numofsamples)-1,0], 
                                              endpt[1] -intermcells_m[int(numofsamples)-1,1])))
                #print('l=%f errx=%f erry=%f\n'%(l, errorxy[0], errorxy[1]))
                interpfactor = np.array(np.hstack((np.arange(0., 1.+(1./(numofsamples)), 1./(numofsamples-1)))))

                #print "intermcells_m=",intermcells_m
                #print "interp'=",interpfactor.conj().T

                intermcells_m[:,0] = intermcells_m[:,0]+errorxy[0]*interpfactor.conj().T
                intermcells_m[:,1] = intermcells_m[:,1]+errorxy[1]*interpfactor.conj().T

    print "Finished primitive definitions!"
    return primitive_definitions

def defineMotionPrimitivesExtended(args):

    visualize = visualize_plt and args['showPrimitives']  # Plot the primitives

    # Arguments
    resolution              = args['resolution']
    numberofangles          = args['numberofangles']
    numberofprimsperangle   = args['numberofprimsperangle']
    
    # Cost multipliers (multiplier is used as costmult*cost)
    forwardcostmult         = args['forwardcostmult']
    forwardandturncostmult  = args['forwardandturncostmult']
    backwardcostmult        = args['backwardcostmult']
    turninplacecostmult     = args['turninplacecostmult']
    sidestepcostmult        = args['sidestepcostmult']


    primitive_generator_definitions = None

    if (args['input'] != ""):
        print "Load primitive generator definitions from file <",args['input'],"> ..."
        primitive_generator_definitions = readPrimitiveGeneratorsDefinitions(args)
    else:
        print "Use default primitives ..."
        primitive_generator_definitions = getDefaultPrimitiveGenerators(args)

    #print "primitive_generator_definitions=",primitive_generator_definitions

    primitive_definitions = generateMotionPrimitives(primitive_generator_definitions, args)

    if (args['output'] != ""):
        try:
            writeMotionPrimitiveFile(primitive_definitions,args)
        except OSError as err:
            print "Failed to write primitives to file!"
            print("OS error: {0}".format(err))
        except ValueError as err:
            print "Failed to write primitives to file!"
            print("Value error: {0}".format(err))
        except TypeError as err:
            print "Failed to write primitives to file!"
            print("Type error: {0}".format(err))
        except NameError as err:
            print "Failed to write primitives to file!"
            print("Name error: {0}".format(err))
        except:
            print "Failed to write primitives to file!"
            print("Unexpected error:", sys.exc_info()[0])            

    if (args['showPrimitives'] ):
        try:
            visualizeMotionPrimitives(primitive_definitions,args)
        except OSError as err:
            print "Failed to write primitives to file!"
            print("OS error: {0}".format(err))
        except NameError as err:
            print "Failed to write primitives to file!"
            print("Name error: {0}".format(err))
        except:
            print "Cannot visualize the primitives"            
            print("     Unexpected error:", sys.exc_info()[0])            

          
    print "Finished motion primitive generation!"      
    return


# Helper for printing default arguments
def default(str):
    return str + ' [Default: %default]'

# Parse the command line arguments
def readCommand( argv ):
    """
    Processes the command used to run genmprim_extended from the command line.
    """
    from optparse import OptionParser

    usageStr = """
    USAGE:      python genmprim_extended.py <options>
                     - creates the designated motion primitives

    """
    parser = OptionParser(usageStr)

    parser.add_option('-o', '--output', dest='output', type='string',
                      help=default('The output file name (including path information)'), 
                      metavar='FILE', default="extended.mprim")
    parser.add_option('-i', '--input', dest='input', type='string',
                      help=default('Optional input file name with primitive terminal grid points'), 
                      metavar='PRIM_DEFINITION_FILE', default="")
    parser.add_option('-s', '--showPrimitives', action='store_true', dest='showPrimitives',
                      help='Show the generated primitives', default=False)
    parser.add_option('-g', '--showPotentialGrid', action='store_true', dest='showPotential',
                      help='Show the valid grid values', default=False)
    parser.add_option('-r', '--resolution', type='float', dest='resolution',
                      help=default('Translational resolution of the lattice grid'), default=0.025)
    parser.add_option('-m', '--minTurnRadius', type='float', dest='minTurnRadius',
                      help=default('Minimum forward turning radius (same units as resolution)'), default=0.1)
    parser.add_option('-n', '--numberOfAngles', type='int', dest='numberOfAngles',
                      help=default('Number of angles in the lattice grid (preferably a power of 2, definitely multiple of 8)'), 
                      metavar='ANGLES',  default=16)
    parser.add_option('-p', '--primitivesPerAngle', type='int', dest='primPerAngle',
                      help=default('Number of primitives per angle in the lattice grid'),  default=16)
    parser.add_option('-f', '--fwdCostMult', type='int', dest='fwdCost',
                      help=default('Forward cost multiplier'), 
                      metavar='FWDCOST',  default=1)
    parser.add_option('-t', '--fwdTurnCostMult', type='int', dest='fwdTurnCost',
                      help=default('Forward cost multiplier'),  default=2)
    parser.add_option('-b', '--backCostMult', type='int', dest='backCost',
                      help=default('Backward cost multiplier'),  default=5)
    parser.add_option('-z', '--turnInPlaceCostMult', type='int', dest='turnInPlaceCost',
                      help=default('Turn in place (zero radius turn) cost multiplier'), default=5)
    parser.add_option('-j', '--sideStepCostMult', type='int', dest='sideStepCost',
                      help=default('Side step (jump) cost multiplier'), default=10)
    parser.add_option('-d', '--numberOfSamples', type='int', dest='numSamples',
                      help=default('Number of samples along motion primitive'), default=10)
     
    options, otherjunk = parser.parse_args(argv)
    if len(otherjunk) != 0:
        raise Exception('Command line input not understood: ' + str(otherjunk))

    # Create a dictionary of arguments
    args = dict()


    # Choose a layout
    args['output'] = options.output
    if args['output'] == None or args['output'] == "": 
        raise Exception("The output file is not defined !")

    args['input'] = options.input

    # Lattice definition
    args['resolution']              = options.resolution
    args['numberofangles']          = options.numberOfAngles
    args['numberofprimsperangle']   = options.primPerAngle
    args['minTurnRadius']           = options.minTurnRadius
    args['numberOfSamples']         = options.numSamples

    # Cost multipliers (multiplier is used as costmult*cost)
    args['forwardcostmult']         = options.fwdCost
    args['forwardandturncostmult']  = options.fwdTurnCost
    args['backwardcostmult']        = options.backCost
    args['turninplacecostmult']     = options.turnInPlaceCost
    args['sidestepcostmult']        = options.sideStepCost
    
    # Display options
    args['showPrimitives'] = options.showPrimitives
    args['showPotential'] = options.showPotential

    if (options.showPrimitives and not visualize_plt):
        raise Exception("Cannot show primitives because MatPlotLib is not available!")

    if (options.showPotential and not visualize_plt):
        raise Exception("Cannot show potential grid points because MatPlotLib is not available!")

    if (options.numberOfAngles%8 != 0):
        raise Exception(str("Construction requires that number of angles is a multiple of 8 (Preferably a power of 2)  numAngles=%d" % (options.numberOfAngles)))

    return args

if __name__ == '__main__':
    """
    The main function called when genmprim_extended.py is run from the command line:

    > python genmprim_extended.py

    See the usage string for more details.

    > python genmprim_extended.py --help
    """
    args = readCommand( sys.argv[1:] ) # Get options from command line
    print "Running genmprim_extended with the following arguments:"
    print args

    if (args['showPotential']):
        print "Display acceptable end points in lattice given parameters"
    else:
        print "Define the motion primitives give generator definitions"
        defineMotionPrimitivesExtended(args)
