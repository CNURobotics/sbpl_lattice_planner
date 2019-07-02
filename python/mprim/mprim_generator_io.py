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

def write_motion_primitives_file(primitive_definitions,args):
    """
        Write the full set of primitives to a file
    """
    try:
      resolution              = args['resolution']
      numberofangles          = args['numberofanglesperquadrant']*4
      numberofprimsperangle   = args['numberofprimsperangle']
      min_turning_radius_m    = args['minTurnRadius']

      print " Write primitives to file <",args['output'],"> ..."
      with open(args['output'], 'w') as fout:

        # write the header
        fout.write('resolution_m: %f\n'%(resolution))
        fout.write('min_turning_radius_m: %f\n'%(min_turning_radius_m))
        fout.write('numberofangles: %d\n'%(numberofangles))
        for ndx in range(0,numberofangles):
            fout.write('angle:%d %f\n'%(ndx,args['angles'][ndx]))

        fout.write('totalnumberofprimitives: %d\n'%(numberofprimsperangle*numberofangles))

        # iterate over angles
        for angleind in range(0, numberofangles):
            currentangle = (angleind*2.*np.pi)/numberofangles

            # iterate over primitives
            for primind in range(0, numberofprimsperangle):
              try:
                fout.write('primID: %d\n'%primind)
                fout.write('startangle_c: %d\n'%angleind)

                prim_def_primitive = primitive_definitions[angleind][primind]

                endpose_c      = prim_def_primitive['endpose']
                intermcells_m  = prim_def_primitive['intermcells']
                actioncostmult = prim_def_primitive['actioncostmult']
                actioncost     = prim_def_primitive['actioncost']
                turning_radius = prim_def_primitive['turning_radius']

                # write out
                fout.write('endpose_c: %d %d %d\n'% ( endpose_c[0], endpose_c[1], endpose_c[2]))
                fout.write('actioncost: %d\n'%( actioncost))
                fout.write('additionalactioncostmult: %d\n'%( actioncostmult))
                fout.write('turning_radius: %f\n' %( turning_radius))
                #fout.write('# %s\n'% prim_def_primitive['style'])
                fout.write('intermediateposes: %d\n'%( matrix_size(intermcells_m, 0)))
                for interind in range(0, matrix_size(intermcells_m, 0)):
                    fout.write('%.4f %.4f %.4f\n'%( intermcells_m[interind,0],
                                                    intermcells_m[interind,1],
                                                    intermcells_m[interind,2]))
              except Exception as err:
                print "Failed to write primitive {} {} to file!".format(angleind,primind)
                print "  Error: {}".format(err)

      print "     Finished writing the motion primitive file!"
    except OSError as err:
        print "Failed to write primitives to file!"
        print("OS error: {0}".format(err))
    except ValueError as err:
        print "Failed to write primitives to file!"
        print("Value error: {0}".format(err))
    except TypeError as err:
        print "Failed to write primitives to file!"
        print("Type error: {0}".format(err))
    except KeyError as err:
        print "Failed to write primitives to file!"
        print("Key error: {0}".format(err))
    except NameError as err:
        print "Failed to write primitives to file!"
        print("Name error: {0}".format(err))
    except:
        print "Failed to write primitives to file!"
        print("Unexpected error:", sys.exc_info()[0])

    return

def read_motion_primitives_file(filename, primitive_definitions,args):
    """
        Read the full set of primitives from a file
    """
    initialized = False
    initial_data = (
            ('resolution_m:', 'resolution'),
            ('min_turning_radius_m:','minTurnRadius'),
            ('numberofangles:', 'numberofangles'),
            ('angle:','angles'),
            ('totalnumberofprimitives:','numberofprimsperangle'))
    initial_cnt = 0

    initial_data_dict = {}
    for tup in initial_data:
        initial_data_dict[tup[0]] = tup[1]

    prim_data_text = ('primID:','startangle_c:',
                      'endpose_c:', 'actioncost:',
                      'additionalactioncostmult:',
                      'turning_radius:',
                      'intermediateposes:')

    prim_data_cnt = 0
    interm_cnt = 0
    intermcells_m = []
    expected_prims = None
    total_prims = 0
    angle_cnt = 0
    angles = []
    try:
      with open(filename,'rt') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip()
            line_vals = line.split(' ')
            if not initialized:
                if not line.startswith(initial_data[initial_cnt][0]):
                    # Mismatch - old style file?
                    if line_vals[0] in initial_data_dict:
                        # The  line is expected laster, so
                        # process missing data
                        if initial_data[initial_cnt][0].startswith('angle:'):
                            # The expected line is angle:, so we must
                            # be using regular angle spacing
                            print "Using fixed number of angles=",args['numberofangles']
                            angles=[0.0]
                            ang_increment = 2.0*math.pi/args['numberofangles']

                            for i in range(1,args['numberofangles']):
                                angles.append(angles[-1]+ang_increment)
                            if len(angles) == args['numberofangles']:
                                args['numberofanglesperquadrant'] = len(angles)//4
                                args[initial_data[initial_cnt][1]] = angles
                            initial_cnt += 1
                        elif initial_data[initial_cnt][0].startswith('min_turning_radius_m'):
                            print "Old style primitive - min radius = 0.0"
                            args['minTurnRadius'] = 0.0
                            initial_cnt += 1
                    else:
                        print "Unknown match at ",initial_cnt
                        print " line_vals:",line_vals
                        print " keys:",initial_data_dict.keys()
                        print "  Prior: >",initial_data[initial_cnt][0],"<"
                        print "  Line : >",line,"<"
                        print "Current: >",initial_data[initial_cnt+1][0],"<"
                        print "unknown primitive style at line=",initial_cnt
                        raise Exception("Unknown primitive style!")


                # Process expected data
                if line.startswith(initial_data[initial_cnt][0]):
                    key = initial_data[initial_cnt][1]
                    len_text = len(initial_data[initial_cnt][0])
                    text = line[len_text:].strip()
                    if initial_cnt < 2:
                        args[key] = float(text)
                        print "Read ",key," = ",args[key]
                        initial_cnt += 1
                    elif initial_cnt == 2:
                        args[key] = int(text)
                        print "Read ",key," = ",args[key]
                        initial_cnt += 1
                    elif initial_cnt == 3:
                        text = text.split(' ')
                        if int(text[0]) != len(angles):
                            print line
                            print text
                            print angles
                            raise Exception(" Invalid angle !")
                        else:
                            angles.append(float(text[1]))

                        if len(angles) == args['numberofangles']:
                            args['numberofanglesperquadrant'] = len(angles)//4
                            args[key] = angles
                            initial_cnt += 1
                    else:
                        print 30*"="," Final check of initialization ",30*"="
                        print args

                        expected_prims = int(text)
                        initial_cnt += 1
                        assert initial_cnt == len(initial_data),"Failed to initialize properly"
                        numberofangles = args['numberofanglesperquadrant']*4
                        args['numberofprimsperangle'] = expected_prims//numberofangles

                        assert expected_prims == args['numberofprimsperangle']*numberofangles, "Failed to load initial primitive data properly : {} {} {}".format(key, expected_prims, numberofprimsperangle*numberofangles)
                        initialized = True
                        prim_data_cnt = 0
                        interm_cnt = 0
                        intermcells_m = []
                        primind=-1
                        angleind=-1
                else:
                    print "Unknown match at ",initial_cnt
                    print "  Prior: >",initial_data[initial_cnt][0],"<"
                    print "  Line : >",line,"<"
                    print "Current: >",initial_data[initial_cnt+1][0],"<"
                    print "unknown primitive style at line=",initial_cnt
                    raise Exception("Unknown primitive style!")

            else:
                # Now we are initialized, need to read primitives
                if prim_data_cnt < len(prim_data_text):
                    if not line.startswith(prim_data_text[prim_data_cnt]):
                        # Check for old style primitive with missing data
                        if prim_data_text[prim_data_cnt] == 'actioncost:' and \
                            line.startswith(prim_data_text[prim_data_cnt+1]):
                            print "Missing action cost for ",primind, angleind
                            actioncost = 1
                            prim_data_cnt += 1


                        if prim_data_text[prim_data_cnt] == 'turning_radius:' and \
                            line.startswith(prim_data_text[prim_data_cnt+1]):
                            print "Missing turning radius for ",primind, angleind
                            turning_radius = None
                            prim_data_cnt += 1


                    if line.startswith(prim_data_text[prim_data_cnt]):
                        len_text = len(prim_data_text[prim_data_cnt])
                        text = line[len_text:].strip()
                        if prim_data_cnt == 0:
                            primind = int(text)
                        if prim_data_cnt == 1:
                            angleind = int(text)
                        if prim_data_cnt == 2:
                            endpose_c = np.array([float(val) for val in text.split(' ')])
                        if prim_data_cnt == 3:
                            actioncost = int(text)
                        if prim_data_cnt == 4:
                            actioncostmult = int(text)
                        if prim_data_cnt == 5:
                            turning_radius = float(text)
                        if prim_data_cnt == 6:
                            interm_cnt = int(text)
                            intermcells_m = []

                        prim_data_cnt += 1
                    else:
                        print "not starts with ",prim_data_text[prim_data_cnt], "<",line,">"
                        prim_data_cnt += 1


                else:

                    intermcells_m.append([float(val) for val in line.split(' ')])

                    if len(intermcells_m) == interm_cnt:
                        print "Complete m-prim for ",primind,angleind

                        prim_def_primitive={}
                        prim_def_primitive['endpose'] = endpose_c
                        prim_def_primitive['intermcells'] = np.array(intermcells_m)
                        prim_def_primitive['endpoint'] = intermcells_m[-1]
                        prim_def_primitive['actioncost']=actioncost
                        prim_def_primitive['actioncostmult']=actioncostmult
                        prim_def_primitive['turning_radius']=turning_radius
                        #prim_def_primitive['style'] = style
                        print "set the primitive for ",primind, angleind
                        if angleind == len(primitive_definitions):
                            primitive_definitions.append([])

                        if primind == len(primitive_definitions[angleind]):
                            primitive_definitions[angleind].append(prim_def_primitive)
                            prim_data_cnt = 0
                            interm_cnt = 0
                            intermcells_m = []
                            total_prims += 1

                            # Reset data for next primitive
                            primind = None
                            angleind = None
                            endpose_c = None
                            actioncost = None
                            actioncostmult = None
                            turning_radius = None
                            interm_cnt = None
                        else:
                            print "Angle, prim",angleind, primind
                            print line
                            print primitive_definitions

        assert total_prims == expected_prims," Invalid number of prims loaded! {} vs. {}".format(total_prims, expected_prims)
        print "     Finished reading the motion primitive file!"
        return True
    except OSError as err:
        print "Failed to read primitives from file!"
        print("OS error: {0}".format(err))
    except ValueError as err:
        print "Failed to read primitives from file!"
        print("Value error: {0}".format(err))
    except TypeError as err:
        print "Failed to read primitives from file!"
        print("Type error: {0}".format(err))
    except KeyError as err:
        print "Failed to read primitives from file!"
        print("Key error: {0}".format(err))
    except NameError as err:
        print "Failed to read primitives from file!"
        print("Name error: {0}".format(err))
    except AssertionError as err:
        print "Failed to read primitives from file!"
        print("Assertion error: {0}".format(err))
    except:
        print "Failed to read primitives from file!"
        print("Unexpected error:", sys.exc_info()[0])

    return False
