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

import sys, math #, types, time, random, os
import numpy as np

# Our custom utility methods
from mprim_generator_utils import *
from mprim_generator_generators import *
from mprim_generator_primitives import *
from mprim_generator_io import *

# Use the original SBPL approach to defining motion primitives
def define_motion_primitives_extended(args):
    ''' Define the motion primitives for all quadrants based on generators '''
    # Arguments
    resolution              = args['resolution']
    numberofangles          = args['numberofanglesperquadrant']*4
    numberofprimsperangle   = args['numberofprimsperangle']

    # Cost multipliers (multiplier is used as costmult*cost)
    forwardcostmult         = args['forwardcostmult']
    forwardandturncostmult  = args['forwardandturncostmult']
    backwardcostmult        = args['backwardcostmult']
    turninplacecostmult     = args['turninplacecostmult']
    sidestepcostmult        = args['sidestepcostmult']

    # Define the primitive generators for the first quadrant
    # Using either default SBPL parameters or from input file
    primitive_generator_definitions = get_primitive_generators(args)

    # Generate the motion primitives for all quadrants given the generators
    primitive_definitions = generate_motion_primitives(primitive_generator_definitions, args)

    if (args['output'] != ""):
        write_motion_primitives_file(primitive_definitions,args)

    if (args['showPrimitives']):
        try:
            visualize_motion_primitives(primitive_definitions,args)
        except OSError as err:
            print "Failed to visualize primitives !"
            print("OS error: {0}".format(err))
        except NameError as err:
            print "Failed to visualize primitives!"
            print("Name error: {0}".format(err))
        except:
            print "Cannot visualize the primitives"
            print("     Unexpected error:", sys.exc_info()[0])


    print "Finished motion primitive generation!"
    return

# Parse the command line arguments
def read_command( argv ):
    """
    Processes the command used to run mprim_generator from the command line.
    """
    from optparse import OptionParser

    usageStr = """
    USAGE:      python mprim_generator.py <options>
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
    parser.add_option('-r', '--resolution', type='float', dest='resolution',
                      help=default('Translational resolution of the lattice grid'), default=0.025)
    parser.add_option('-m', '--minTurnRadius', type='float', dest='minTurnRadius',
                      help=default('Minimum forward turning radius (same units as resolution)'), default=0.05)
    parser.add_option('-n', '--numberOfAnglesPerQuadrant', type='int',
                      dest='numberOfAnglesPerQuadrant',
                      help=default('Number of angles per quadrant the lattice grid (preferably a power of 2)'),
                      metavar='ANGLES',  default=4)
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
                      help=default('Number of samples along motion primitive'), default=20)
    parser.add_option('-a', '--checkArcLineArcGenerator', action='store_true',
                   dest='check_arc_line_arc_generator',
                   help='Check arc-line-arc shortest path generator approach', default=True)
    parser.add_option('-w', '--wheelbase', type='float', dest='wheelbase',
                      help=default('Distance between wheels (wheelbase) in meters'), default=0.230)
    parser.add_option('-c', '--cost_factor', type='float', dest='cost_conversion_factor',
                      help=default('Travel distance to cost'), default=1000)

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
    args['resolution']                  = options.resolution
    args['numberofanglesperquadrant']   = options.numberOfAnglesPerQuadrant
    args['numberofprimsperangle']       = options.primPerAngle
    args['minTurnRadius']               = options.minTurnRadius
    args['numberOfSamples']             = options.numSamples
    args['numberofprimitivegenerators'] = 4 # Place holder (deterimined by loading data)

    # Cost multipliers (multiplier is used as costmult*cost)
    args['forwardcostmult']         = options.fwdCost
    args['forwardandturncostmult']  = options.fwdTurnCost
    args['backwardcostmult']        = options.backCost
    args['turninplacecostmult']     = options.turnInPlaceCost
    args['sidestepcostmult']        = options.sideStepCost
    args['cost_conversion_factor']  = options.cost_conversion_factor
    args['wheelbase']               = options.wheelbase
    args['check_arc_line_arc_generator'] = options.check_arc_line_arc_generator
    
    print "Cost conversion factor = ",args['cost_conversion_factor']


    # Display options
    args['showPrimitives']    = options.showPrimitives

    if (options.showPrimitives and not can_visualize_plt()):
        raise Exception("Cannot show primitives because MatPlotLib is not available!")

    if (options.numberOfAnglesPerQuadrant < 2):
        raise Exception(str("Construction requires that number of angles is a multiple of 8 (Preferably a power of 2)  numAngles=%d" % (options.numberOfAngles)))

    return args

if __name__ == '__main__':
    """
    The main function called when mprim_generator.py is run from the command line:

    > python mprim_generator.py

    See the usage string for more details.

    > python mprim_generator.py --help
    """
    args = read_command( sys.argv[1:] ) # Get options from command line
    print "Running mprim_generator with the following arguments:"
    print args

    print "Generate the motion primitives given definitions ..."
    define_motion_primitives_extended(args)
