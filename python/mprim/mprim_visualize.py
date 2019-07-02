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

import math, sys #, types, time, random
import numpy as np

# Our custom utility methods
from mprim_generator_utils import *
from mprim_generator_generators import *
from mprim_generator_primitives import *
from mprim_generator_io import *

# Parse the command line arguments
def read_command( argv ):
    """
    Processes the command used to run mprim_generator from the command line.
    """
    from optparse import OptionParser

    usageStr = """
    USAGE:      python mprim_visualize.py <options>
                     - Loads and displayes the defined motion primitives

    """
    parser = OptionParser(usageStr)

    parser.add_option('-o', '--output', dest='output', type='string',
                      help=default('The output file name (including path information)'),
                      metavar='FILE',
                      default="")

    parser.add_option('-i', '--input', dest='filename', type='string',
                      help=default('Optional input file name with primitive terminal grid points'),
                      metavar='PRIM_DEFINITION_FILE',
                      default="")

    parser.add_option('-s', '--showPrimitives', action='store_true',
                      dest='showPrimitives',
                      help='Show the generated primitives',
                      default=True)
    options, otherjunk = parser.parse_args(argv)

    # Create a dictionary of arguments
    args = dict()
    args['output'] = options.output
    args['filename'] = options.filename
    args['showPrimitives'] = options.showPrimitives

    if len(otherjunk) == 1 and args['filename'] == "":
        args['filename'] = otherjunk[0]
    elif args['filename'] == "":
        print args
        print otherjunk
        raise Exception('Command line input not understood: ' + str(otherjunk))

    if args['filename'] == None or args['filename'] == "":
        raise Exception("The input file is not defined !")

    if (options.showPrimitives and not can_visualize_plt()):
        raise Exception("Cannot show primitives because MatPlotLib is not available!")

    return args

if __name__ == '__main__':
    """
    The main function called when mprim_visualize.py is run from the command line:

    > python mprim_visualize.py <filename>

    See the usage string for more details.

    > python mprim_visualize.py --help
    """
    args = read_command( sys.argv[1:] ) # Get options from command line
    print "Running mprim_generator with the following arguments:"
    print args

    print "Generate the motion primitives given definitions ..."
    primitive_definitions = []
    if not read_motion_primitives_file(args['filename'],
                                       primitive_definitions,
                                       args):
        print "Failed to read the primitive file <",args['filename'],'>!'
        sys.exit(-1)

    if args['output'] != "":
        print "Re-write output for testing to <",args['output'],"> ..."
        write_motion_primitives_file(primitive_definitions, args)
        print "Done re-write!"


    print args
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
            print "     Unexpected error:", sys.exc_info()[0]

    print "Finished!"
