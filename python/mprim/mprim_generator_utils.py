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
# if available import pylab (from matlibplot)
visualize_plt = False
try:
    import matplotlib.pylab as plt
    visualize_plt = True
    print "Can visualize plots if needed !"

except ImportError as e:
    print e
    pass

def can_visualize_plt():
    return visualize_plt

#
def matrix_size(mat, elem=None):
    """ get the size of numpy matrix along given dimension """
    if elem is None:
        return mat.shape
    else:
        return mat.shape[elem]

def wrap_angle(rotation_angle):
    """keep angles in +/- pi bounds"""
    if (rotation_angle > np.pi):
        rotation_angle = rotation_angle - 2.0*np.pi
    elif (rotation_angle < -np.pi):
        rotation_angle = rotation_angle + 2.0*np.pi
    return rotation_angle

def normalize_angle(angle):
    """Normalize an angle to the range [0, 2*pi)."""
    if abs(angle) > 2.0 * np.pi:
        angle = fmod(angle, 2.0 * np.pi)
    if angle < 0.0:
        angle += 2.0 * np.pi
    return angle

def shortest_angle_dist(a1, a2):
    """Compute the minor arc distance between two angles."""
    a1_norm = normalize_angle(a1)
    a2_norm = normalize_angle(a2)
    return min(abs(a1_norm - a2_norm), 2.0 * np.pi - abs(a2_norm - a1_norm))

def shortest_angle_diff(a1, a2):
    """Compute the minor arc difference between two angles."""
    a1_norm = normalize_angle(a1)
    a2_norm = normalize_angle(a2)
    dist = shortest_angle_dist(a1, a2)
    if shortest_angle_dist(a1_norm + dist, a2_norm) < shortest_angle_dist(a1_norm - dist, a2_norm):
        return -dist
    else:
        return dist


def get_normal_from_heading_vectors(v0,v1,vl):
    """
      Define normal vectors for turning circle calculations
      Given:
            v0 and v1 the robot heading vectors, and
            vl, the connecting vector between the start and end point positions
      Returns:
           vc0 & vc1 - the unit orthongonal vectors
           DT     - Define local coordinate frame
           DTO    - orthogonal to DT to define right hand coordinate frame
           DTx    - x coordinate along the DT vector
           DTy    - y coordinate in the DTO vector
           DTquad - which DT frame quadrant the terminal position is in
    """

    # vector orthogal to current headings
    vc0 = np.array((-v0[1],v0[0])) # assume positive rotation, then correct
    vc1 = np.array((-v1[1],v1[0]))

    # Define local coordinate frame for the normal calculations
    # DT - line of y=-x*Cot[dT/2] rotated by T1  (need reference to write up)
    DT = np.array((vc0[0]-vc1[0],vc0[1]-vc1[1] ))
    if (np.linalg.norm(DT) < 1.e-12):
        DT = 1.*v0

    DT = DT/np.linalg.norm(DT)
    if (DT[0] < 0.0):
        #print "reverse DT"
        DT = -DT

    # Define orthogonal for right handed coordinate frame
    DTO = np.array((-DT[1],DT[0]))

    # Final configuration in DT coordinates
    DTx = np.dot(vl,DT)
    DTy = np.dot(vl,DTO)
    DTquad = 0
    if (DTx >= 0.0):
        if (DTy >= 0.0):
            # First quadrant
            # No change to normals
            DTquad = 1
        else:
            # Fourth quadrant
            vc0 = -vc0
            vc1 = -vc1
            DTquad = 4
    else:
        if (DTy >= 0.0):
            # Second quadrant
            # No change to normals
            DTquad = 2
        else:
            # Third quadrant
            DTquad = 3
            vc0 = -vc0
            vc1 = -vc1

    # Normalize
    vc0 = vc0/np.linalg.norm(vc0)
    vc1 = vc1/np.linalg.norm(vc1)
    #print "DT (",DTx,", ",DTy,") ==",DTquad
    return (vc0,vc1,DT,DTO,DTx,DTy,DTquad)

def line_intersection(p0,v0,p1,v1):
    """
        Find the intersection parameter values (s,t) of two parameterized lines
            p0(s) = p0 + s*v0
            p1(t) = p1 + t*v1
        Returns (s,t) if intersection, otherwise (None, None)
    """
    denom = v0[1]*v1[0] - v0[0]*v1[1]
    if (abs(denom) > 0.0):
        # Find intersection of the two lines
        s = -(-(v1[1]*p0[0]) + v1[0]*p0[1] + v1[1]*p1[0] - v1[0]*p1[1])/denom
        t = -(-(v0[1]*p0[0]) + v0[0]*p0[1] + v0[1]*p1[0] - v0[0]*p1[1])/denom
        return (s,t)
    else:
        return (None, None)

def visualize_motion_primitives(primitive_definitions,args):
    """
        Visualize the primitives using MatPlotLib
    """
    if not visualize_plt:
        print "Cannot visualize motion primitives - no matplotlib!"
        return

    print args
    resolution              = args['resolution']
    numberofangles          = args['numberofanglesperquadrant']*4
    numberofprimsperangle   = args['numberofprimsperangle']
    angles                  = args['angles']

    # iterate over angles
    for angleind in range(numberofangles):
        print "angleid=",angleind
        currentangle = angles[angleind]

        plt.figure(angleind)  # use angleind to plot each primitive set on different window
        plt.hold(True)
        plt.grid(True)
        plt.axis('equal')
        plt.axis([-10*resolution,10*resolution,-10*resolution,10*resolution])
        plt.title(str("%d (%.3f)"%(angleind, currentangle*180./np.pi)))

        # major ticks every 5, minor ticks every 1 resolution
        major_ticks = np.arange(-10*resolution,10*resolution,5*resolution)
        minor_ticks = np.arange(-10*resolution,10*resolution,  resolution)
        ax = plt.gca()
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(minor_ticks, minor=True)
        ax.set_xticks(major_ticks)
        ax.set_yticks(major_ticks)
        ax.grid(which='both')
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)
        plt.axis('equal')

        # iterate over primitives
        for primind in range(0, numberofprimsperangle):
            prim_def_primitive = primitive_definitions[angleind][primind]
            #print prim_def_primitive
            endpose_c      = prim_def_primitive['endpose']
            intermcells_m  = prim_def_primitive['intermcells']
            actioncostmult = prim_def_primitive['actioncostmult']
            endpt          = prim_def_primitive['endpoint']
            print "primind=",primind, "  endpt=",endpt," "#,prim_def_primitive['style']

            hndl = plt.plot(intermcells_m[:,0], intermcells_m[:,1],linestyle="-")
            plt.plot([intermcells_m[0,0], intermcells_m[-1,0]], [intermcells_m[0, 1], intermcells_m[-1, 1]],color=hndl[0].get_color(),linestyle='None', marker="o")#
            plt.text(endpt[0]+0.075*resolution, endpt[1]+0.075*resolution*math.copysign(1.0, endpose_c[2]), str(primind)+" ("+str(endpose_c[2])+")")

        plt.waitforbuttonpress()  # uncomment to plot each primitive set one at a time

    #print "Hold windows open until a button is pressed"
    #plt.waitforbuttonpress()  # Uncomment to hold until button pressed
    print "Hold windows open until all closed"
    plt.show()  # Un comment to keep windows open
                #    until the program is terminated or all individually closed
    return

def visualize_plt_circle(plt, startpt, vc0, endpt,vc1,radius,circle,color_value):
    """
     Plot circles of given radius around to points (given vector of unit circle points)
        circle is a vector of (x,y) points along unit circle
    """
    pc0 = startpt[:2]+vc0*radius
    pc1 = endpt[:2]  +vc1*radius
    plt.plot(pc0[0], pc0[1],  marker="+",color=color_value)
    plt.plot(pc0[0]+radius*circle[0,:], pc0[1]+radius*circle[1,:],  color=color_value)
    plt.plot(pc1[0], pc1[1],  marker="+",color=color_value)
    plt.plot(pc1[0]+radius*circle[0,:], pc1[1]+radius*circle[1,:],  color=color_value)
    return

def visualize_primitive_geometry_figure(name, startpt,endpt,
                                        v1,vc0,vc1, ptp, ptp1, radius):
    """
    Plot geometry relative to start and end points based on
    turning radius and common tangents
    """
    if not visualize_plt:
        print "Cannot visualize primitive geometry - no matplotlib!"
        return
    #if (name != "valid"):
    #    print "skip drawing invalid tangent!"
    #    return

    fig = plt.figure()
    fig.hold(True)
    plt.grid(True)
    plt.axis('equal')
    plt.title(str("%s (%s)->(%s)"%(name, startpt,endpt)))
    #plt.plot([startpt[0],pi[0],endpt[0],pe[0]], [startpt[1], pi[1], endpt[1],pe[1]],color="cyan")
    plt.plot(startpt[0], startpt[1],linewidth="2",marker="o",color="green")
    plt.plot(endpt[0],   endpt[1],  linewidth="2",marker="o",color="red")

    # Tangent circles
    theta = np.arange(0,2.03*np.pi,np.pi/50.)
    circle = np.array(np.vstack((np.cos(theta),np.sin(theta))))
    visualize_plt_circle(plt, startpt,  vc0, endpt, vc1,radius, circle,"green")
    visualize_plt_circle(plt, startpt, -vc0, endpt,-vc1,radius, circle,"red")

    if (name == "valid"):
        plt.plot([ ptp[0], ptp1[0] ],[ ptp[1],   ptp1[1] ],  linewidth="2",marker="*",color="orange")
    else:
        plt.plot([ ptp[0], ptp1[0] ],[ ptp[1],   ptp1[1] ],  linewidth="2",marker="*",color="cyan")

    print "pt0=",ptp,"  pt1=",ptp1
    print "Wait for button press for geometry figure"
    plt.waitforbuttonpress()  # uncomment to plot each primitive set one at a time

    return

def get_max_radius_circles(startpt,vc0,endpt,vc1,min_turning_radius_m):
    """
        Find turning radius for CSC based on min turning radius and
        max circles that share a common tangent
    """

    # terms in quadratic equation                                     line wrap\
    ac = -4 + vc0[0]*vc0[0] + vc0[1]*vc0[1] - 2*vc0[0]*vc1[0]                  \
            + vc1[0]*vc1[0] - 2*vc0[1]*vc1[1] + vc1[1]*vc1[1]
    bc = 2*startpt[0]*vc0[0] - 2*endpt[0]*vc0[0] + 2*startpt[1]*vc0[1]         \
         - 2*endpt[1]*vc0[1] - 2*startpt[0]*vc1[0] + 2*endpt[0]*vc1[0]         \
         - 2*startpt[1]*vc1[1] + 2*endpt[1]*vc1[1]
    cc = startpt[0]*startpt[0] + startpt[1]*startpt[1] - 2*startpt[0]*endpt[0] \
         + endpt[0]*endpt[0] - 2*startpt[1]*endpt[1] + endpt[1]*endpt[1]
    disc = bc*bc - 4.0*ac*cc
    if (disc < 0.0):
        # Cannot solve for max radius
        vl = endpt[:2] - startpt[:2]
        shortest = np.linalg.norm(vl)
        #print "Warning: Cannot solve for max radius of Arc-Line-Arc "
        #print "     start =",startpt," v0=",v0," vc0=",vc0
        #print "     end   =",endpt," v1=",v1," vc1=",vc1
        #print "     dist  =",shortest, " vl =",vl
        #print "ac=",ac,"  bc=",bc," cc=",cc,"  disc=",disc
        #print " Using min turning radius!"
        return min_turning_radius_m

    radius = 1.*min_turning_radius_m # initialize as separate variable
    if (np.math.fabs(ac) > 1.e-6):
        #print "ac=",ac,"  bc=",bc," cc=",cc,"  disc=",disc
        r1 = (-bc + np.math.sqrt(disc))/(2.0*ac)
        r2 = (-bc - np.math.sqrt(disc))/(2.0*ac)
        rmin = np.min((r1,r2))
        rmax = np.max((r1,r2))
        if (rmax < min_turning_radius_m):
            vl = endpt[:2] - startpt[:2]
            shortest = np.linalg.norm(vl)
            #print "Warning: Arc-Line-Arc - max radius is less than min_turning_radius_m"
            #print "     start =",startpt," vc0=",vc0
            #print "     end   =",endpt  ," vc1=",vc1
            #print "     dist  =",shortest, " vl =",vl
            #print "     ac=",ac,"  bc=",bc," cc=",cc,"  disc=",disc
            #print "     r1=",r1,"  r2=",r2," min_turning_radius_m=",min_turning_radius_m
            #print "   Using min turning radius!"
            return min_turning_radius_m

        if (rmin > 0.0):
            #print "CSC: rmin =",rmin," rmax=",rmax
            #print "     start =",startpt," v0=",v0," vc0=",vc0
            #print "     end   =",endpt," v1=",v1," vc1=",vc1
            ##print "     dist  =",shortest, " vl =",vl
            #print "     ac=",ac,"  bc=",bc," cc=",cc,"  disc=",disc
            #print "     r1=",r1,"  r2=",r2
            rmax = rmin

        #print "     CSC: rmin =",rmin," rmax=",rmax
        return min_turning_radius_m*0.5 + rmax*0.5 # Chose a radius for CSC (@todo - make 0.5 a parameter)
    else:
        #print "   Warning: ac=0 ... Using min turning radius!"
        return min_turning_radius_m

def get_potential_internal_tangents(startpt,  vc0, endpt,  v1, vc1,radius, visualize=False):
    """
    Get potential tangents for CSC motion between start and end points
    consider two circles adjacent to each terminal point and all common
    internal tangents between them.
    """

    # Center points of the arcs
    pc0 = startpt[:2] + radius*vc0
    pc1 = endpt[:2]   + radius*vc1

    # Internal tangent line intersects the midpoint between the circles (for same radius)
    pi = 0.5*(pc0 + pc1)

    #print "internal tangent: radius=",radius," pc0=",pc0," pc1=",pc1," pi=",pi

    # Caclulate the angle between radius and vector to pi (midpoint)
    #      (hypotenus of rt. triangle)
    vi = pi - pc0
    di = np.linalg.norm(vi)

    if (di < radius):
        # Coincident circles - no internal tangent!
        #print " coincident circles - no internal tangent!"
        sys.stdout.flush()
        return (False,None,None,None,None,None,None,None,None,None,None)

    #print "vi=",vi,"  di=",di
    ai = np.math.acos(radius/di)
    vn = vi/np.linalg.norm(vi)
    #print "     vn=",vn,"  |vn|=",np.linalg.norm(vn)," ai=",ai

    rot = np.array(np.vstack((np.hstack(( np.cos(ai), -np.sin(ai))),
                             np.hstack(( np.sin(ai),  np.cos(ai))))))

    # Vector toward tangent points
    #print "rot=",rot
    vtp1 = np.dot(rot,vn)
    vtp2 = np.dot(vn,rot)  # opposite rotation treated as transpose

    vs0 = -vc0 # vector from center to start

    #print " vn=",vn
    #print "vs0=",vs0
    #print "vtp1=",vtp1
    #print "vtp2=",vtp2

    # Angles between tangent points and the starting point
    a1 = np.dot(vtp1,vs0)
    a2 = np.dot(vtp2,vs0)
    #print "a1=",a1
    #print "a2=",a2

    # Only one can be a potential tangent
    vtp = vtp1 # chose point 1 as likely tangent
    atp = np.math.acos(a1)
    if (a2 > a1):
        # Point 2 is potential tangent
        vtp = vtp2
        atp = np.math.acos(a2)

    # Assume atp is positive and check for negative rotation to correct
    rotv = np.cross(vs0,vtp)
    rad0 = radius
    if (rotv < 0.0):
        atp  = -atp
        rad0 = -radius

    # Point on the first circle
    ptp = pc0 + vtp*radius
    #print "     ** ptp0=",ptp

    dS = ptp - startpt[:2]
    start_is_valid = dS[0] > 0.0  # moving forward

    # vector to the intersection point
    vtv = pi - ptp
    # Tangent on the second circle
    ptp1 = ptp + 2.0*vtv

    if (np.linalg.norm(vtv) < 1.e-5):
        rot = np.cross(-vc0,vtp)
        vtv2 = np.cross(np.array((0,0,rot)), np.array((vtp[0],vtp[1],0)))[:2]
        if (np.linalg.norm(vtv2) < 1.e-6):
            #print " coincident circles - no internal tangent!"
            sys.stdout.flush()
            return (False,None,None,None,None,None,None,None,None,None,None)
        else:

            print "\n\n----------------------------"
            print "Invalid vtv=",vtv
            print " internal tangent: radius=",radius
            print "    pc0=",pc0," pc1=",pc1
            print "     pi=",pi, " ptp=",ptp
            print "startpt=",startpt
            print "  endpt=",endpt

            print "   -vc0=",vc0," X vtp=",vtp
            print " Updated vtv=",vtv2, " rot=",rot
            print "----------------------------"
            print "\n\n\n\n\n\n"
            sys.stdout.flush()
            vtv = vtv2



    #print "     ** ptp1=",ptp1

    vs1 = -vc1 # Center to end point
    vtc = ptp1-pc1
    vtcn = vtc/np.linalg.norm(vtc)
    # Determine angle between tangent point and the endpt
    ai1 = np.math.acos(np.dot(vs1,vtcn))


    # Determine if this circle has a positive or negative rotation about the Center
    rotv = np.cross(vs1,v1)
    #print " a X b = ",rotv

    rott = np.array(np.vstack((np.hstack(( np.cos(ai1), -np.sin(ai1))),
                               np.hstack(( np.sin(ai1),  np.cos(ai1))))))
    vt1 = np.dot(rott,v1)
    atp1 = -ai1
    rad1 = -radius
    if (rotv > 0.0):
        # positive rotation about center, so rotate heading by negative (transpose)
        vt1 = np.dot(v1,rott)
        atp1 = ai1
        rad1 = radius

    vtn = vtv/np.linalg.norm(vtv)
    test = np.dot(vt1,vtn)

    sys.stdout.flush()
    if (start_is_valid and np.fabs(test-1.0) < 1e-6):
        if (visualize):
            visualize_primitive_geometry_figure("valid",startpt,endpt,v1,vc0,vc1, ptp, ptp1, radius)
        return (True,vc0,ptp,atp,rad0,vtn,pi,vc1,ptp1,atp1,rad1)
    else:
        #print " invalid tangent = ",test
        if (visualize):
            visualize_primitive_geometry_figure("invalid",startpt,endpt,v1,vc0,vc1, ptp, ptp1, radius)
        return (False,vc0,ptp,atp,rad0,vtn,pi,vc1,ptp1,atp1,rad1)

def get_potential_external_tangents(startpt,  vc0, endpt,  v1, vc1,radius, visualize=False):
    """
    Get potential external tangents for CSC motion between start and end points
    Consider two circles adjacent to each terminal point and all common
    external tangents between them.
    """

    # Center points of the arcs
    pc0 = startpt[:2] + radius*vc0
    pc1 = endpt[:2]   + radius*vc1
    #print "radius=",radius," pc0=",pc0," pc1=",pc1

    # vector connecting the center points of two circles
    vc = pc1-pc0

    if (np.linalg.norm(vc) < 1.e-6):
        # Coincident circles - no specific tangent point
        # No linear segment required
        # Calculate the midpoint of transition and return data
        print "Coincident circles - "
        atp = (endpt[2] - startpt[2])*0.5
        atp1 = atp
        if (atp < 0.0):
            rad0 = -radius
            rad1 = -radius
        else:
            rad0 = radius
            rad1 = radius

        rott = np.array(np.vstack((np.hstack(( np.cos(atp), -np.sin(atp))),
                                   np.hstack(( np.sin(atp),  np.cos(atp))))))

        vt1 = -vc1 # vector from center to end point

        vt1 = np.dot(-vc1,rott) # post multiply for negative rotation as row vector

        ptp  = pc1 + vt1
        ptp1 = pc1 + vt1
        pi   = pc1
        vcn = vc  # This is (0,0)
        sys.stdout.flush()
        return (True,vc0,ptp,atp,rad0,vcn,pi,vc1,ptp1,atp1,rad1)

    # Regular calculation

    # Normal vector along line connecting circle centers
    vcn = vc/np.linalg.norm(vc)

    # Midpoint along that line
    pi = pc0 + 0.5*vc

    # Vector toward tangent points
    rot = np.array(np.vstack((np.hstack((  0.0, -1.0)),
                              np.hstack((  1.0,  0.0)))))
    vtp1 = np.dot(rot,vcn)
    vtp2 = np.dot(vcn,rot)  # opposite rotation treated as transpose

    vs0 = -vc0 # vector from center to start

    # Angles between tangent points and the starting point
    a1 = np.dot(vtp1,vs0)
    a2 = np.dot(vtp2,vs0)
    #print "a1=",a1
    #print "a2=",a2

    # Only one can be a potential tangent
    vtp = vtp1 # chose point 1 as likely tangent
    atp = np.math.acos(a1)
    if (a2 > a1):
        # Point 2 is potential tangent
        vtp = vtp2
        atp = np.math.acos(a2)

    # Assume atp is positive and check for negative rotation to correct
    rotv = np.cross(vs0,vtp)
    rad0 = radius
    if (rotv < 0.0):
        atp  = -atp
        rad0 = -radius


    # Point on the first circle
    ptp = pc0 + vtp*radius
    #print "     ** ptp0=",ptp

    dS = ptp - startpt[:2]
    start_is_valid = dS[0] > 0.0  # moving forward

    # Tangent on the second circle
    ptp1 = ptp + vc
    #print "     ** ptp1=",ptp1

    vs1 = -vc1 # Center to end point

    vtc = ptp1-pc1
    vtcn = vtc/np.linalg.norm(vtc)
    # Determine angle between tangent point and the endpt
    ai1 = np.math.acos(np.dot(vs1,vtcn))
    #print " vc1=",vc1,"  vs1=",vs1
    #print "   dP=",ptp1-pc1

    # Determine if this circle has a positive or negative rotation about the Center
    rotv = np.cross(vs1,v1)
    #print " a X b = ",rotv

    rott = np.array(np.vstack((np.hstack(( np.cos(ai1), -np.sin(ai1))),
                               np.hstack(( np.sin(ai1),  np.cos(ai1))))))
    vt1 = np.dot(rott,v1)
    atp1 = -ai1 # Rotation angle of curve after tangent point
    rad1 = -radius
    if (rotv > 0.0):
        # positive rotation about center, so rotate heading by negative (transpose)
        vt1 = np.dot(v1,rott)
        atp1 = ai1
        rad1 = radius

    test  = np.dot(vt1 ,vcn)
    test2 = np.dot(vtcn,vt1)

    sys.stdout.flush()
    if (start_is_valid and np.fabs(test-1.0) < 1e-6):

        if (visualize):
            visualize_primitive_geometry_figure("valid",startpt,endpt,v1,vc0,vc1, ptp, ptp1, radius)

        return (True,vc0,ptp,atp,rad0,vcn,pi,vc1,ptp1,atp1,rad1)
    else:
        #print " invalid tangent = ",test,"  test2=",test2
        #print "     vt1=",vt1," angle=",atp1
        #print "      v1=",v1
        #print "     vcn=",vcn
        #print "     rotv=",rotv
        if (visualize):
            visualize_primitive_geometry_figure("invalid",startpt,endpt,v1,vc0,vc1, ptp, ptp1,radius)

        return (False,vc0,ptp,atp, rad0, vcn,pi,vc1,ptp1,atp1, rad1)


def default(str):
    """
       Helper for printing default arguments for parser help message
    """
    return str + ' [Default: %default]'
