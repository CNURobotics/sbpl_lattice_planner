#% /*
#%  * Copyright (c) 2016, David Conner (Christopher Newport University)
#%  * Based on genmprim.m 
#%  * Copyright (c) 2008, Maxim Likhachev
#%  * All rights reserved.
#%  * converted by libermat utility (https://github.com/awesomebytes/libermate)
#%  *
#%  * Redistribution and use in source and binary forms, with or without
#%  * modification, are permitted provided that the following conditions are met:
#%  * 
#%  *     * Redistributions of source code must retain the above copyright
#%  *       notice, this list of conditions and the following disclaimer.
#%  *     * Redistributions in binary form must reproduce the above copyright
#%  *       notice, this list of conditions and the following disclaimer in the
#%  *       documentation and/or other materials provided with the distribution.
#%  *     * Neither the name of the Carnegie Mellon University nor the names of its
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

import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def genmprim(outfilename):

    #%
    #%generates motion primitives and saves them into file
    #%
    #%written by Maxim Likhachev
    #%---------------------------------------------------
    #%
    #%defines
    LINESEGMENT_MPRIMS = 1.
    #%set the desired type of motion primitives
    UNICYCLE_MPRIMS = 0.
    if LINESEGMENT_MPRIMS == 1.:
        resolution = 0.01
        numberofangles = 32
        #%preferably a power of 2, definitely multiple of 8
        numberofprimsperangle = 16
        #%multipliers (multiplier is used as costmult*cost)
        forwardcostmult = 1.
        backwardcostmult = 5.
        forwardandturncostmult = 1.
        sidestepcostmult = 50.
        turninplacecostmult = 50.
        #%note, what is shown x,y,theta changes (not absolute numbers)
        #%0 degreees
        basemprimendpts0_c = np.zeros((numberofprimsperangle, 4))
        #%x,y,theta,costmult 
        #%x aligned with the heading of the robot, angles are positive
        #%counterclockwise
        #%0 theta change
        basemprimendpts0_c[0,:] = np.array(np.hstack((1., 0., 0., forwardcostmult)))
        basemprimendpts0_c[1,:] = np.array(np.hstack((4., 0., 0., forwardcostmult)))
        basemprimendpts0_c[2,:] = np.array(np.hstack((8., 0., 0., forwardcostmult)))
        basemprimendpts0_c[3,:] = np.array(np.hstack((6., 2., 0., sidestepcostmult)))
        basemprimendpts0_c[4,:] = np.array(np.hstack((6., -2., 0., sidestepcostmult)))
        basemprimendpts0_c[5,:] = np.array(np.hstack((2., 3., 0., sidestepcostmult)))
        basemprimendpts0_c[6,:] = np.array(np.hstack((2., -3., 0., sidestepcostmult)))
        basemprimendpts0_c[7,:] = np.array(np.hstack((-5., 0., 0., backwardcostmult)))
		
		#%1/32 theta change
        basemprimendpts0_c[8,:] = np.array(np.hstack((6., 2., 1., forwardandturncostmult)))
        basemprimendpts0_c[9,:] = np.array(np.hstack((6., -2., -1., forwardandturncostmult)))
        #%2/32 theta change
        basemprimendpts0_c[10,:] = np.array(np.hstack((4., 3., 2., forwardandturncostmult)))
        basemprimendpts0_c[11,:] = np.array(np.hstack((4., -3., -2., forwardandturncostmult)))
        #%turn in place
        basemprimendpts0_c[12,:] = np.array(np.hstack((0., 0., 1., turninplacecostmult)))
        basemprimendpts0_c[13,:] = np.array(np.hstack((0., 0., -1., turninplacecostmult)))
        basemprimendpts0_c[14,:] = np.array(np.hstack((0., 0., 3., turninplacecostmult)))
        basemprimendpts0_c[15,:] = np.array(np.hstack((0., 0., -3., turninplacecostmult)))
        print basemprimendpts0_c
        #%45 degrees
        basemprimendpts45_c = np.zeros((numberofprimsperangle, 4)) #%x,y,theta,costmult (multiplier is used as costmult*cost)
        #%x aligned with the heading of the robot, angles are positive
        #%counterclockwise
        #%0 theta change 
        basemprimendpts45_c[0,:] = np.array(np.hstack((1., 1., 0., forwardcostmult)))
        basemprimendpts45_c[1,:] = np.array(np.hstack((3., 3., 0., forwardcostmult)))
        basemprimendpts45_c[2,:] = np.array(np.hstack((6., 6., 0., forwardcostmult)))
        basemprimendpts45_c[3,:] = np.array(np.hstack((2., 6., 0., sidestepcostmult)))
        basemprimendpts45_c[4,:] = np.array(np.hstack((6., 2., 0., sidestepcostmult)))
        basemprimendpts45_c[5,:] = np.array(np.hstack((0., 4., 0., sidestepcostmult)))
        basemprimendpts45_c[6,:] = np.array(np.hstack((4., 0., 0., sidestepcostmult)))
        basemprimendpts45_c[7,:] = np.array(np.hstack((-4., -4., 0., backwardcostmult)))
        #%1/32 theta change
        basemprimendpts45_c[8,:] = np.array(np.hstack((2., 6., 1., forwardandturncostmult)))
        basemprimendpts45_c[9,:] = np.array(np.hstack((6., 2., -1., forwardandturncostmult)))
        #%2/32 theta change
        basemprimendpts45_c[10,:] = np.array(np.hstack((1., 5., 2., forwardandturncostmult)))
        basemprimendpts45_c[11,:] = np.array(np.hstack((5., 1., -2., forwardandturncostmult)))
        #%turn in place
        basemprimendpts45_c[12,:] = np.array(np.hstack((0., 0., 1., turninplacecostmult)))
        basemprimendpts45_c[13,:] = np.array(np.hstack((0., 0., -1., turninplacecostmult)))
        basemprimendpts45_c[14,:] = np.array(np.hstack((0., 0., 3., turninplacecostmult)))
        basemprimendpts45_c[15,:] = np.array(np.hstack((0., 0., -3., turninplacecostmult)))
        print basemprimendpts45_c
		
        #%22.5 degrees
        basemprimendpts22p5_c = np.zeros((numberofprimsperangle, 4))
        #%x,y,theta,costmult (multiplier is used as costmult*cost)
        #%x aligned with the heading of the robot, angles are positive
        #%counterclockwise
        #%0 theta change     
        basemprimendpts22p5_c[0,:] = np.array(np.hstack((2., 1., 0., forwardcostmult)))
        basemprimendpts22p5_c[1,:] = np.array(np.hstack((4., 2., 0., forwardcostmult)))
        basemprimendpts22p5_c[2,:] = np.array(np.hstack((6., 3., 0., forwardcostmult)))
        basemprimendpts22p5_c[3,:] = np.array(np.hstack((4., 4., 0., sidestepcostmult)))
        basemprimendpts22p5_c[4,:] = np.array(np.hstack((6., 2., 0., sidestepcostmult)))
        basemprimendpts22p5_c[5,:] = np.array(np.hstack((0., 3., 0., sidestepcostmult)))
        basemprimendpts22p5_c[6,:] = np.array(np.hstack((4., -1., 0., sidestepcostmult)))
        basemprimendpts22p5_c[7,:] = np.array(np.hstack((-4., -2., 0., backwardcostmult)))
        #%1/32 theta change
        basemprimendpts22p5_c[8,:] = np.array(np.hstack((4., 4., 1., forwardandturncostmult)))
        basemprimendpts22p5_c[9,:] = np.array(np.hstack((6., 2., -1., forwardandturncostmult)))
        #%2/32 theta change
        basemprimendpts22p5_c[10,:] = np.array(np.hstack((2., 4., 2., forwardandturncostmult)))
        basemprimendpts22p5_c[11,:] = np.array(np.hstack((6., 0., -2., forwardandturncostmult)))
        #%turn in place
        basemprimendpts22p5_c[12,:] = np.array(np.hstack((0., 0., 1., turninplacecostmult)))
        basemprimendpts22p5_c[13,:] = np.array(np.hstack((0., 0., -1., turninplacecostmult)))
        basemprimendpts22p5_c[14,:] = np.array(np.hstack((0., 0., 3., turninplacecostmult)))
        basemprimendpts22p5_c[15,:] = np.array(np.hstack((0., 0., -3., turninplacecostmult)))
        #%11.25 degrees
        basemprimendpts11p25_c = np.zeros((numberofprimsperangle, 4))
        #%x,y,theta,costmult (multiplier is used as costmult*cost)
        #%x aligned with the heading of the robot, angles are positive
        #%counterclockwise
        #%0 theta change     
        basemprimendpts11p25_c[0,:] = np.array(np.hstack((3., 1., 0., forwardcostmult)))
        basemprimendpts11p25_c[1,:] = np.array(np.hstack((6., 2., 0., forwardcostmult)))
        basemprimendpts11p25_c[2,:] = np.array(np.hstack((9., 3., 0., forwardcostmult)))
        basemprimendpts11p25_c[3,:] = np.array(np.hstack((4., 3., 0., sidestepcostmult)))
        basemprimendpts11p25_c[4,:] = np.array(np.hstack((6., 0., 0., sidestepcostmult)))
        basemprimendpts11p25_c[5,:] = np.array(np.hstack((1., 3., 0., sidestepcostmult)))
        basemprimendpts11p25_c[6,:] = np.array(np.hstack((3., -2., 0., sidestepcostmult)))
        basemprimendpts11p25_c[7,:] = np.array(np.hstack((-6., -2., 0., backwardcostmult)))
        #%1/32 theta change
        basemprimendpts11p25_c[8,:] = np.array(np.hstack((4., 3., 1., forwardandturncostmult)))
        basemprimendpts11p25_c[9,:] = np.array(np.hstack((6., 0., -1., forwardandturncostmult)))
        #%2/32 theta change
        basemprimendpts11p25_c[10,:] = np.array(np.hstack((2., 4., 2., forwardandturncostmult)))
        basemprimendpts11p25_c[11,:] = np.array(np.hstack((5., -1., -2., forwardandturncostmult)))
        #%turn in place
        basemprimendpts11p25_c[12,:] = np.array(np.hstack((0., 0., 1., turninplacecostmult)))
        basemprimendpts11p25_c[13,:] = np.array(np.hstack((0., 0., -1., turninplacecostmult)))
        basemprimendpts11p25_c[14,:] = np.array(np.hstack((0., 0., 3., turninplacecostmult)))
        basemprimendpts11p25_c[15,:] = np.array(np.hstack((0., 0., -3., turninplacecostmult)))
        #%33.75 degrees
        basemprimendpts33p75_c = np.zeros((numberofprimsperangle, 4))
        #%x,y,theta,costmult 
        #%x aligned with the heading of the robot, angles are positive
        #%counterclockwise
        #%0 theta change
        basemprimendpts33p75_c[0,:] = np.array(np.hstack((3., 2., 0., forwardcostmult)))
        basemprimendpts33p75_c[1,:] = np.array(np.hstack((6., 4., 0., forwardcostmult)))
        basemprimendpts33p75_c[2,:] = np.array(np.hstack((9., 6., 0., forwardcostmult)))
        basemprimendpts33p75_c[3,:] = np.array(np.hstack((4., 5., 0., sidestepcostmult)))
        basemprimendpts33p75_c[4,:] = np.array(np.hstack((6., 2., 0., sidestepcostmult)))
        basemprimendpts33p75_c[5,:] = np.array(np.hstack((0., 4., 0., sidestepcostmult)))
        basemprimendpts33p75_c[6,:] = np.array(np.hstack((3., -2., 0., sidestepcostmult)))
        basemprimendpts33p75_c[7,:] = np.array(np.hstack((-6., -4., 0., backwardcostmult)))
        #%1/32 theta change
        basemprimendpts33p75_c[8,:] = np.array(np.hstack((4., 5., 1., forwardandturncostmult)))
        basemprimendpts33p75_c[9,:] = np.array(np.hstack((6., 2., -1., forwardandturncostmult)))
        #%2/32 theta change    
        basemprimendpts33p75_c[10,:] = np.array(np.hstack((1., 5., 2., forwardandturncostmult)))
        basemprimendpts33p75_c[11,:] = np.array(np.hstack((3., -2., -2., forwardandturncostmult)))
        #%turn in place
        basemprimendpts33p75_c[12,:] = np.array(np.hstack((0., 0., 1., turninplacecostmult)))
        basemprimendpts33p75_c[13,:] = np.array(np.hstack((0., 0., -1., turninplacecostmult)))
        basemprimendpts33p75_c[14,:] = np.array(np.hstack((0., 0., 3., turninplacecostmult)))
        basemprimendpts33p75_c[15,:] = np.array(np.hstack((0., 0., -3., turninplacecostmult)))
        print basemprimendpts33p75_c
			
    elif UNICYCLE_MPRIMS == 1.:
        print( 'ERROR: unsupported mprims type\n')
        return []
        
    else:
        print( 'ERROR: undefined mprims type\n')
        return []
        
    
    fout = open(outfilename, 'w')
    #%write the header
    fout.write( 'resolution_m: %f\n' %( resolution))
    fout.write( 'numberofangles: %d\n'%( numberofangles))
    fout.write( 'totalnumberofprimitives: %d\n'%( numberofprimsperangle*numberofangles))
    #%iterate over angles
    for angleind in np.arange(1, (numberofangles)+1):
        #%current angle
        currentangle = (angleind-1.)*2.*np.pi/numberofangles
        currentangle_36000int = int(np.round((angleind-1)*36000./numberofangles))

        plt.figure(angleind)
        plt.hold(True)
        plt.grid(True)
        plt.axis('equal')
        plt.axis([-7*resolution,7*resolution,-7*resolution,7*resolution])
        plt.title(str(angleind)+" "+str(currentangle_36000int/100.))
        #%iterate over primitives    
        for primind in np.arange(1., (numberofprimsperangle)+1):
            fout.write( 'primID: %d\n'% (primind-1))
            fout.write( 'startangle_c: %d\n'% (angleind-1))
			
            print "angleID=",angleind-1
            print "primID=",primind-1
            print "current_angle=",currentangle
            print "currentangle_36000int=",currentangle_36000int
			
            #%compute which template to use
            if (currentangle_36000int%9000) == 0:
                print "(currentangle_36000int%9000) == 0"
                basemprimendpts_c = basemprimendpts0_c[int(primind)-1,:]
                angle = currentangle
            elif (currentangle_36000int%4500) == 0:
                print "(currentangle_36000int%4500)"
                basemprimendpts_c = basemprimendpts45_c[int(primind)-1,:]
                angle = currentangle-45.*np.pi/180.
                
            elif ((currentangle_36000int-7875)%9000) == 0:
                print "((currentangle_36000int-7875.)%9000)"
                basemprimendpts_c = 1*basemprimendpts33p75_c[int(primind)-1,:] # 1* to force deep copy to avoid reference update below
                basemprimendpts_c[0] = basemprimendpts33p75_c[int(primind)-1,1]
                #%reverse x and y
                basemprimendpts_c[1] = basemprimendpts33p75_c[int(primind)-1,0]
                basemprimendpts_c[2] = -basemprimendpts33p75_c[int(primind)-1,2]
                #%reverse the angle as well
                angle = currentangle-(78.75*np.pi)/180.                
            elif ((currentangle_36000int-6750)%9000) == 0:
                print "((currentangle_36000int-6750.)%9000)"
                basemprimendpts_c = 1*basemprimendpts22p5_c[int(primind)-1,:] # 1* to force deep copy to avoid reference update below
                basemprimendpts_c[0] = basemprimendpts22p5_c[int(primind)-1,1]
                #%reverse x and y
                basemprimendpts_c[1] = basemprimendpts22p5_c[int(primind)-1,0]
                basemprimendpts_c[2] = -basemprimendpts22p5_c[int(primind)-1,2]
                #%reverse the angle as well
                #%fprintf(1, '%d %d %d onto %d %d %d\n', basemprimendpts22p5_c(1), basemprimendpts22p5_c(2), basemprimendpts22p5_c(3), ...
                #%    basemprimendpts_c(1), basemprimendpts_c(2), basemprimendpts_c(3));
                angle = currentangle-(67.5*np.pi)/180.
                
            elif ((currentangle_36000int-5625)%9000) == 0:
                print "((currentangle_36000int-5625.)%9000)"
                basemprimendpts_c = 1*basemprimendpts11p25_c[int(primind)-1,:] # 1* to force deep copy to avoid reference update below
                basemprimendpts_c[0] = basemprimendpts11p25_c[int(primind)-1,1]
                #%reverse x and y
                basemprimendpts_c[1] = basemprimendpts11p25_c[int(primind)-1,0]
                basemprimendpts_c[2] = -basemprimendpts11p25_c[int(primind)-1,2]
                #%reverse the angle as well
                angle = currentangle-(56.25*np.pi)/180.
                
            elif ((currentangle_36000int-3375)%9000) == 0:
                print "((currentangle_36000int-3375.)%9000)"
                basemprimendpts_c = basemprimendpts33p75_c[int(primind)-1,:]
                angle = currentangle-(33.75*np.pi)/180.
                
            elif ((currentangle_36000int-2250)%9000) == 0:
                print "((currentangle_36000int-2250.)%9000)"
                basemprimendpts_c = basemprimendpts22p5_c[int(primind)-1,:]
                angle = currentangle-(22.5*np.pi)/180.
                
            elif ((currentangle_36000int-1125)%9000) == 0:
                print "((currentangle_36000int-1125.)%9000)"
                basemprimendpts_c = basemprimendpts11p25_c[int(primind)-1,:]
                angle = currentangle-(11.25*np.pi)/180.
                
            else:
                print( 'ERROR: invalid angular resolution. angle = %d\n', currentangle_36000int)
                return []
                
            print "basemprimendpts_c=",basemprimendpts_c
            print "angle=",angle
			
			
            #%now figure out what action will be        
            baseendpose_c = basemprimendpts_c[0:3]
            additionalactioncostmult = basemprimendpts_c[3]
            endx_c = np.round((baseendpose_c[0]*np.cos(angle))-(baseendpose_c[1]*np.sin(angle)))
            endy_c = np.round((baseendpose_c[0]*np.sin(angle))+(baseendpose_c[1]*np.cos(angle)))
            endtheta_c = (angleind-1+baseendpose_c[2])%numberofangles
            endpose_c = np.array(np.hstack((endx_c, endy_c, endtheta_c)))
            print "endpose_c=",endpose_c            
            print( 'rotation angle=%f\n'% (angle*180./np.pi))
            #if (baseendpose_c[1] == 0.) and ( baseendpose_c[2] == 0.):
            #    print('endpose=%d %d %d\n' % ( endpose_c[0], endpose_c[1], endpose_c[2]))
            
            
            #%generate intermediate poses (remember they are w.r.t 0,0 (and not centers of the cells)
            numofsamples = 10
            intermcells_m = np.zeros((numofsamples, 3))
            if LINESEGMENT_MPRIMS == 1.:
                startpt = np.array(np.hstack((0., 0., currentangle)))
                endpt = np.array(np.hstack(((endpose_c[0]*resolution), 
                                            (endpose_c[1]*resolution), 
                                            (( ((angleind-1+baseendpose_c[2])%numberofangles)*2.*np.pi)/numberofangles))))
				
                print "startpt =",startpt
                print "endpt   =",endpt
				
                intermcells_m = np.zeros((numofsamples, 3))
                for iind in np.arange(1, (numofsamples)+1):
                    fraction = float(iind-1)/(numofsamples)					
                    intermcells_m[int(iind)-1,:] = np.array((startpt[0]+(endpt[0]-startpt[0])*fraction, startpt[1]+(endpt[1]-startpt[1])*fraction, 0))
                    rotation_angle = baseendpose_c[2]*(2.*np.pi/numberofangles)
                    intermcells_m[int(iind)-1,2] = (startpt[2]+rotation_angle*fraction)%(2.*np.pi)
                    #print " ",iind,"  of ",numofsamples," fraction=",fraction," rotation=",rotation_angle

            print "  intermcells=",intermcells_m
            #%write out
            fout.write( 'endpose_c: %d %d %d\n'% (endpose_c[0], endpose_c[1], endpose_c[2]))
            fout.write( 'additionalactioncostmult: %d\n'%(additionalactioncostmult))
            fout.write( 'intermediateposes: %d\n'%( matcompat.size(intermcells_m, 1.)))
            for interind in np.arange(1., (matcompat.size(intermcells_m, 1.))+1):
                fout.write( '%.4f %.4f %.4f\n' %( intermcells_m[int(interind)-1,0], intermcells_m[int(interind)-1,1], intermcells_m[int(interind)-1,2]))
                
            plt.plot(intermcells_m[:,0], intermcells_m[:,1],linestyle="-",marker="o")
            plt.text(endpt[0], endpt[1], str(endpose_c[2]))
            plt.hold(True)
        plt.waitforbuttonpress()

 
    fout.close()
    #plt.close()
    return []
