
from rawe.collocation import trajectory
import numpy

def interpolate(traj,pps,times,names):
    ############# interpolate ###########
    sys.stdout.write('reticulating splines... ')

    h = (traj.tgrid[-1,0,0] - traj.tgrid[0,0,0])/float(traj.dvMap._nk*traj.dvMap._nicp)
    h *= traj.dvMap._nk*traj.dvMap._nicp/float(self.nk*self.nicp)
    h *= numLoops

    t0 = 0.0
    x_interp_times = []
    for timestepIdx in range(self.nk):
        for nicpIdx in range(self.nicp):
            for degIdx in range(self.deg+1):
                time = t0 + h*self.lagrangePoly.tau_root[degIdx]
                if time > traj.tgrid[-1,0,0]:
                    time -= traj.tgrid[-1,0,0]
                x_interp_times.append(time)
            t0 += h
            if t0 > traj.tgrid[-1,0,0]:
                t0 -= traj.tgrid[-1,0,0]
    x_interp_times.append(t0)

if __name__=='__main__':
    traj = trajectory.load_traj('data/transition_backup.dat')
    print 'making piecewise polys ...'
    pps = trajectory.make_pps(traj)
    t0 = traj.tgrid[0,0,0]
    tF = traj.tgrid[-1,0,0]
    times = numpy.linspace(t0,tF,tF*24)
    for name in ['x','y','z','aileron']:
        vals = pps[name](times)
        print '# '+name
        print ' '.join([('[%d]%.15f' % (k,val)) for k,val in enumerate(vals)])
    sin_delta = pps['sin_delta'](times)
    cos_delta = pps['cos_delta'](times)
    vals = [numpy.degrees(numpy.arctan2(s,c)) for (s,c) in zip(sin_delta,cos_delta)]
    print '# delta'
    print ' '.join([('[%d]%.15f' % (k,val)) for k,val in enumerate(vals)])
