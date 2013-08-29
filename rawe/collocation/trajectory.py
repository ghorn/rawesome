# Copyright 2012-2013 Greg Horn
#
# This file is part of rawesome.
#
# rawesome is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rawesome is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with rawesome.  If not, see <http://www.gnu.org/licenses/>.

import matplotlib.pyplot as plt
import pickle
import numpy
import scipy.io
from scipy.interpolate import PiecewisePolynomial

import casadi as C
import collmaps


def load_traj(filename):
    '''
    take in a filename
    return an unpicked trajectory object
    '''
    print "loading trajectory from "+filename+" ..."
    f=open(filename,'r')
    traj = pickle.load(f)
    f.close()
    assert isinstance(traj,Trajectory), "the file \""+filename+"\" doesn't have a pickled Trajectory"
    return traj

def make_pps(traj):
    '''
    take a trajectory object as input
    return a dictionary of piecewise polynomials, one for each name
    '''
    pps = {}

    # differential states
    for name in traj.dvMap._xNames:
        # make piecewise poly
        pps[name] = None
        for timestepIdx in range(traj.dvMap._nk):
            for nicpIdx in range(traj.dvMap._nicp):
                ts = []
                ys = []
                for degIdx in range(traj.dvMap._deg+1):
                    ts.append(traj.tgrid[timestepIdx,nicpIdx,degIdx])
                    ys.append([traj.dvMap.lookup(name,timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx)])
                if pps[name] is None:
                    pps[name] = PiecewisePolynomial(ts,ys)
                else:
                    pps[name].extend(ts,ys)
        pps[name].extend([traj.tgrid[-1,0,0]],[[traj.dvMap.lookup(name,timestep=-1,nicpIdx=0,degIdx=0)]])

    # algebraic variables
    for name in traj.dvMap._zNames:
        # make piecewise poly
        pps[name] = None
        for timestepIdx in range(traj.dvMap._nk):
            for nicpIdx in range(traj.dvMap._nicp):
                ts = []
                ys = []
                for degIdx in range(1,traj.dvMap._deg+1):
                    ts.append(traj.tgrid[timestepIdx,nicpIdx,degIdx])
                    ys.append([traj.dvMap.lookup(name,timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx)])
                if pps[name] is None:
                    pps[name] = PiecewisePolynomial(ts,ys)
                else:
                    pps[name].extend(ts,ys)

    # controls
    for name in traj.dvMap._uNames:
        # make piecewise poly
        ts = []
        ys = []
        for timestepIdx in range(traj.dvMap._nk):
            ts.append(traj.tgrid[timestepIdx,0,0])
            ys.append([traj.dvMap.lookup(name,timestep=timestepIdx)])
        pps[name] = PiecewisePolynomial(ts,ys)

    return pps


class Trajectory(object):
    """
    Trajectory contains x/z/u/p collmap, and output map and quadrature map.
    It can lookup values from these things.
    It is the thing which is saved and loaded.
    """
    def __init__(self,ocp,v_opt):
        self.dvMap = collmaps.VectorizedReadOnlyCollMap(ocp,'devectorized design vars',v_opt)
        self.outputMap = collmaps.OutputMap(ocp._outputMapGenerator, v_opt)
        self.quadratureMap = collmaps.QuadratureMap(ocp._quadratureManager, v_opt)

        self.nk = ocp.nk
        self.nicp = ocp.nicp
        self.deg = ocp.deg
        self.collPoly = ocp.collPoly

        # make time grid
        ocp.hfun.setInput(v_opt)
        ocp.hfun.evaluate()
        h = float(ocp.hfun.output())

        self.tgrid = numpy.resize([],(ocp.nk+1,ocp.nicp,ocp.deg+1))
        tf = 0.0
        for k in range(ocp.nk):
            for i in range(ocp.nicp):
                self.tgrid[k,i,:] = tf + h*numpy.array(ocp.lagrangePoly.tau_root)
                tf += h
        self.tgrid[ocp.nk,0,0] = tf

    def getDvs(self):
        return self.dvMap.vectorize()

    def lookup(self,name,timestep=None,nicpIdx=None,degIdx=None):
        try:
            return self.dvMap.lookup(name,timestep=timestep,nicpIdx=nicpIdx,degIdx=degIdx)
        except NameError:
            pass
        try:
            return self.outputMap.lookup(name,timestep=timestep,nicpIdx=nicpIdx,degIdx=degIdx)
        except NameError:
            pass
        try:
            return self.quadratureMap.lookup(name,timestep,nicpIdx,degIdx)
        except NameError:
            pass
        raise NameError("lookup fail, unrecognized name \""+name+"\"")

    def save(self,filename):
        assert isinstance(filename,str), "filename must be a string"

        print "saving trajectory as \"%s\"" % filename
        f=open(filename,'w')
        pickle.dump(self,f)
        f.close()

    def saveMat(self,filename,dataname='rawesomeTrajectory'):
        assert isinstance(filename,str), "filename must be a string"
        ret = {}
        for name in self.dvMap._xNames + \
            self.dvMap._zNames + \
            self.dvMap._uNames + \
            self.outputMap._outputNames + \
            self.quadratureMap._quadMap.keys():
            (ts,ys) = self.getTimeSeries(name)
            ts = numpy.array(ts)
            ys = numpy.array(ys)
            ret[name] = {'time':ts, 'value':ys}

        print "saving trajectory as matlab file \"%s\"" % filename
        scipy.io.savemat(filename, {dataname:ret})

    def getTimeSeries(self,name):
        # if it's a differential state
        if name in self.dvMap._xNames:
            ys = []
            ts = []
            for i in range(self.dvMap._nk):
                for k in range(self.dvMap._nicp):
                    for j in range(self.dvMap._deg+1):
                        ys.append(self.dvMap.lookup(name,timestep=i,nicpIdx=k,degIdx=j))
                        ts.append(self.tgrid[i,k,j])
                    ys.append(numpy.nan*ys[-1])
                    ts.append(ts[-1])
            ys.append(self.dvMap.lookup(name,timestep=-1,nicpIdx=0,degIdx=0))
            ts.append(self.tgrid[-1,0,0])
            return (ts,ys)

        # if it's an algebraic var
        elif name in self.dvMap._zNames:
            ys = []
            ts = []
            for i in range(self.dvMap._nk):
                for k in range(self.dvMap._nicp):
                    for j in range(1,self.dvMap._deg+1):
                        ys.append(self.dvMap.lookup(name,timestep=i,nicpIdx=k,degIdx=j))
                        ts.append(self.tgrid[i,k,j])
                    ys.append(numpy.nan*ys[-1])
                    ts.append(ts[-1])
            return (ts,ys)

        # if it's a control
        elif name in self.dvMap._uNames:
            ys = []
            ts = []
            for i in range(self.dvMap._nk):
                y = self.dvMap.lookup(name,timestep=i)
                t0 = self.tgrid[i,0,0]
                t1 = self.tgrid[i+1,0,0]
                ys.extend([y,y,numpy.nan*y])
                ts.extend([t0,t1,t1])
            return (ts,ys)

        # if it's an output defined everywhere
        elif name in self.outputMap._outputs0:
            ys = []
            ts = []
            for i in range(self.dvMap._nk):
                for k in range(self.dvMap._nicp):
                    for j in range(self.dvMap._deg+1):
                        ys.append(self.outputMap.lookup(name,timestep=i,nicpIdx=k,degIdx=j))
                        ts.append(self.tgrid[i,k,j])
                    ys.append(numpy.nan*ys[-1])
                    ts.append(ts[-1])
            return (ts, ys)

        # if it's an output defined only on collocation points
        elif name in self.outputMap._outputs:
            ys = []
            ts = []
            for i in range(self.dvMap._nk):
                for k in range(self.dvMap._nicp):
                    for j in range(1,self.dvMap._deg+1):
                        ys.append(self.outputMap.lookup(name,timestep=i,nicpIdx=k,degIdx=j))
                        ts.append(self.tgrid[i,k,j])
                    ys.append(numpy.nan*ys[-1])
                    ts.append(ts[-1])
            return (ts, ys)

        # if it's a quadrature state
        elif name in self.quadratureMap._quadMap:
            ys = []
            ts = []
            for i in range(self.quadratureMap._nk):
                for k in range(self.quadratureMap._nicp):
                    for j in range(self.quadratureMap._deg+1):
                        ys.append(self.quadratureMap.lookup(name,timestep=i,nicpIdx=k,degIdx=j))
                        ts.append(self.tgrid[i,k,j])
                    ys.append(numpy.nan*ys[-1])
                    ts.append(ts[-1])
            return (ts, ys)

        # throw error on parameter
        elif name in self.dvMap._pNames:
            raise Exception("can't plot a parameter (\""+name+"\")")

        # throw error on unrecognized
        else:
            raise Exception("unrecognized name \""+name+"\"")


# thing which facilitates plotting
class TrajectoryPlotter(Trajectory):
    """
    A Trajectory that can plot itself
    """
    def subplot(self,names,title=None,style=None):
        assert isinstance(names,list)

        fig = plt.figure()
        if title is None:
            if isinstance(names,str):
                title = names
            else:
                assert isinstance(names,list)
                if len(names) == 1:
                    title = names[0]
                else:
                    title = str(names)
        fig.canvas.set_window_title(str(title))

        plt.clf()
        n = len(names)
        if style is None:
            style = [None]*n
        for k,name in enumerate(names):
            plt.subplot(n,1,k+1)
            if k==0:
                self._plot(name,title,style=style[k])
            else:
                self._plot(name,None,style=style[k])

    def plot(self,names,title=None,style=None):
        fig = plt.figure()
        if title is None:
            if isinstance(names,str):
                title = names
            else:
                assert isinstance(names,list)
                if len(names) == 1:
                    title = names[0]
                else:
                    title = str(names)

        fig.canvas.set_window_title(str(title))
        plt.clf()
        self._plot(names,title,style=style)

    def _plot(self,names,title,style=None,showLegend=True):
        if isinstance(names,str):
            names = [names]
        assert isinstance(names,list)

        legend = []
        for name in names:
            assert isinstance(name,str)
            legend.append(name)
            (ts,ys) = self.getTimeSeries(name)
            if style is None:
                plt.plot(ts,ys)
            else:
                plt.plot(ts,ys,style)

        if title is not None:
            assert isinstance(title,str), "title must be a string"
            plt.title(title)
        plt.xlabel('time [s]')
        if showLegend is True:
            plt.legend(legend)
        plt.grid()
