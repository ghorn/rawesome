
integratorTypes = \
    ['INT_EX_EULER',  # Explicit Euler method.
     'INT_RK2',       # Explicit Runge-Kutta integrator of order 2.
     'INT_RK3',       # Explicit Runge-Kutta integrator of order 3.
     'INT_RK4',       # Explicit Runge-Kutta integrator of order 4.
     'INT_IRK_GL2',   # Gauss-Legendre integrator of order 2 (Continuous output Implicit Runge-Kutta).
     'INT_IRK_GL4',   # Gauss-Legendre integrator of order 4 (Continuous output Implicit Runge-Kutta).
     'INT_IRK_GL6',   # Gauss-Legendre integrator of order 6 (Continuous output Implicit Runge-Kutta).
     'INT_IRK_GL8',   # Gauss-Legendre integrator of order 8 (Continuous output Implicit Runge-Kutta).
     
     'INT_IRK_RIIA1', # Radau IIA integrator of order 1 (Continuous output Implicit Runge-Kutta).
     'INT_IRK_RIIA3', # Radau IIA integrator of order 3 (Continuous output Implicit Runge-Kutta).
     'INT_IRK_RIIA5', # Radau IIA integrator of order 5 (Continuous output Implicit Runge-Kutta).
     
     'INT_DIRK3',     # Diagonally Implicit 2-stage Runge-Kutta integrator of order 3 (Continuous output).
     'INT_DIRK4',     # Diagonally Implicit 3-stage Runge-Kutta integrator of order 4 (Continuous output).
     'INT_DIRK5',     # Diagonally Implicit 5-stage Runge-Kutta integrator of order 5 (Continuous output).
     
     'INT_DT',        # An algorithm which handles the simulation and sensitivity generation for a discrete time state-space model.
     'INT_NARX']       # An algorithm which handles the simulation and sensitivity generation for a NARX model.


class RtIntegratorOptions(object):
    def __init__(self, optionList=[]):
        # make sure optionList :: [(string, something)]
        usage = '''\
RtIntegratorOptions takes a list of options on construction
Each option must be a (key, value) tuple where key is a string.
'''
        assert isinstance(optionList,list), usage
        for opt in optionList:
            assert isinstance(opt, tuple), usage
            assert len(opt) == 2, usage

        # loop through user options and make sure they're valid and un-repeated
        self._options = {}
        self._defaultOptions = {'NUM_INTEGRATOR_STEPS':1,
                                'LINEAR_ALGEBRA_SOLVER':'GAUSS_LU'}
        for (name,value) in optionList:
            self[name] = value

    def __setitem__(self,name,value):
        assert type(name) == str, 'integrator option key must be a string, you gave: '+str(name)
        assert name not in self._options, \
            'Integrator option "'+name+'" set twice, '+\
            'old val: '+str(self._options[name])+', new val: '+str(value)
        
        def assertType(desiredType):
            assert type(value) == desiredType, 'integrator option "'+name+'" must be '+\
                str(desiredType)+', you gave "'+str(value)+'" which is '+str(type(value))
        if name == 'INTEGRATOR_TYPE':
            assertType(str)
            assert value in integratorTypes, \
                name+' must be one of: '+str(integratorTypes)+'\nyou gave: '+str(value)
            self._options[name] = value

        elif name in ['NUM_INTEGRATOR_STEPS','IMPLICIT_INTEGRATOR_NUM_ITS',
                      'IMPLICIT_INTEGRATOR_NUM_ITS_INIT']:
            assertType(int)
            self._options[name] = value

#        elif name == 'MEASUREMENT_GRID':
#            assertType(str)
#            legalGrids = ['EQUIDISTANT_SUBGRID','EQUIDISTANT_GRID','ONLINE_GRID']
#            assert value in legalGrids, name+' must be one of '+str(legalGrids)+\
#                ', you gave: '+str(value)
#            self._options[name] = value

        elif name == 'LINEAR_ALGEBRA_SOLVER':
            assertType(str)
            legalLAS = ['GAUSS_LU', 'HOUSEHOLDER_QR']
            assert value in legalLAS, name+' must be one of '+str(legalLAS)+\
                ', you gave: '+str(value)
            self._options[name] = value

        else:
            raise Exception('unrecognized integrator option "'+name+'"')

    def __contains__(self, name):
        return name in self._options

    def __getitem__(self,name):
        if name in self._options:
            return self._options[name]
        elif name in self._defaultOptions:
            return self._defaultOptions[name]
        raise Exception('integrator option "'+str(name)+'" has not been set and does not have a default')
