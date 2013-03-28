import casadi as C

class AlgorithmWriter(object):
    def __init__(self):
        self.inputStr = 'input'
        self.outputStr = 'output'
        self.realStr = 'real_t'
        self.spacesStr = '    '

    def writePrototype(self, funname):
        self.replace0 = {'spaces':self.spacesStr, 'real':self.realStr, 'work':'work','init':'',
                         'input':self.inputStr, 'output':self.outputStr, 'funname':funname}
        return 'void %(funname)s( %(real)s * %(input)s, %(real)s * %(output)s ){' % self.replace0

    def convertAlgorithm(self, f):
        # error dictionary
        errorMap = {}
        for key,value in C.__dict__.iteritems():
            if key[0:3] == "OP_":
                errorMap[value] = key
    
        self.replace0 = {'spaces':self.spacesStr, 'real':self.realStr, 'work':'work','init':'',
                         'input':self.inputStr, 'output':self.outputStr}

        algStrings = []
        initializedWorkVars = set();
        def write(blah):
            algStrings.append(blah)
        def makeUnary(op, replace):
            replace['op'] = op
            write( '%(spaces)s%(init)s%(work)s_%(i1)d = %(op)s( %(work)s_%(i2)d );' %  replace)
        def makeInfixBinary(op, replace):
            replace['op'] = op
            write( '%(spaces)s%(init)s%(work)s_%(i1)d = %(work)s_%(i2)d %(op)s %(work)s_%(i3)d;' %  replace)
        def makePrefixBinary(op, replace):
            replace['op'] = op
            write( '%(spaces)s%(init)s%(work)s_%(i1)d = %(op)s( %(work)s_%(i2)d, %(work)s_%(i3)d );' %  replace)

        # Loop over the algorithm
        for i in range(f.getAlgorithmSize()):
          
            # Get the atomic operation
            op = f.getAtomicOperation(i)
            i1 = f.getAtomicOutput(i)
    
            replace = dict(self.replace0.items() + {'i1':i1}.items())
            if op != C.OP_OUTPUT:
                if i1 not in initializedWorkVars:
                    initializedWorkVars.add(i1)
                    replace['init'] = replace['real']+' '
            if(op==C.OP_CONST):
                replace['const'] = repr(f.getAtomicInputReal(i))
                write( '%(spaces)s%(init)s%(work)s_%(i1)d = %(const)s;' % replace )
            else:
                i2,i3 = f.getAtomicInput(i)
                replace['i2'] = i2
                replace['i3'] = i3
                if op==C.OP_INPUT:
                    assert i2==0, "oh noes, INPUT IS MULTIDIMENSIONAL!!!"
                    write( '%(spaces)s%(init)s%(work)s_%(i1)d = %(input)s[%(i3)d];' % replace)
                elif op==C.OP_OUTPUT:
                    assert i1==0, "oh noes, OUTPUT IS MULTIDIMENSIONAL!!!"
                    rowidx = f.output(i1).sparsity().getRow()[i3]
                    colidx = f.output(i1).sparsity().col()[i3]
                    assert colidx == 0, 'colidx != 0'
                    replace['i3'] = rowidx
                    write( '%(spaces)s%(output)s[%(i3)d] = %(work)s_%(i2)d;' % replace )
                
                ########## BINARY ########
                elif op==C.OP_ADD:
                    makeInfixBinary('+',replace)
                elif op==C.OP_SUB:
                    makeInfixBinary('-',replace)
                elif op==C.OP_MUL:
                    makeInfixBinary('*',replace)
                elif op==C.OP_DIV:
                    makeInfixBinary('/',replace)
                elif op==C.OP_ATAN2:
                    makePrefixBinary('atan2',replace)
                elif op in [C.OP_POW,C.OP_CONSTPOW]:
                    makePrefixBinary('pow',replace)
    
                ########## UNARY ########
                elif op==C.OP_ACOS:
                    makeUnary('acos',replace)
                elif op==C.OP_ATAN:
                    makeUnary('atan',replace)
                elif op==C.OP_COS:
                    makeUnary('cos',replace)
                elif op==C.OP_EXP:
                    makeUnary('exp',replace)
                elif op==C.OP_INV:
                    makeUnary('1.0/',replace)
                elif op==C.OP_SIN:
                    makeUnary('sin',replace)
                elif op==C.OP_ACOSH:
                    makeUnary('acosh',replace)
                elif op==C.OP_COSH:
                    makeUnary('cosh',replace)
                elif op==C.OP_FABS:
                    makeUnary('fabs',replace)
                elif op==C.OP_SINH:
                    makeUnary('sinh',replace)
                elif op==C.OP_ATANH:
                    makeUnary('tanh',replace)
                elif op==C.OP_NEG:
                    makeUnary('-',replace)
                elif op==C.OP_LOG:
                    makeUnary('log',replace)
                elif op==C.OP_SQRT:
                    makeUnary('sqrt',replace)
                elif op==C.OP_ASIN:
                    makeUnary('asin',replace)
                elif op==C.OP_ASINH:
                    makeUnary('asinh',replace)
                elif op==C.OP_TAN:
                    makeUnary('tan',replace)
                elif op==C.OP_ERF:
                    makeUnary('erf',replace)
                elif op==C.OP_ERFINV:
                    makeUnary('erfinv',replace)
                elif op==C.OP_SIGN:
                    makeUnary('sign',replace)
                elif op==C.OP_TANH:
                    makeUnary('tanh',replace)
                elif op==C.OP_ASSIGN:
                    makeUnary('',replace)
                elif op==C.OP_PARAMETER:
                    raise KeyError('Oh man, there is a free parameter in your SXFunction')
                else:
                    raise KeyError('Unknown operation: '+ errorMap[op])

        return algStrings


def generateCModel(dae,ag):
    writer = AlgorithmWriter()

    inputs = C.veccat([ag['x'], ag['z'], ag['u'], ag['p'], ag['xdot']])

    # dae residual
    f = ag['f']
    rhs = C.SXFunction( [inputs], [f] )
    rhs.init()
    # handle time scaling
    [f] = rhs.eval([C.veccat([ag['x'], ag['z'], ag['u'], ag['p'], ag['xdot']/ag['timeScaling']])])
    rhs = C.SXFunction( [inputs], [f] )
    rhs.init()
    rhs_string = [writer.writePrototype('rhs')]
    rhs_string.extend(writer.convertAlgorithm(rhs))
    rhs_string.append('}')

    # dae residual jacobian
    jf = C.veccat( [ C.jacobian(f,inputs).T ] )
    rhs_jacob = C.SXFunction( [inputs], [jf] )
    rhs_jacob.init()
    rhs_jacob_string = [writer.writePrototype('rhs_jac')]
    rhs_jacob_string.extend(writer.convertAlgorithm(rhs_jacob))
    rhs_jacob_string.append('}')

    # outputs
    o = C.veccat( [dae[outname] for outname in dae.outputNames()] )
    outputs = C.SXFunction( [inputs], [o] )
    outputs.init()
    outputs_string = [writer.writePrototype('out')]
    outputs_string.extend(writer.convertAlgorithm(outputs))
    outputs_string.append('}')

    # outputs jacobian
    jo = C.veccat( [ C.jacobian(o,inputs).T ] )
    outputs_jacob = C.SXFunction( [inputs], [jo] )
    outputs_jacob.init()
    outputs_jacob_string = [writer.writePrototype('out_jac')]
    outputs_jacob_string.extend(writer.convertAlgorithm(outputs_jacob))
    outputs_jacob_string.append('}')

    # model file
    modelFile = ['#include "acado.h"']
    modelFile.append('')
    modelFile.extend(rhs_string)
    modelFile.append('')
    modelFile.extend(rhs_jacob_string)
    modelFile.append('')
    modelFile.append('')
    modelFile.extend(outputs_string)
    modelFile.append('')
    modelFile.extend(outputs_jacob_string)
    return {'modelFile':'\n'.join(modelFile),
            'rhs':rhs,
            'rhsJacob':rhs_jacob}



def generateSimExport(dae):
    return \
    '\n'.join(['DifferentialState '+xn+';' for xn in dae.xNames()]) + \
    '\n'+\
    '\n'.join(['DifferentialState '+pn+'; /* really a parameter */' for pn in dae.pNames()]) + \
    '\n'+\
    '\n'.join(['AlgebraicState '+zn+';' for zn in dae.zNames()]) + \
    '\n'+\
    '\n'.join(['Control '+un+';' for un in dae.uNames()]) + \
    '''

SIMexport sim();
sim.set( INTEGRATOR_TYPE, INT_IRK_RIIA3 );
sim.set( NUM_INTEGRATOR_STEPS, 4 );
sim.set( MEASUREMENT_GRID, EQUIDISTANT_GRID );

sim.setModel( "model", "rhs", "rhs_jac" );
sim.setDimensions( %d, %d, %d, %d );

sim.addOutput( "out", "out_jac", 2 );
sim.setMeasurements( Meas );
sim.setTimingSteps( 10000 );
sim.exportAndRun( "externModel_export", "init_externModel.txt", "controls_externModel.txt" );
''' % (len(dae.xNames()+dae.pNames()), len(dae.xNames()+dae.pNames()), len(dae.zNames()), len(dae.uNames()))

def simExport(dae, ag):
    modelFile = generateCModel(dae, ag)
    simExportFile = generateSimExport(dae)

    return (modelFile, simExportFile)
