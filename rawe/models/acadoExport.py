
import casadi as C

def convertAlgorithm(funname, f):
    # error dictionary
    errorMap = {}
    for key,value in C.__dict__.iteritems():
        if key[0:3] == "OP_":
            errorMap[value] = key

    inputStr = 'input'
    outputStr = 'output'
    realStr = 'real_t'
    spacesStr = '    '
    replace0 = {'spaces':spacesStr, 'real':realStr, 'work':'work','init':'',
                'input':inputStr, 'output':outputStr, 'funname':funname}

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

    write( 'void %(funname)s( %(real)s * %(input)s, %(real)s * %(output)s ){' % replace0 )
    # Loop over the algorithm
    for i in range(f.getAlgorithmSize()):
      
        # Get the atomic operation
        op = f.getAtomicOperation(i)
        i1 = f.getAtomicOutput(i)

        replace = dict(replace0.items() + {'i1':i1}.items())
        if op != C.OUTPUT:
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
                assert i1==0, "oh noes, INPUT IS MULTIDIMENSIONAL!!!"
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
    write( '}' )
    return algStrings

def generateModel(dae,ag):

    dummyParamDots = C.veccat( [C.ssym(p+'__PARAMDOT_DUMMY_ACADO_DOESNT_HANDLE_PARAMS_UUUUGH')
                                for p in dae.pNames()] )
    inputs = C.veccat([ag['x'], ag['p'], ag['z'], ag['u'], ag['xdot'], dummyParamDots])

    # dae residual
    f = C.veccat([ag['f'], dummyParamDots])
    rhs = C.SXFunction( [inputs], [f] )
    rhs.init()
    rhs_string = convertAlgorithm('rhs', rhs)

    # dae residual jacobian
    jf = C.veccat( [ C.jacobian(f,inputs) ] )
    rhs_jacob = C.SXFunction( [inputs], [jf] )
    rhs_jacob.init()
    rhs_jacob_string = convertAlgorithm('rhs_jac', rhs_jacob)

    # outputs
    o = C.veccat( [dae[outname] for outname in dae.outputNames()] )
    outputs = C.SXFunction( [inputs], [o] )
    outputs.init()
    outputs_string = convertAlgorithm('out', outputs)

    # outputs jacobian
    jo = C.veccat( [ C.jacobian(o,inputs) ] )
    outputs_jacob = C.SXFunction( [inputs], [jo] )
    outputs_jacob.init()
    outputs_jacob_string = convertAlgorithm('out_jac', outputs_jacob)

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
    return '\n'.join(modelFile)

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
    modelFile = generateModel(dae, ag)
    simExportFile = generateSimExport(dae)

    return (modelFile, simExportFile)
