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

import casadi as C

replace0 = {'real':'IntermediateState',
            'work':'_work',
            'init':'',
            'output':'_output'
            }
def writeAcadoAlgorithm(dae, f):
    xdot = C.veccat([dae.ddt(name) for name in dae.xNames()])
    inputs = [dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec(), xdot]
    outputs = [dae[name] for name in dae.outputNames()]
    f = C.SXFunction( inputs, [f]+outputs )
    f.init()

    inputNames = [dae.xNames(), dae.zNames(), dae.uNames(), dae.pNames(), ['dot( '+name+' )' for name in dae.xNames()]]

    # error dictionary
    errorMap = {}
    for key,value in C.__dict__.iteritems():
        if key[0:3] == "OP_":
            errorMap[value] = key

    algStrings = []
    initializedWorkVars = set();
    def write(blah):
        algStrings.append(blah)
    def makeUnary(op, replace):
        replace['op'] = op
        write( '%(init)s%(work)s_%(i1)d = %(op)s( %(work)s_%(i2)d );' %  replace)
    def makeInfixBinary(op, replace):
        replace['op'] = op
        write( '%(init)s%(work)s_%(i1)d = %(work)s_%(i2)d %(op)s %(work)s_%(i3)d;' %  replace)
    def makePrefixBinary(op, replace):
        replace['op'] = op
        write( '%(init)s%(work)s_%(i1)d = %(op)s( %(work)s_%(i2)d, %(work)s_%(i3)d );' %  replace)
    # Loop over the algorithm
    for i in range(f.getAlgorithmSize()):

        # Get the atomic operation
        op = f.getAtomicOperation(i)
        i1 = f.getAtomicOutput(i)

        replace = dict(replace0.items() + {'i1':i1}.items())
        if op != C.OP_OUTPUT:
            if i1 not in initializedWorkVars:
                initializedWorkVars.add(i1)
                replace['init'] = replace['real']+' '
        if(op==C.OP_CONST):
            replace['const'] = repr(f.getAtomicInputReal(i))
            write( '%(init)s%(work)s_%(i1)d = %(const)s;' % replace )
        else:
            i2,i3 = f.getAtomicInput(i)
            replace['i2'] = i2
            replace['i3'] = i3
            if op==C.OP_INPUT:
                #assert i2==0, "oh noes, INPUT IS MULTIDIMENSIONAL!!!"
                replace['input'] = inputNames[i2][i3]
                write( '%(init)s%(work)s_%(i1)d = %(input)s;' % replace)
            elif op==C.OP_OUTPUT:
#                assert i1==0, "oh noes, OUTPUT IS MULTIDIMENSIONAL!!!"
#                write( '%(spaces)sf << 0 == %(work)s_%(i2)d;' % replace )
                rowidx = f.output(i1).sparsity().getRow()[i3]
                colidx = f.output(i1).sparsity().col()[i3]
                if i1==0:
                    assert colidx==0
                    replace['row'] = rowidx
                    write( '_f << 0 == %(work)s_%(i2)d; /* %(row)dth residual */' % replace )
                else:
                    replace['i1'] = i1-1
                    replace['outputname'] = dae.outputNames()[i1-1]
                    write( '%(real)s %(output)s_%(i1)d_%(i3)d = %(work)s_%(i2)d; /* Output "%(outputname)s" */' % replace )

            ########## BINARY ########
            elif op==C.OP_ADD:
                makeInfixBinary('+',replace)
            elif op==C.OP_SUB:
                makeInfixBinary('-',replace)
            elif op==C.OP_MUL:
                makeInfixBinary('*',replace)
            elif op==C.OP_SQ:
                replace['i3'] = replace['i2']
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

def generateAcadoCodegenModel(dae,f):
    lines = []
    lines.append('/* Differential States */')
    for name in dae.xNames():
        lines.append('DifferentialState '+name+';')
    lines.append('')
    lines.append('/* Parameters posing as Differential States */')
    for name in dae.pNames():
        lines.append('DifferentialState '+name+'; /* REALLY A PARAMETER */')
    lines.append('')
    lines.append('/* Algebraic States */')
    for name in dae.zNames():
        lines.append('AlgebraicState '+name+';')
    lines.append('')
    lines.append('/* Control Inputs */')
    for name in dae.uNames():
        lines.append('Control '+name+';')
    lines.append('')

    lines.append('/* The Differential Equation */')
    lines.append('DifferentialEquation _f;')

    # dae residual
    lines.extend( writeAcadoAlgorithm(dae, f) )
    lines.append('')
    lines.append('/* hack to make parameters differential states: */')
    for p in dae.pNames():
        lines.append('_f << 0 == dot( '+p+' );')
    lines.append('')
    lines.append('/* the outputs */')
    lines.append('OutputFcn _h;')
    for k,name in enumerate(dae.outputNames()):
        size = 1
        if hasattr(dae[name], 'size'):
            size = dae[name].size()
        for j in range(size):
            replace = dict(replace0.items() + {'k':k,'j':j,'outputname':name}.items())
            lines.append('_h << %(output)s_%(k)d_%(j)d; /* %(outputname)s (%(j)d)*/' % replace)
    lines.append('')

    return ('\n'.join(lines))
