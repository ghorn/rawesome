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


class Option(object):
    def __init__(self, name, valueType, default):
        assert type(name) == str, 'option name must be a string'
        self._name = name
        self._valueType = valueType
        self._defaultValue = default

    def getOutputVal(self):
        '''
        Tries to returns the user set value.
        If no user value, tries to return default.
        If no default, throws exception.
        '''
        if hasattr(self, '_value'):
            return self._value
        elif self._defaultValue is not None:
            return self._defaultValue
        raise Exception(self.prettyName()+' is unset')

    def prettyName(self):
        return self._parentName+' option "'+self._name+'"'

    def name(self):
        return self._name

    def _assertType(self,value):
        assert type(value) == self._valueType, self.prettyName()+' must be ' +\
            str(self._valueType)+', you gave '+str(value)+' which is '+str(type(value))

    def setVal(self,value):
        '''
        Assert value is the right type, and that this option has not been set yet.
        Then store the value.
        '''
        self._assertType(value)
        assert not hasattr(self,'_value'), self.prettyName()+' has already been set, '+\
            'old value: '+str(self._value)+', new value: '+str(value)
        self._value = value

class OptStr(Option):
    def __init__(self, name, opts, default=None):
        self._legalOptions = opts
        if default is not None:
            self._assertLegalOption(default)
        Option.__init__(self,name,str,default)

    def _assertLegalOption(self,value):
        assert value in self._legalOptions, self.prettyName()+' must be in '+\
            str(self._legalOptions)+', you gave: '+str(value)

    def setVal(self,value):
        self._assertLegalOption(value)
        Option.setVal(self,value)

    def toAcadoStr(self):
        return self.getOutputVal()

class OptInt(Option):
    def __init__(self, name, default=None):
        Option.__init__(self,name,int,default)

    def toAcadoStr(self):
        return repr(self.getOutputVal())

class OptBool(Option):
    def __init__(self, name, default=None):
        Option.__init__(self,name,bool,default)

    def toAcadoStr(self):
        if self.getOutputVal():
            return 'YES'
        else:
            return 'NO'
        

class Options(object):
    '''
    To use this, subclass Options and call .add() with a bunch of options.
    Then the use __setitem__ to set options.
    '''
    def __init__(self,name):
        self._name = name
        self._options = {}

    def _assertValidName(self,name):
        assert type(name) == str, 'option key must be a string, you gave: '+str(name)
        assert name in self._options, '"'+name+'" is not a valid option\n'+\
            'valid options: '+str(sorted(self._options.keys()))

    def add(self,opt):
        assert isinstance(opt, Option)
        name = opt.name()
        assert name not in self._options, self._name+' already has Option '+str(name)
        opt._parentName = self._name
        self._options[name] = opt

    def __setitem__(self,name,value):
        self._assertValidName(name)
        self._options[name].setVal(value)

    def __getitem__(self,key):
        self._assertValidName(key)
        return self._options[key].getOutputVal()

    def getAcadoOpts(self):
        ret = {}
        for name,opt in self._options.iteritems():
            ret[name] = opt.toAcadoStr()
        return ret
