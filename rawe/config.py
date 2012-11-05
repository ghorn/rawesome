from configobj import ConfigObj,flatten_errors
from validate import Validator

def readConfig(filename,specfilename):
    assert(isinstance(filename,str))
    assert(isinstance(specfilename,str))
    
    conf = ConfigObj(filename, configspec=specfilename)
    results = conf.validate(Validator())
    if results != True:
        print 'Configuration validation error'
        print conf
        print results
        def niceErr(confdict,resultsdict):
            for (section_list, key, _) in flatten_errors(confdict, resultsdict):
                print (section_list,key)
                if key is not None:
                    print 'The "%s" key in the section "%s" failed validation' % (key, ', '.join(section_list))
                elif isinstance(confdict[key],dict):
                    niceErr(confdict[key],resultsdict[key])
                else:
                    print 'The following section was missing:%s ' % key
        niceErr(conf,results)
        raise ValueError('error loading config file "'+filename+'" with specfile "'+specfilename+'"')

    return conf
