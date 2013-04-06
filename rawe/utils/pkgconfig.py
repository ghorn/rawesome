import subprocess
import sys
import os

def call(args):
    # make sure pkg-config is available
    try:
        subprocess.check_call(["pkg-config", '--version'], \
                                  stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    except (subprocess.CalledProcessError, OSError):
        sys.stderr.write('`pkg-config` has not been found, please install it\n')
        sys.exit(os.EX_CONFIG)

    # make sure pkg-config can find acado and ocg2
    try:
        p = subprocess.Popen(['pkg-config']+args, \
                                 stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (output, error) = p.communicate()
        if error != '':
            raise OSError
    except (subprocess.CalledProcessError, OSError):
        raise Exception('\n'+error)

    return output.strip()
