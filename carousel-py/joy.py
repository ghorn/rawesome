import pygame
import sys
import os
import time

class Joy():
    def __init__(self,idx=0):
        pygame.init()
        pygame.joystick.init()
        self.js = pygame.joystick.Joystick(idx)
        self.js.init()

        self.oldstdout = os.dup(1)
        self.devnull = os.open('/dev/null', os.O_WRONLY)

    def getAxes(self):
        return self.redirectStdout(self._getAxes)
    def getButtons(self):
        return self.redirectStdout(self._getButtons)
    def getHats(self):
        return self.redirectStdout(self._getHats)

    def _getAxes(self):
        pygame.event.pump()
        return [self.js.get_axis(k) for k in range(0,self.js.get_numaxes())] 

    def _getButtons(self):
        pygame.event.pump()
        return [self.js.get_button(k) for k in range(0,self.js.get_numbuttons())] 
    
    def _getHats(self):
        pygame.event.pump()
        return [self.js.get_hat(k) for k in range(0,self.js.get_numhats())] 
    

    def redirectStdout(self,action):
        sys.stdout.flush() # <--- important when redirecting to files
        
        os.dup2(self.devnull, 1)
#        os.close(devnull)
        ret = action()
        os.dup2(self.oldstdout,1)
#        sys.stdout = os.fdopen(newstdout, 'w')


        return ret

if __name__=='__main__':
    pygame.init()
    j = Joy(0)
    try:
        while True:
            print "axes:    "+str(j.getAxes())
            print "buttons: "+str(j.getButtons())
            print "hats:    "+str(j.getHats())
            print ""
            time.sleep(0.5)
    except KeyboardInterrupt:
        print ""
        pygame.quit()
        pass
