import pygame
import sys
import os
import time

def redirect_stdout():
    print "Joy: Redirecting C stdout because current pygame was compiled with DEBUG flags on :("
    sys.stdout.flush() # <--- important when redirecting to files
    newstdout = os.dup(1)
    devnull = os.open('/dev/null', os.O_WRONLY)
    os.dup2(devnull, 1)
    os.close(devnull)
    sys.stdout = os.fdopen(newstdout, 'w')

class Joy():
    def __init__(self,idx=0):
        redirect_stdout()
        pygame.init()
        pygame.joystick.init()
        self.js = pygame.joystick.Joystick(idx)
        self.js.init()

    def getAxes(self):
        pygame.event.pump()
        return [self.js.get_axis(k) for k in range(0,self.js.get_numaxes())] 

    def getButtons(self):
        pygame.event.pump()
        return [self.js.get_button(k) for k in range(0,self.js.get_numbuttons())] 
    
    def getHats(self):
        pygame.event.pump()
        return [self.js.get_hat(k) for k in range(0,self.js.get_numhats())] 
    
#    def update(self):
#    def getAll(self):

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
