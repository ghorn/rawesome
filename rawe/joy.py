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

    def _getAxes(self):
        return [self.js.get_axis(k) for k in range(0,self.js.get_numaxes())]

    def _getButtons(self):
        return [self.js.get_button(k) for k in range(0,self.js.get_numbuttons())]

    def _getHats(self):
        return [self.js.get_hat(k) for k in range(0,self.js.get_numhats())]

    def _getAll(self):
        js = {}
        js['buttonsDown'] = set()
        js['buttonsUp'] = set()
        for e in pygame.event.get():
            if pygame.event.event_name(e.type) =="JoyButtonDown":
                js['buttonsDown'].add(e.button)
            elif pygame.event.event_name(e.type) =="JoyButtonUp":
                js['buttonsUp'].add(e.button)
        js['axes'] = self._getAxes()
        return js

    def getAll(self):
        return self.redirectStdout(self._getAll)

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
