#!/usr/bin/env python

import vacuumms as v
import os

g = v.Configuration( os.environ["TestDataPath"] + "/fcc.gfg")
g.setBoxDimensions([4.242640687119285,4.242640687119285,4.242640687119285])
g.cram()
p = v.Parameters(["-n", "5"])
c=v.DDX(g, p)
c.execute()
print(c.getResult())


