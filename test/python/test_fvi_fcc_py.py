#!/usr/bin/env python 

import vacuumms as v
import os

g = v.Configuration( os.environ["TestDataPath"] + "/fcc.gfg")
g.setBoxDimensions([4.242640687119285,4.242640687119285,4.242640687119285])
g.cram()

fvi=v.FVIX(g)
fvi.setDimensions([4,4,4])
fvi.execute()

print(fvi.getFVI())

print("done")
