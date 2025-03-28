#!/usr/bin/env python

import vacuumms as v
import os

g = v.LAMMPSConfiguration( os.environ["TestDataPath"] + "/PMP.lmps")
g.cram()
p = v.Parameters(["-n", "5"])
c=v.DDX(g, p)
c.execute()
print(c.getResult())


