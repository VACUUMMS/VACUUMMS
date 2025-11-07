#!/usr/bin/env python

import vacuumms as v
import os

gfg_path = os.environ["TestDataPath"] + "/PMP.lmps"
print("gfg_path: " + gfg_path)

gfg = v.LAMMPSConfiguration(gfg_path)
print(gfg)

#g.cram()
#p = v.Parameters(["-n", "5"])
#c=v.DDX(g, p)
#c.execute()
#print(c.getResult())


