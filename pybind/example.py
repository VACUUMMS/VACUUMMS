#!/usr/bin/env python

import vacuumms as v
c = v.Configuration('fcc.gfg')
p = v.Parameters(["-n", "100", "-box", "4.24264", "4.24264", "4.24264"])
o = v.DDX(c, p)
o.execute()
result=o.getResult()
csdp=v.Parameters(['-width', '0.1', '-n_bins', '50'])
csd=v.CavitySizeDistribution(result, csdp)
csd

csd.writeToFile('foo.bar')
