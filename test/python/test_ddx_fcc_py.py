#!/usr/bin/env python

import vacuumms as v
import os

gfg = v.Configuration( os.environ["TestDataPath"] + "/fcc_unit_cell.gfg")

gfg.setBoxDimensions([1.0, 1.0, 1.0])

print("Dumping gfg:")
print(gfg)

gfg.replicate([3,3,3])

print("Dumping gfg again:")
print(gfg)

ddx=v.DDX(gfg)
ddx.setNumberOfSamples(5)

ddx.setNumberOfSteps(50)  # Maximum number of steps before giving up and accepting result without explicit convergence
ddx.setLearningRate(0.01) # Default is 0.01, sets learning rate for gradient descent of location of center
ddx.setTolerance(10.0)    # 10.0 is default, comparison value to derivative to determine conversion
ddx.setRNGSeed(555)       # Set random number generator seed
ddx.setVerletCutoff(64.0) # default is 10.0 * 10.0 = 100.0, cutoff radius^2 when building Verlet list
ddx.setVerletExtent(3)    # how many levels of mirror boxes to search when building Verlet list, default is 1

ddx.execute()
print(ddx.getResult())

