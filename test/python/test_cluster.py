#!/usr/bin/env python

import vacuumms as v
import os

cavs = v.Configuration( os.environ["TestDataPath"] + "/test.cav")
cavs
cavs.setBoxDimensions([5,5,5])
cavs.generateClusters()
