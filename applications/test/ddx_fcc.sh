#!/bin/sh

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

${BinaryPath}/cmake/applications/ddx < ${TestDataPath}/applications/ddx/fcc.gfg -box 4.242640687119285 4.242640687119285 4.242640687119285 -n 10 > fcc.out

diff fcc.out ${TestDataPath}/applications/test/ddx_fcc.cav

