#!/bin/sh

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

${BinaryPath}/cmake/applications/ddx < ${TestDataPath}/applications/ddx/fcc.gfg -box 4.242640687119285 4.242640687119285 4.242640687119285 -n 10 > ddx_fcc.out

diff ddx_fcc.out ${TestDataPath}/applications/test/expected_ddx_fcc.out

