#!/bin/sh

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

${BinaryPath}/applications/ddx < ${TestDataPath}/fcc.gfg -box 4.242640687119285 4.242640687119285 4.242640687119285 -n 10 > test_ddx_fcc.out

diff test_ddx_fcc.out ${TestDataPath}/expected_ddx_fcc.out

