#!/bin/sh

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

${BinaryPath}/cmake/utils/uniq < ${TestDataPath}/applications/ddx/fcc.gfg -box 4.242640687119285 4.242640687119285 4.242640687119285 | sort > uniq.out

diff uniq.out ${TestDataPath}/utils/test/uniq.cav

