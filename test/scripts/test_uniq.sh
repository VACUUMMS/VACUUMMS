#!/bin/sh

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

${BinaryPath}/utils/uniq < ${SourcePath}/applications/ddx/fcc.gfg -box 4.242640687119285 4.242640687119285 4.242640687119285 | sort > test_uniq.out

diff test_uniq.out ${SourcePath}/test/data/expected_uniq.out



