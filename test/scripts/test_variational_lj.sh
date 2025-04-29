#!/bin/bash

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

${BinaryPath}/libraries/vacuumms_variational/test/test_variational_lj -filename ${TestDataPath}/ljx.gfg > test_variational_lj.out

diff test_variational_lj.out ${TestDataPath}/expected_variational_lj.out


