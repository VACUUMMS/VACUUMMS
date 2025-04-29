#!/bin/bash

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

${BinaryPath}/libraries/vacuumms_variational/test/test_variational_3D -filename ${TestDataPath}/ljx.gfg > test_variational_3D.out

diff test_variational_3D.out ${TestDataPath}/expected_variational_3D.out


