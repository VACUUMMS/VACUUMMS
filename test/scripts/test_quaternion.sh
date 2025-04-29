#!/bin/bash

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

${BinaryPath}/libraries/vacuumms_variational/test/test_quaternion > test_quaternion.out

diff test_quaternion.out ${TestDataPath}/expected_quaternion.out


