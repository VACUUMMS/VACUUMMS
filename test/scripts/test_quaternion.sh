#!/bin/bash

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

${BinaryPath}/modules/variational/test/test_quaternion > test_quaternion.out

diff test_quaternion.out ${TestDataPath}/expected_quaternion.out


