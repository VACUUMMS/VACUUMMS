#!/bin/bash

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

${BinaryPath}/modules/variational/test/test_variational_lj -filename ${TestDataPath}/ljx.gfg > test_variational_lj.out

diff test_variational_lj.out ${TestDataPath}/expected_variational_lj.out


