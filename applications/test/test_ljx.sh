#!/bin/sh

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

${BinaryPath}/cmake/applications/ljx -ng -N 27 -end_mcs 10 -energy_report_frequency 1 |grep -v "started" | grep -v "finished" > test_ljx.out

# If the files match, the exit code is zero, and the shell script will return that code
diff test_ljx.out ${TestDataPath}/applications/test/expected_ljx.out

