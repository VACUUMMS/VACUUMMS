#!/bin/sh

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

export TestDataPath

${SourcePath}/test/python/test_ddx_fcc_py.py > test_ddx_fcc_py.out

# If the files match, the exit code is zero, and the shell script will return that code
diff test_ddx_fcc_py.out ${TestDataPath}/expected_ddx_fcc_py.out

