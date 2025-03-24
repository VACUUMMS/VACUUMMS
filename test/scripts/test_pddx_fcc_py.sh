#!/bin/sh

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

export TestDataPath

${SourcePath}/test/python/test_pddx_fcc_py.py | sort > test_pddx_fcc_py.out

# If the files match, the exit code is zero, and the shell script will return that code
diff test_pddx_fcc_py.out ${TestDataPath}/expected_pddx_fcc_py.out

