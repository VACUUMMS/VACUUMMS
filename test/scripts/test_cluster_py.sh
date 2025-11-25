#!/bin/sh

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

export TestDataPath

${SourcePath}/test/python/test_cluster.py > test_cluster_py.out

# If the files match, the exit code is zero, and the shell script will return that code
diff test_cluster_py.out ${TestDataPath}/expected_cluster_py.out

