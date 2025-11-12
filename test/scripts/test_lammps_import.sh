#!/bin/sh

echo "SourcePath>>"${SourcePath}
echo "TestDataPath>>"${TestDataPath}

export TestDataPath
export SourcePath

${SourcePath}/test/python/test_lammps_import.py > test_lammps_import.out

# If the files match, the exit code is zero, and the shell script will return that code
diff test_lammps_import.out ${TestDataPath}/expected_lammps_import.out

