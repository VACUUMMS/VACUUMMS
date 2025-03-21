#!/bin/sh

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

export TestDataPath

${SourcePath}/test/python/test_lammps_import.py > test_lammps_import.out

# If the files match, the exit code is zero, and the shell script will return that code
diff test_lammps_import.out ${TestDataPath}/expected_lammps_import.out

