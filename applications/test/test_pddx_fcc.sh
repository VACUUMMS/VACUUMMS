#!/bin/sh

echo "source>>"${BinaryPath}
echo "project>>"${TestDataPath}

${BinaryPath}/cmake/applications/pddx < ${TestDataPath}/applications/ddx/fcc.gfg -box 4.242640687119285 4.242640687119285 4.242640687119285 -n_samples 10 -n_threads 2 -n_steps 1000 | sort > pddx_fcc.out

diff pddx_fcc.out ${TestDataPath}/applications/test/expected_pddx_fcc.out

