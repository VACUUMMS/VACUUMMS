#!/bin/bash

echo "Example output to file" > test_compare_output.out

# diff returns 0 when files match, and shell returns code from last process.
#diff  ${TestDataPath}/test/expected_compare_output.out
diff  test_compare_output.out ../../test/expected_compare_output.out

