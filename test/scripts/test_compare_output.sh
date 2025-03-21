#!/bin/bash

echo "Example output to file" > test_compare_output.out

diff  test_compare_output.out ${TestDataPath}/expected_compare_output.out

