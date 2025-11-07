import os

# Get the value as a string (or None if not set)
pythonpath = os.environ.get('PYTHONPATH')

if pythonpath:
    print(f"PYTHONPATH: {pythonpath}")
    # It's a colon-separated list of paths on Unix/macOS, semicolon-separated on Windows
    paths = pythonpath.split(os.pathsep)
    print(f"Paths: {paths}")
else:
    print("PYTHONPATH is not set.")

ld_library_path = os.environ.get('LD_LIBRARY_PATH')

if ld_library_path:
    print(f"LD_LIBRARY_PATH: {ld_library_path}")
    # It's a colon-separated list of paths on Unix/macOS, semicolon-separated on Windows
    paths = ld_library_path.split(os.pathsep)
    print(f"Paths: {paths}")
else:
    print("LD_LIBRARY_PATH is not set.")

import vacuumms as v
print("import vacuumms success.")
