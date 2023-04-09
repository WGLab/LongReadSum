"""
__main__.py:
Call the command-line interface.
"""

import os
import sys
from os.path import dirname, abspath

# Set the library path
lib_dir = os.path.join(os.environ['CONDA_PREFIX'], "lib")
sys.path.append(lib_dir)

# Set the HDF5 plugin path
parent_dir = dirname(dirname(abspath(__file__)))
os.environ['HDF5_PLUGIN_PATH'] = os.path.join(parent_dir, "lib")

# Import the command-line interface
import cli


def main():
    cli.main()


if __name__ == '__main__':
    main()
