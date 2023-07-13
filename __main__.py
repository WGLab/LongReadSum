"""
__main__.py:
Call the command-line interface.
"""

import os

from src import cli


from os.path import dirname, abspath

# Get the parent directory
parent_dir = dirname(dirname(abspath(__file__)))

# # Set the HDF5 plugin path
# os.environ['HDF5_PLUGIN_PATH'] = os.path.join(parent_dir, "lib")
print("HDF5_PLUGIN_PATH is " + os.environ.get('HDF5_PLUGIN_PATH', ''))


def main():
    cli.main()


if __name__ == '__main__':
    main()
