import os
import subprocess

# Set the package name, subdirectory, and output directory
subdir = '..\w4h'
output_dir = '.'

# Run the pdoc command
subprocess.run(['pdoc', subdir, '-o', output_dir])
