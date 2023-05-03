import os
import subprocess

# Set the package name, subdirectory, and output directory
subdir = '.\w4h'
output_dir = 'docs'

os.environ['PYTHONPATH'] = '..' + os.pathsep + os.environ.get('PYTHONPATH', '')

# Run the pdoc command
subprocess.run(['pdoc', '--html', subdir, '-o', output_dir, '--force'])


#This script generates simple documentation using the pdoc3 library that can easily be published on github pages
#import os
#import subprocess

# Set the package name, subdirectory, and output directory
#mod = 'invplot.py'
#output_dir = 'docs'

#os.environ['PYTHONPATH'] = '..' + os.pathsep + os.environ.get('PYTHONPATH', '')

# Run the pdoc command
#subprocess.run(['pdoc3', '--html', mod, '-o', output_dir, '--force'])
#os.chdir('./docs')
#os.rename(mod[:-3]+'.html', 'index.html')