import os
import subprocess
import pathlib 

# Set the package name, subdirectory, and output directory
subdir = '.\w4h'
output_dir = 'docs'

os.environ['PYTHONPATH'] = '..' + os.pathsep + os.environ.get('PYTHONPATH', '')

# Run the pdoc command
subprocess.run(['pdoc', '--html', subdir, '-o', output_dir, '--force'])

workDir = pathlib.Path(os.getcwd())
if workDir.stem == 'docs':
    pass
elif 'docs' in str(workDir):
    pass
else:
    for p in pathlib.Path('.').rglob('*'): 
        if 'docs\\' in str(p):
            docsPath = pathlib.Path(str(p)[:str(p).find('docs')+ 4])
            os.chdir(docsPath)
            break

src_path = pathlib.Path('./w4h')
trg_path = src_path.parent # gets the parent of the folder 

for each_file in src_path.glob('*.*'): # grabs all files
    each_file.rename(trg_path.joinpath(each_file.name)) # moves to parent folder.

os.rmdir('./w4h')

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