import os
import pathlib 
import shutil
import subprocess

#Whether to convert_md using markdown library (True), or let github do it (False)
convert_md=True

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
    destFilePath = trg_path.joinpath(each_file.name)
    if destFilePath.is_file():
        os.remove(destFilePath)
    each_file = each_file.rename(destFilePath) # moves to parent folder.
    if each_file.name == 'index.html':
        mainhtmlFPath = each_file.parent.parent.joinpath('main.html')
        if mainhtmlFPath.is_file():
            os.remove(mainhtmlFPath)
        each_file.rename(mainhtmlFPath)
os.rmdir('./w4h')

repo_path = pathlib.Path('..')
for each_file in repo_path.iterdir():
    if each_file.name == 'README.md':
        if convert_md:
            import markdown

            with open(each_file, 'r') as f:
                markdown_text = f.read()

            html = markdown.markdown(markdown_text)
            
            dst = pathlib.Path('index.html')
            with open(dst, 'w') as f:
                f.write("<head>\n")
                f.write('\t<script src="https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.min.js"></script>\n')
                f.write('\t<script>console.log(mermaid.version);</script>\n')
                f.write("</head>\n")     
                f.write(html)
            print(dst)
            break
        else:
            #Copy main readme file into docs so github pages will read it
            shutil.copy(src=str(each_file), dst='.')
            break