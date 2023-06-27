import os
import pathlib
import shutil
import subprocess
import sys

#Whether to convert_md using markdown library (True), or let github do it (False)
convert_md=True
rtd_theme=False #Not currently working

# Set the package name, subdirectory, and output directory
subdir = '.\w4h'
output_dir = 'docs'

venvPath = pathlib.Path(sys.executable).parent.parent

os.environ['PYTHONPATH'] = '..' + os.pathsep + os.environ.get('PYTHONPATH', '')

# Run the pdoc command
if rtd_theme:
    themePath = venvPath.as_posix()+'/lib/site-packages/sphinx_rtd_theme/'
    subprocess.run(['pdoc', '--html', '-o', output_dir, '--force', '--template-dir', themePath, subdir])
else:
    subprocess.run(['pdoc', '--html', '-o', output_dir, '--force', subdir ])

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

keepList = ['generate_docs', 'conf']

for f in trg_path.iterdir():
    if f.stem in keepList or f.is_dir():
        pass
    else:
        os.remove(f)

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
            
            html = html.replace('```mermaid', '<pre class="mermaid">')
            html = html.replace('```', '</pre>')
            html = html.replace("direction RL</p>", "direction RL")
            html = html.replace("<pre><code>    A0((file_setup))", "    A0((file_setup))")
            html = html.replace("</code></pre>", "</pre>")
            html = html.replace("<p></pre></p>", '')
            html = html.replace("&gt;", '>')

            dst = pathlib.Path('index.html')
            with open(dst, 'w') as f:
                f.write("<head>\n")
                #f.write('\t<script src="https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.min.js"></script>\n')
                f.write('<script type="module">')
                f.write("\timport mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';\n")
                f.write("\tmermaid.initialize({ startOnLoad: true });\n")
                f.write("</script>\n")
                f.write('<script>console.log(mermaid.version);</script>\n')
                f.write("</head>\n")
                f.write(html)
            print(dst)
            break
        else:
            #Copy main readme file into docs so github pages will read it
            shutil.copy(src=str(each_file), dst='.')
            break
print('docs updated')