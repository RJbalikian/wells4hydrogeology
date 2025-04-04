"""This file is meant to be run on its own 
to generate documentation from docstrings and other repository files/variables"""

import os
import pathlib
import re
import shutil
import subprocess
import sys

import markdown

RELEASE_VERSION = "0.1.2"

VERBOSE = True

# Whether to CONVERT_MD using markdown (True), or github (False)
RTD_DOCS = True
GITHUB_PAGES = False  # Don't think I  need this anymore, and it still works

CONVERT_MD = True
RTD_THEME = False  # Not currently working

RUN_TESTS = True
LINT_IT = False

if VERBOSE:
    print('Specifying filepaths')

# Set the filepaths
currentDir = pathlib.Path((__file__)).parent
docsDir = currentDir
repoDir = docsDir.parent
w4hDir = repoDir.joinpath('w4h')

classifyPath = w4hDir.joinpath('classify.py')
cleanPath = w4hDir.joinpath('clean.py')
corePath = w4hDir.joinpath('core.py')
exportPath = w4hDir.joinpath('export.py')
layersPath = w4hDir.joinpath('layers.py')
mappingPath = w4hDir.joinpath('mapping.py')
readPath = w4hDir.joinpath('read.py')

OUTPUT_DIR = docsDir

# Update setup file(s) with release version number
setupFPath = repoDir.joinpath('setup.py')
pyprojectFPath = repoDir.joinpath('pyproject.toml')
confFPath = docsDir.joinpath('conf.py')
condaFPath = repoDir.joinpath('conda/meta.yaml')
initFPath = w4hDir.joinpath('__init__.py')

if VERBOSE:
    print('Updating version/release numbers')
confFilePaths = [setupFPath, pyprojectFPath, condaFPath, confFPath, initFPath]
for cFile in confFilePaths:
    if cFile.exists():
        with open(cFile.as_posix(), mode='r', encoding='utf-8') as f:
            cFileText = f.read()

        # Update which file is being analyzed for version number
        # Intended for setup.py
        VERTEXT = r'version=".*?"'
        NEWVERTEXT = r'version="'+RELEASE_VERSION+'"'
        cFileText = re.sub(VERTEXT, NEWVERTEXT, cFileText, flags=re.DOTALL)

        # intended for pyproject.toml
        VERTEXT = r'version:\s+\d+\.\d+\.\d+[^\n]*'
        NEWVERTEXT = r'version: '+RELEASE_VERSION
        cFileText = re.sub(VERTEXT, NEWVERTEXT, cFileText, flags=re.DOTALL)

        # intended for conda/meta.yaml, if used
        VERTEXT = r'git_tag:\s+v+\d+\.\d+\.\d+[^\n]*'
        NEWVERTEXT = r'git_tag: v'+RELEASE_VERSION
        cFileText = re.sub(VERTEXT, NEWVERTEXT, cFileText, flags=re.DOTALL)

        # intended for conf.py
        VERTEXT = r"release\s+=+\s+'+\d+\.\d+\.\d+'+[^\n]*"
        NEWVERTEXT = r"release = '"+RELEASE_VERSION+"'"
        cFileText = re.sub(VERTEXT, NEWVERTEXT, cFileText, flags=re.DOTALL)

        # intended for __init__.py
        VERTEXT = r'__version__ = ".*?"'
        NEWVERTEXT = r'__version__ = "'+RELEASE_VERSION+'"'
        cFileText = re.sub(VERTEXT, NEWVERTEXT, cFileText, flags=re.DOTALL)

        with open(cFile.as_posix(), mode='w', encoding='utf-8') as f:
            f.write(cFileText)


venvPath = pathlib.Path(sys.executable).parent.parent
os.environ['PYTHONPATH'] = '..' + os.pathsep + os.environ.get('PYTHONPATH', '')

if RTD_DOCS:
    if VERBOSE:
        print('Creating Read the Docs files')
    keepList = ['_generate_docs', 'conf', 'requirements']

    if VERBOSE:
        print('Saving conf.py and index.rst so they do not get overwritten')
    # It seems apidoc rewrites conf.py and index.rst file (don't want that),
    # so save it first and rewrite after
    confFilePath = docsDir.joinpath('conf.py')
    indFilePath = docsDir.joinpath('index.rst')

    with open(confFilePath.as_posix(), mode='r', encoding='utf-8') as f:
        confFileText = f.read()
    with open(indFilePath.as_posix(), mode='r', encoding='utf-8') as f:
        indFileText = f.read()

    if VERBOSE:
        print('Running sphinx-apidoc')
    # Run apidoc to update api documentation from docstrings
    subprocess.run(['sphinx-apidoc', '-F', '-M', '-e', '-f', '-o',
                    docsDir.as_posix(), w4hDir.as_posix()])

    if VERBOSE:
        print('sphinx-apidoc complete')
        print('Putting original conf.py/index.rst files back where I got them')
    with open(confFilePath.as_posix(), mode='w', encoding='utf-8') as f:
        f.write(confFileText)
    with open(indFilePath.as_posix(), mode='w', encoding='utf-8') as f:
        f.write(indFileText)

    if VERBOSE:
        print('Cleaning make file and sources')
    subprocess.run([docsDir.joinpath('make.bat').as_posix(), 'clean'],
                   check=False)
    if VERBOSE:
        print('Creating make file (for html)')
    subprocess.run([docsDir.joinpath('make.bat').as_posix(), 'html'],
                   check=False)

    if VERBOSE:
        print('Organizing documentation files')
    buildDir = docsDir.joinpath('_build')
    htmlDir = buildDir.joinpath('html')

    copyList = ['documentation_options', 'doctools',
                'sphinx_highlight', 'theme', 'pygments']

    for f in htmlDir.iterdir():
        if f.name[0] != '_':
            shutil.copy(f, docsDir.joinpath(f.name))
        else:
            if f.name == '_static':
                for f2 in f.iterdir():
                    if f2.stem in copyList:
                        shutil.copy(f2, docsDir.joinpath(f.name))
                    elif f2.stem == 'js':
                        shutil.copy(f2.joinpath('theme.js'),
                                    docsDir.joinpath('theme.js'))
                    elif f2.stem == 'css':
                        shutil.copy(f2.joinpath('theme.css'),
                                    docsDir.joinpath('theme.css'))

    for file in docsDir.iterdir():
        if file.suffix == '.html':
            with open(file.as_posix(), mode='r', encoding='utf-8') as f:
                htmlFileText = f.read()
            htmlFileText = htmlFileText.replace('src="_static/', 'src="')
            htmlFileText = htmlFileText.replace('src="js/', 'src="')

            htmlFileText = htmlFileText.replace('href="_static/', 'href="')
            htmlFileText = htmlFileText.replace('href="css/', 'href="')
            with open(file.as_posix(), mode='w', encoding='utf-8') as f:
                f.write(htmlFileText)

if GITHUB_PAGES:
    if VERBOSE:
        print('Running the pdoc command (for Github pages documentation)')
    # Run the pdoc command
    if RTD_THEME:
        themePath = venvPath.as_posix()+'/lib/site-packages/sphinx_RTD_THEME/'
        subprocess.run(['pdoc', '--html', '-o', OUTPUT_DIR, '--force',
                        '--template-dir', themePath, w4hDir], check=False)
    else:
        subprocess.run(['pdoc', '--html', '-o', OUTPUT_DIR.as_posix(),
                        '--force', w4hDir.as_posix()], check=False)

    # Not used anymore, probably
    workDir = repoDir
    if workDir.stem == 'docs':
        pass
    elif 'docs' in str(workDir):
        pass
    else:
        for p in pathlib.Path(workDir).rglob('*'):
            if 'docs\\' in str(p):
                docsPath = pathlib.Path(str(p)[:str(p).find('docs') + 4])
                os.chdir(docsPath)
                break

    src_path = docsDir.joinpath('w4h')
    trg_path = docsDir

    if VERBOSE:
        print('Reading .html files from', src_path.absolute())
        print('Placing html files in', trg_path.absolute())

    # Move items back into main docs folder
    keepList = ['_generate_docs', 'conf', 'make',
                'Makefile', 'index', '.readthedocs']

    for f in trg_path.iterdir():
        if f.stem in keepList or f.is_dir():
            if f.stem == 'resources':
                os.remove(f)
        else:
            os.remove(f)

    for each_file in src_path.glob('*.*'):  # grabs all files
        destFilePath = trg_path.joinpath(each_file.name)
        if destFilePath.is_file():
            os.remove(destFilePath)
        each_file = each_file.rename(destFilePath)  # moves to parent folder.
        if each_file.name == 'index.html':
            mainhtmlFPath = docsDir.joinpath('main.html')
            if mainhtmlFPath.is_file():
                os.remove(mainhtmlFPath)
            each_file.rename(mainhtmlFPath)
    os.rmdir(src_path)

    readmePath = repoDir.joinpath("README.md")

    if CONVERT_MD:
        if VERBOSE:
            print('Converting markdown to html')
        with open(readmePath.as_posix(), mode='r', encoding='utf-8') as f:
            markdown_text = f.read()

        html = markdown.markdown(markdown_text)
       
        html = html.replace('```mermaid', '<pre class="mermaid">')
        html = html.replace('```', '</pre>')
        html = html.replace("direction RL</p>", "direction RL")
        html = html.replace("<pre><code>    A0((file_setup))", "    A0((file_setup))")
        html = html.replace("</code></pre>", "</pre>")
        html = html.replace("<p></pre></p>", '')
        html = html.replace("&gt;", '>')

        dst = docsDir.joinpath('index.html')
        with open(dst, mode='w', encoding='utf-8') as f:
            f.write("<head>\n")
            f.write('<script type="module">')
            f.write("\timport mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';\n")
            f.write("\tmermaid.initialize({ startOnLoad: true });\n")
            f.write("</script>\n")
            f.write('<script>console.log(mermaid.version);</script>\n')
            f.write("</head>\n")
            f.write(html)
        if VERBOSE:
            print('HTML landing page:', dst)

    else:
        # Copy main readme file into docs so github pages will read it
        shutil.copy(src=str(readmePath.as_posix()), dst='.')

if LINT_IT:
    if VERBOSE:
        print('Running linting')
    fileList = [classifyPath, cleanPath, corePath, exportPath,
                layersPath, mappingPath, readPath]
    for fileP in fileList:
        if VERBOSE:
            print(f'\nLINTING {fileP.as_posix()}')
        ignoreList = ['E501']
        STR_IGNORE_LIST = "--ignore="+str(str(ignoreList)[1:-1].replace(' ', '').replace("'", ""))
        result = subprocess.run(['flake8', STR_IGNORE_LIST, fileP.as_posix(),],
                                stdout=subprocess.PIPE, check=False)
        if VERBOSE:
            print(result.stdout.decode('utf-8'))

if RUN_TESTS:
    if VERBOSE:
        print('Testing w4h')
    SHELL_TYPE = True
    if sys.platform == 'linux':
        SHELL_TYPE = False
    try:
        subprocess.run(["python", "-m", "pytest", repoDir.as_posix()],
                       shell=SHELL_TYPE, check=False)
    except Exception:
        subprocess.run(["pytest", repoDir.as_posix()], shell=SHELL_TYPE, check=False)


print('docs updated')
