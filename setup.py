from setuptools import setup

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="wells4hydrogeology",
    author= "Riley Balikian",
    author_email = "balikian@illinois.edu",
    version="0.0.14",
    install_requires=["geopandas", "rioxarray", "owslib", "scipy", "matplotlib", "pandas", "numpy"],
    long_description = long_description,
    long_description_content_type="text/markdown",
    description="A package to read in geology data from wells and create a layered, gridded hydrogeologic model of a study region, all within a python environment, automating and performing tasks often carried out in a dedicated GIS software.",
    package_data={'w4h': ['resources/*', 'resources/sample_data/*', 'resources/sample_data/DictionaryTerms/*','resources/sample_data/LithologyInterpretations/*']}
    )