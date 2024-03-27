from setuptools import setup, find_packages

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="wells4hydrogeology",
    author= "Riley Balikian",
    author_email = "balikian@illinois.edu",
    version="0.0.22-dev",
    package_data={'w4h': ['resources/*', 'resources/sample_data/*', 'resources/sample_data/statewide_sample_data/*',
                          'resources/sample_data/DictionaryTerms/*','resources/sample_data/LithologyInterpretations/*']},
    install_requires=["geopandas", "rioxarray", "owslib", "scipy", "matplotlib", "pandas", "numpy"],
    long_description_content_type="text/markdown",
    long_description = long_description,
    packages=find_packages(),
    )