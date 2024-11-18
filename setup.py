"""Setup project"""

from pathlib import Path
from setuptools import setup, find_packages

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
        name="wells4hydrogeology",
        author="Riley Balikian",
        author_email="balikian@illinois.edu",
        version="0.0.25",
        package_data={'w4h':
                      ['resources/*', 'resources/sample_data/*',
                       'resources/sample_data/statewide_sample_data/*',
                       'resources/sample_data/DictionaryTerms/*',
                       'resources/sample_data/LithologyInterpretations/*']},
        install_requires=["geopandas", "rioxarray", "owslib", "scipy",
                          "matplotlib", "pandas", "numpy"],
        long_description_content_type="text/markdown",
        long_description=long_description,
        packages=find_packages(),
        )
