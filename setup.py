from setuptools import setup

setup(
    name="wells4hydrogeology",
    author= "Riley Balikian",
    author_email = "balikian@illinois.edu",
    version="0.0.1",
    install_requires=["geopandas", "rioxarray", "owslib", "scipy", "matplotlib", "pandas", "numpy"],
    description="A package to read in geology data from wells and create a layered, gridded hydrogeologic model of a study region, all within a python environment, automating and performing tasks often carried out in a dedicated GIS software."
    )
