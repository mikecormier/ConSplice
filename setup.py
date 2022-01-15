from setuptools import setup
import re

with open("consplice/__init__.py", "r") as fd:
    version = version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',fd.read(), re.MULTILINE).group(1)

requires = ["pyyaml"]

setup(
    name="consplice",
    version=version,
    description="CLI for Constrained Splice (ConSplice) model",
    long_description=open("README.md").read(),
    author="Mike Cormier",
    author_email="cormiermichaelj@gmail.com",
    url="https://github.com/mikecormier/ConSplice",
    packages=["consplice",
              "consplice.constraint",
              "consplice.ML", 
              "consplice.config", 
              "consplice.config.ConSpliceML_Model"],
    package_dir={"consplice":"consplice",
                 "consplice.constraint": "consplice/constraint",
                 "consplice.ML":"consplice/ML"
                 "consplice.config":"consplice/config"
                 "consplice.config.ConSpliceML_Model": "consplice/config/ConSpliceML_Model"},
    package_data={"":["LICENSE","README.md"]},
    include_package_data=True,
    install_requires=requires,
    license="MIT",
    zip_safe=False,

    entry_points={
        "console_scripts": [ "consplice = consplice.__main__:main" ]
    },

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bioinformatics']
)
