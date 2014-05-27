from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "python-fit",
    version = "1.0",
    author = "Michael Woods",
    author_email = "physicsmichael@gmail.com",
    description = ("A python module using scipy's orthogonal distance regression that makes fitting data easy."),
    license = "MIT",
    keywords = "fitting curve",
    url = "https://github.com/vgm64/python-fit",
    packages=['fit', 'tests'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering", 
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License"
    ],
)
