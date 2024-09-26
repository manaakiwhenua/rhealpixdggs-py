#!/usr/bin/env python

from setuptools import setup, find_packages
import codecs
import os

VERSION = "0.5.4"


def open_local(paths, mode="r", encoding="utf8"):
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), *paths)
    return codecs.open(path, mode, encoding)


with open_local(["README.rst"], encoding="utf-8") as readme:
    long_description = readme.read()


with open_local(["requirements.txt"]) as req:
    install_requires = req.read().split("\n")


setup(
    name="rHEALPixDGGS",
    version=VERSION,
    author="Alexander Raichev",
    author_email="alex@raichev.net",
    packages=["rhealpixdggs"],
    url="https://github.com/manaakiwhenua/rhealpixdggs-py",
    download_url="https://github.com/manaakiwhenua/rhealpixdggs-py/archive/v{:s}.tar.gz".format(
        VERSION
    ),
    license="LICENSE.txt",
    description="An implementation of the rHEALPix discrete global grid system",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    keywords=[
        "DGGS",
        "Discrete Global Grid System",
        "reference system",
        "spatial",
        "geospatial",
    ],
    project_urls={
        "Bug Reports": "https://github.com/manaakiwhenua/rhealpixdggs-py/issues",
        "Source": "https://github.com/manaakiwhenua/rhealpixdggs-py",
    },
    install_requires=install_requires,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
