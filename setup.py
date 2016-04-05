"""Standard Python setup script for GPkit Models"""
from __future__ import print_function
import sys
from distutils.core import setup

LONG_DESCRIPTION = """
GPkit Models is a library of geometric programming and signomial programming
models that can be manipulated and solved using 
`GPkit <https://github.com/hoburg/gpkit/>`_.

`Documentation <http://gpkit.rtfd.org/>`_

`Citing GPkit <http://gpkit.rtfd.org/en/latest/citinggpkit.html>`_
"""

LICENSE = """The MIT License (MIT)

Copyright (c) 2015 MIT Hoburg Research Group

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE."""

setup(
    name="gpkitmodels",
    description="Library of geometric and signomial programming models "
                "that can be manipulated and solved using GPkit.",
    author="MIT Department of Aeronautics and Astronautics",
    author_email="gpkit@mit.edu",
    url="https://www.github.com/hoburg/gpkit-models",
    install_requires=["numpy", "scipy", "pint"],
    version="0.0.0.0",
    packages=["gpkit-models"],
    package_data={},
    license=LICENSE,
    long_description=LONG_DESCRIPTION,
)
