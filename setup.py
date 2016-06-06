"""
Setup configuration file for the Enhanced Sampling Toolkit.
"""

from distutils.core import setup

setup(name="est",
      version="0.1a",
      description="A toolkit for rapid prototyping of enhanced sampling algorithms",
      author="Jeremy Tempkin",
      author_email="jtempkin@uchicago.edu",
      url="https://github.com/jtempkin/enhanced-sampling-toolkit",
      requires={},
      package_dir={'applications': 'neus'},
      packages=['neus', 'walker'],
      )
