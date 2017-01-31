from setuptools import setup, find_packages

setup(name="est",
      version="0.1a1",

      #package_dir={'applications': 'neus'},
      packages=find_packages(),
      install_requires=[
        'numpy'
      ],

      # metadata
      author="Jeremy O. B. Tempkin",
      author_email="jtempkin@uchicago.edu",
      url="https://github.com/jtempkin/enhanced_sampling_toolkit",
      description="A toolkit for rapid prototyping of enhanced sampling algorithms.",
      license="MIT",
      keywords="enhanced sampling toolkit molecular dynamics"
      # should add some classifiers
      # also setup a unit testing hook into this
      )
