# Makefile for installing package
#

install:
	pip install .

uninstall:
	pip uninstall est

reinstall:
	pip uninstall est
	pip install .

install-dev:
	pip install -e .

