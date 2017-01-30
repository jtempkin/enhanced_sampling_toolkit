# Makefile for installing package
#

install:
	pip install .

uninstall:
	pip uninstall jester

reinstall:
	pip uninstall jester
	pip install .

install-dev:
	pip install -e .

test: 
	py.test tests
