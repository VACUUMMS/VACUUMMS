all: upload conda

pypi: setup.py
	python setup.py sdist

conda: pypi conda-recipe/meta.yaml
	conda-build conda-recipe

upload: pypi
	twine upload dist/*

clean:
	rm -rf build dist vacuumms.egg-info


