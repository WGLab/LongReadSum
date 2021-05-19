all:
	swig -c++ -python lrst.i
	python setup.py build_ext --inplace

