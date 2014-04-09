.PHONY : clean inplace build test

PYTHON ?= python

EXTS = prep_drm/ssw.so prep_drm/rle.so prep_drm/pwalign.so prep_drm/samutil.so
#CYTHON_COMPILED  = prep_drm/cstrutils.c prep_drm/codons.cpp prep_drm/ssw.cpp

build:
	$(PYTHON) setup.py build

clean:
	$(PYTHON) setup.py clean
	rm -rf build dist
	#rm -f $(EXTS) $(CYTHON_COMPILED)

inplace: $(EXTS)

test:
	$(PYTHON) setup.py test

%.so: %.pyx
	$(PYTHON) setup.py build_ext --inplace
