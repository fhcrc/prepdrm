

test:
	R --slave --vanilla -e "library(devtools); library(testthat); load_all('.'); test_package('prepdrm')"

doc:
	R --slave --vanilla -e "library(devtools); document('.');"

check:
	R --slave --vanilla -e "library(devtools); check('.');"

.PHONY: test doc
