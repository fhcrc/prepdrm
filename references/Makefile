# Makefile for preparing reference sequences
# Not part of main SConstruct - this should only be done once

.PHONY: all

PYTHON=../venv/bin/python

BASE = prep_refs_MB2059

all: $(BASE).fasta $(BASE).bounds.json

$(BASE).fasta: prep_refs.fasta MB2059_genotypes.fasta ../bin/choose_reference_per_subject.py
	$(PYTHON) $(word 3, $^) < $< | muscle -profile -in2 - -in1 $(word 2,$^) | seqmagick convert --squeeze - $@

$(BASE).bounds.json: ../bin/reference_bounds.py prep_refs_MB2059.fasta
	$(PYTHON) $^ -o $@
