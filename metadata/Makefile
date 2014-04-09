
.PHONY: all

CSVS = $(wildcard plate*.csv)
META = $(CSVS:.csv=_metadata.json)

all: $(META)


%_metadata.json: %_plateinfo.json %.csv
	../venv/bin/make_sample_json.py -o $@ $^
