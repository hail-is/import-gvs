TMPDIR := build/tmp/

current_dir = $(shell pwd)

PYTHON3 := PYTHONPATH="$(current_dir):$(PYTHONPATH)" python3

clean:
	rm -rf build/*

build/avro/.sentinel: test/resources/avro.tar.gz
	mkdir -p build/tmp/
	gcloud storage cp gs://hail-common/gvs-resources/avro.tar.gz build/
	tar -xzvf build/avro.tar.gz -C build
	touch $@

build/gvs-vcfs/.sentinel: test/resources/gvs-vcfs.tar.gz
	mkdir -p build/tmp/
	gcloud storage cp gs://hail-common/gvs-resources/gvs-vcfs.tar.gz build/
	tar -xzvf test/resources/gvs-vcfs.tar.gz -C build
	touch $@

build/out.vds/.sentinel: build/avro/.sentinel build/gvs-vcfs/.sentinel
	$(PYTHON3) scripts/run_import.py build/out.vds $(TMPDIR)
	touch $@

build/dense.mt/.sentinel: build/out.vds/.sentinel
	$(PYTHON3) scripts/make_dense.py build/out.vds build/dense.mt $(TMPDIR)
	touch $@

build/joined.mt/.sentinel: build/dense.mt/.sentinel
	$(PYTHON3) scripts/join_with_vcfs.py build/dense.mt build/joined.mt $(TMPDIR)
	touch $@

clean-temp:
	rm -rf $(TMPDIR)

build-resources: build/joined.mt/.sentinel clean-temp
