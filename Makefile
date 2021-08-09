SHELL := '/bin/bash'
FRAME_COUNT=5

all: build run

reset: clean build

build:
	/usr/local/cuda-11.0/bin/nvcc particlePusher.cu -o particlePusher -Xptxas -v 

run:
	./particlePusher

graph:
	python3 graph.py

clean: 
	rm -f particlePusher

rung: run graph


loop: 
	for frames in {0..$(FRAME_COUNT)} ; do \
		./particlePusher $$frames ; \
	done

profile:
	nvprof --metrics all --log-file ./log_profiles/log-100-$$(date +%h_%m_%s) ./particlePusher 2
