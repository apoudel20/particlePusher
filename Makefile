SHELL := '/bin/bash'

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

loop5: 
	for frames in 0 1 2 3 4 5 ; do \
		./particlePusher $$frames ; \
	done

loop10: 
	for frames in 0 1 2 3 4 5 6 7 8 9 10 ; do \
		./particlePusher $$frames ; \
	done

loop15: 
	for frames in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ; do \
		./particlePusher $$frames ; \
	done