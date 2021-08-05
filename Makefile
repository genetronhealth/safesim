COMMIT_VERSION=$(shell git rev-parse HEAD | head -c 7)
COMMIT_DIFF_SH=$(shell git diff HEAD --shortstat)
COMMIT_DIFF_FULL=$(shell echo "R\"ZXF_specQUOTE(\n $$(git diff HEAD | sed 's/ZXF_specQUOTE/ZXF_specquote/g') \n)ZXF_specQUOTE\"" > gitdiff.txt)
VERFLAGS=-DCOMMIT_VERSION="\"$(COMMIT_VERSION)\"" -DCOMMIT_DIFF_SH="\"$(COMMIT_DIFF_SH)\"" -DCOMMIT_DIFF_FULL="\"$(COMMIT_DIFF_FULL)\""
CXXFLAGS=-static-libstdc++ umi-mut-sim.cpp ext/htslib-1.11-lowdep/libhts.a -I ext/htslib-1.11-lowdep/ -pthread -lm -lz -lbz2 -llzma

all: umi-mut-sim.out umi-mut-sim.debug.out
	
umi-mut-sim.out : umi-mut-sim.cpp Makefile
	g++ -o umi-mut-sim.out -O2 $(CXXFLAGS) $(VERFLAGS)
umi-mut-sim.debug.out : umi-mut-sim.cpp Makefile
	g++ -o umi-mut-sim.debug.out -O0 -g -p -fsanitize=address $(CXXFLAGS) $(VERFLAGS)
	
.PHONY: clean deploy
	
clean:
	rm *.out
deploy:
	cp umi-mut-sim.out bin/

