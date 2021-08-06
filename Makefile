COMMIT_VERSION=$(shell git rev-parse HEAD | head -c 7)
COMMIT_DIFF_SH=$(shell git diff HEAD --shortstat)
COMMIT_DIFF_FULL=$(shell echo "R\"ZXF_specQUOTE(\n $$(git diff HEAD | sed 's/ZXF_specQUOTE/ZXF_specquote/g') \n)ZXF_specQUOTE\"" > gitdiff.txt)
VERFLAGS=-DCOMMIT_VERSION="\"$(COMMIT_VERSION)\"" -DCOMMIT_DIFF_SH="\"$(COMMIT_DIFF_SH)\"" -DCOMMIT_DIFF_FULL="\"$(COMMIT_DIFF_FULL)\""
CXXFLAGS=-static-libstdc++ safesim.cpp ext/htslib-1.11-lowdep/libhts.a -I ext/htslib-1.11-lowdep/ -pthread -lm -lz -lbz2 -llzma

all: safesim safesim.debug
	
safesim : safesim.cpp Makefile
	g++ -o safesim -O2 $(CXXFLAGS) $(VERFLAGS)
safesim.debug : safesim.cpp Makefile
	g++ -o safesim.debug -O0 -g -p -fsanitize=address $(CXXFLAGS) $(VERFLAGS)
	
.PHONY: clean deploy
	
clean:
	rm safesim safesim.debug bin/safesim
deploy:
	cp safesim bin/

