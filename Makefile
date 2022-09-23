COMMIT_VERSION=$(shell git rev-parse HEAD | head -c 7)
COMMIT_DIFF_SH=$(shell git diff HEAD --shortstat)
COMMIT_DIFF_FULL=$(shell echo "R\"ZXF_specQUOTE(\n $$(git diff HEAD | sed 's/ZXF_specQUOTE/ZXF_specquote/g') \n)ZXF_specQUOTE\"" > gitdiff.txt)

CXX=g++
# For some systems, libhts.so is not available or cannot be used. Therefore, static linking to htslib is used here. 
CXXFLAGS=-I ext/htslib-1.11-lowdep/ -pthread -lm -lz -lbz2 -llzma -static-libstdc++ -Bstatic -lhts
LDFLAGS=-L ext/htslib-1.11-lowdep/ 
VERFLAGS=-DCOMMIT_VERSION="\"$(COMMIT_VERSION)\"" -DCOMMIT_DIFF_SH="\"$(COMMIT_DIFF_SH)\"" -DCOMMIT_DIFF_FULL="\"$(COMMIT_DIFF_FULL)\""

all: safemut safemut.debug safemix safemix.debug
	
safemut : safemut.cpp Makefile
	$(CXX) -o safemut -O2 safemut.cpp $(CXXFLAGS) $(LDFLAGS) $(VERFLAGS)
safemut.debug : safemut.cpp Makefile
	$(CXX) -o safemut.debug -O0 -g -p -fsanitize=address safemut.cpp $(CXXFLAGS) $(LDFLAGS) $(VERFLAGS)
safemix : safemix.cpp Makefile
	$(CXX) -o safemix -O2 safemix.cpp $(CXXFLAGS) $(LDFLAGS) $(VERFLAGS)
safemix.debug : safemix.cpp Makefile
	$(CXX) -o safemix.debug -O0 -g -p -fsanitize=address safemix.cpp $(CXXFLAGS) $(LDFLAGS) $(VERFLAGS)

.PHONY: clean deploy
	
clean:
	rm safemut safemut.debug safemix safemix.debug bin/safemut bin/safemix
deploy:
	cp safemut safemix bin/

