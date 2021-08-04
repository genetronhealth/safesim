g++ -o umi-mut-sim.out -O2 -static-libstdc++ umi-mut-sim.cpp ../uvc/ext/htslib-1.11-lowdep/libhts.a -I ../uvc/ext/htslib-1.11-lowdep/ -pthread -lm -lz -lbz2 -llzma
# -fsanitize=address 
