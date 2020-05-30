SLATEDIR = $(HOME)/files/codes/inq/slate
CUDADIR = /usr/local/cuda
CXX = mpic++
CXXFLAGS = -Ofast -fopenmp \
	-Wall -Wextra -Wno-unused-parameter -Wno-cast-function-type \
	-I$(SLATEDIR)/include          \
	-I$(SLATEDIR)/blaspp/include   \
	-I$(SLATEDIR)/lapackpp/include \
	-I$(CUDADIR)/include

minq: minq.cpp

clean:
	rm -f minq
