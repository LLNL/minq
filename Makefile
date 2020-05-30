SLATEDIR = $(HOME)/files/codes/inq/slate
CUDADIR = /usr/local/cuda
CXX = mpic++
CXXFLAGS = -g -Ofast -fopenmp \
	-Wall -Wextra -Wno-unused-parameter -Wno-cast-function-type \
	-I$(SLATEDIR)/include          \
	-I$(SLATEDIR)/blaspp/include   \
	-I$(SLATEDIR)/lapackpp/include \
	-I$(CUDADIR)/include

LDFLAGS = -L$(SLATEDIR)/lib -L$(CUDADIR)/lib64/
LIBS = -lslate -lcublas -lcudart

minq: minq.cpp
	$(CXX) $(CXXFLAGS) minq.cpp -o minq $(LDFLAGS) $(LIBS)

clean:
	rm -f minq

