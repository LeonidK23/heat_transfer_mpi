all:
	mpicxx -O3 -o stencil main.cpp src/transfer.cpp
