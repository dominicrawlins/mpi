stencil: stencil.c
	icc -std=c99 -Wall $^ -o $@

