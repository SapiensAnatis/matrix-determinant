# PHYM0004 Assignment 1 / Jay Malhotra / 10/10/2021
# The makefile is used to build the project. The 'main' target builds the
# primary executable of the program which gathers performance data, whereas the
# 'test' target uses the main() defined in test.c instead to run some unit
# tests, and generates a binary of a different name. The 'test' target requires
# cmocka to be installed on the system.

SOURCES = matrix_core.c determinant.c
HEADERS = matrix_core.c determinant.c

default: main

main: main.c $(SOURCES) $(HEADERS)
	gcc main.c $(SOURCES) -o matrix -Wall -Ofast

test: test.c $(SOURCES) $(HEADERS)
	gcc test.c $(SOURCES) -o testmatrix -Wall -g -lcmocka

clean:
	rm -f *.o
	rm -f matrix
	rm -f testmatrix