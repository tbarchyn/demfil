# Makefile for compiling in Linux

make: main.cpp
	g++ main.cpp -Wall -pedantic -fopenmp -no-stack-protector -O3 -o filter.exe
	