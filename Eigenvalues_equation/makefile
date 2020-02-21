all: test dotest
test: main.o eigen.o
	g++ -o prog main.o eigen.o --std=c++11 
main.o: main.cpp eigen.h
	g++ -c  main.cpp --std=c++11

eigen.o: eigen.cpp eigen.h
	g++ -c  eigen.cpp --std=c++11
	
dotest: test
	./prog
