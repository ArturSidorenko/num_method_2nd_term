all: test dotest
test: main.o diffur.o
	g++ -o prog main.o diffur.o --std=c++11 
main.o: main.cpp diffur.h
	g++ -c  main.cpp --std=c++11

diffur.o: diffur.cpp diffur.h
	g++ -c  diffur.cpp --std=c++11
	
dotest: test
	./prog
