INCLUDES= -I../FangOost -I../FunctionalUtilities -I../GaussNewton -I../TupleUtilities -I../AutoDiff -I../CharacteristicFunctions -I../RungeKutta
GCCVAL=g++

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	GCCVAL=g++-7
endif
test:test.o 
	$(GCCVAL) -std=c++14 -O3 -pthread --coverage -g  test.o $(INCLUDES) -o test -fopenmp

test.o: test.cpp CFDistUtilities.h 
	$(GCCVAL) -std=c++14 -O3 -pthread --coverage -D VERBOSE_FLAG=1  -g  -c test.cpp   $(INCLUDES) -fopenmp

clean:
	-rm *.o test *.out



