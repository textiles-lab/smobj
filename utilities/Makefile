.PHONY : all

#From: https://stackoverflow.com/questions/714100/os-detecting-makefile
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CPP=g++ -std=c++17 -Wall -Werror -g
else
	CPP=clang++ -std=c++17 -Wall -Werror -g
endif

all : knitout-to-smobj smobj-to-yarns

knitout-to-smobj : knitout-to-smobj.cpp
	$(CPP) -o '$@' '$<'

smobj-to-yarns : smobj-to-yarns.cpp sm.hpp sm.o
	$(CPP) -o '$@' '$<' sm.o

sm.o : sm.cpp sm.hpp
	$(CPP) -c -o '$@' '$<'
