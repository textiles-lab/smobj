
#From: https://stackoverflow.com/questions/714100/os-detecting-makefile
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CPP=g++ -std=c++17 -Wall -Werror -g
else
	CPP=clang++ -std=c++17 -Wall -Werror -g
endif

knitout-to-smobj : knitout-to-smobj.cpp
	$(CPP) -o '$@' '$<'
