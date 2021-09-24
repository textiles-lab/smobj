.PHONY : all

#From: https://stackoverflow.com/questions/714100/os-detecting-makefile
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CPP=g++ -std=c++17 -Wall -Werror -g -DGLM_ENABLE_EXPERIMENTAL -O2
else
	CPP=clang++ -std=c++17 -Wall -Werror -g -DGLM_ENABLE_EXPERIMENTAL -O2
endif

all : knitout-to-smobj smobj-to-yarns text-to-smobj simplify-yarns merge-faces

clean :
	rm -f knitout-to-smobj text-to-smobj smobj-to-yarns simplify-yarns merge-faces sm.o

knitout-to-smobj : knitout-to-smobj.cpp
	$(CPP) -o '$@' '$<'


text-to-smobj : text-to-smobj.cpp sm.hpp sm.o
	$(CPP) -o '$@' '$<' sm.o

smobj-to-yarns : smobj-to-yarns.cpp sm.hpp sm.o
	$(CPP) -o '$@' '$<' sm.o

simplify-yarns : simplify-yarns.cpp sm.hpp sm.o
	$(CPP) -o '$@' '$<' sm.o

sm.o : sm.cpp sm.hpp
	$(CPP) -c -o '$@' '$<'

hinters.o : hinters.cpp hinters.hpp
	$(CPP) -c -o '$@' '$<'

merge-faces : merge-faces.cpp sm.hpp sm.o
	$(CPP) -o '$@' '$<' sm.o

test-code : test-code.cpp hinters.hpp sm.hpp sm.o hinters.o
	$(CPP) -o '$@' '$<'  hinters.o sm.o

#generic rules for knitout -> smobj -> yarns

%.smobj : %.knitout knitout-to-smobj
	./knitout-to-smobj '$<' '$@'

%.yarns : %.smobj ../faces/knitout.sf smobj-to-yarns
	./smobj-to-yarns '$<' '../faces/knitout.sf' '$@'
