CXX       = g++
CXXFLAGS  = -Wall
CXXFLAGS += -g
CXXFLAGS += -O2
CXXFLAGS += -fPIC
CXXFLAGS += -I..
CXXFLAGS += -Wno-unused-function
#CXXFLAGS += -std=c++11

LDFLAGS = -lstdc++

CXXFLAGS += `pkg-config --cflags acado`
LDFLAGS  += `pkg-config --libs   acado`

CPP_SRC = rienIntegratorInterface.cpp

.PHONY: clean all rienIntegratorInterface.so
all : $(OBJ) rienIntegratorInterface.so

%.o : %.cpp
	@echo CPP $@ #: $(CXX) $(CXXFLAGS) -c $< -o $@
	@$(CXX) $(CXXFLAGS) -c $< -o $@

OBJ = $(CPP_SRC:%.cpp=%.o)

rienIntegratorInterface.so : $(OBJ)
	@echo LD $@ #: $(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)
	@$(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)

clean :
	rm -f *.o *.so
