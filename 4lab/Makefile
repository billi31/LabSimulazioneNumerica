SRCDIR := src
HEADIR := include
OBJDIR := obj
BINDIR := bin

CXX = clang++ -g
CXXFLAGS = -Wall -O3 --std=c++17 -I./include

EXE := $(BINDIR)/main
SRC := $(wildcard $(SRCDIR)/*.cpp)
HEA := $(wildcard $(HEADIR)/*.hpp)
OBJ := $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)


./PHONY : all clean

all : $(EXE)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEA)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(EXE): $(OBJ)
	$(CXX) -o $(EXE) $(OBJ)

clean :
	-rm $(OBJ) $(EXE)
