CC = g++
CFLAGS = -Wall -O3 --std=c++11

SRCDIR := src
HEADIR := include
OBJDIR := obj
BINDIR := bin

main1 : $(OBJDIR)/main1.o $(OBJDIR)/lib.o $(OBJDIR)/random.o
	$(CC) $(OBJDIR)/lib.o $(OBJDIR)/random.o $(OBJDIR)/main1.o -o $(BINDIR)/main1
$(OBJDIR)/main1.o : $(SRCDIR)/main1.cpp
	$(CC) -c $(SRCDIR)/main1.cpp -o $(OBJDIR)/main1.o $(CFLAGS)


main2 : $(OBJDIR)/main2.o $(OBJDIR)/lib.o $(OBJDIR)/random.o
	$(CC) $(OBJDIR)/lib.o $(OBJDIR)/random.o $(OBJDIR)/main2.o -o $(BINDIR)/main2
$(OBJDIR)/main2.o : $(SRCDIR)/main2.cpp
	$(CC) -c $(SRCDIR)/main2.cpp -o $(OBJDIR)/main2.o $(CFLAGS)


main3 : $(OBJDIR)/main3.o $(OBJDIR)/lib.o $(OBJDIR)/random.o
	$(CC) $(OBJDIR)/lib.o $(OBJDIR)/random.o $(OBJDIR)/main3.o -o $(BINDIR)/main3
$(OBJDIR)/main3.o : $(SRCDIR)/main3.cpp
	$(CC) -c $(SRCDIR)/main3.cpp -o $(OBJDIR)/main3.o $(CFLAGS)

$(OBJDIR)/lib.o : $(SRCDIR)/lib.cpp $(HEADIR)/lib.hpp
	$(CC) -c $(SRCDIR)/lib.cpp -o $(OBJDIR)/lib.o $(CFLAGS)
$(OBJDIR)/random.o : $(SRCDIR)/random.cpp $(HEADIR)/random.hpp
	$(CC) -c $(SRCDIR)/random.cpp -o $(OBJDIR)/random.o $(CFLAGS)

clean :
	rm $(OBJDIR)/*.o $(BINDIR)/*

