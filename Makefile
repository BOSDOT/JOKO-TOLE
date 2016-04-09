CC          = g++
CFLAGS      = -Wall -ansi -pedantic -msse3 -std=c++0x -g -O3 
LDFLAGS     = -lpthread -flto
CXFLAGS    =  -pedantic -Wextra -Wshadow
OBJS        = board.o general.o zobrist.o hash.o search.o mover.o
ENGINENAME  = Joko

ifeq ($(USE_STATIC), true)
	LDFLAGS += -static -static-libgcc -static-libstdc++
endif

all: uci

uci: $(OBJS) uci.o
	$(CC) -o $(ENGINENAME)$(EXT) $^ $(LDFLAGS) $(CXFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS) -x c++ $< -o $@

clean:
	rm -f *.o $(ENGINENAME)$(EXT).exe $(ENGINENAME)$(EXT)

strip: strip joko.exe
