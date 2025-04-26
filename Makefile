CC = g++
CFLAGS = -std=c++11 -Wall -O2

all: L1simulate

L1simulate: L1_cache_sim.cpp
	$(CC) $(CFLAGS) -o L1simulate L1_cache_sim.cpp

clean:
	rm -f L1simulate *.o