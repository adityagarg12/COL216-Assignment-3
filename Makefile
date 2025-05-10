CC = g++
CFLAGS = -std=c++11 -Wall -O2 -g
all: L1simulate
L1simulate: L1_cache_sim.cpp
	$(CC) $(CFLAGS) -o L1simulate L1_cache_sim.cpp
report:
	@pdflatex report.tex
clean:
	rm -f L1simulate *.o report.pdf report.out report.log report.aux