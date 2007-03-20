CC = gcc
CFLAGS = -g -Wall -O0 -I/usr/include/ufsparse
LIBS = -lgc -lm -lumfpack -lamd

TARGET = a.out
SRCS = tc.c sim.c mem.c solvele.c sparse_matrix.c
OBJS = tc.o sim.o mem.o solvele.o sparse_matrix.o

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS)

depend:
	makedepend $(SRCS)

clean:
	rm -f $(OBJS) $(TARGET)

