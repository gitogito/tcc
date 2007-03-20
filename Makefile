CC = gcc
CFLAGS = -g -Wall -O0
LIBS = -lgc

TARGET = a.out
SRCS = tc.c sim.c mem.c
OBJS = tc.o sim.o mem.o

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS)

depend:
	makedepend $(SRCS)

clean:
	rm -f $(OBJS) $(TARGET)

