YACC	= bison -y
YFLAGS	= -d -v
LEX	= flex

DEP_COMAMND = gcc -MM $(SRCS)

CC = gcc
CFLAGS = -g -Wall -O2
LIBS = -lm -lumfpack -lamd -lblas

TARGET = a.out
SRCS = tc.c sim.c obj.c mem.c solvele.c sparse_matrix.c y.tab.c lex.yy.c
OBJS = tc.o sim.o obj.o mem.o solvele.o sparse_matrix.o y.tab.o lex.yy.o

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS)

y.tab.c y.tab.h: parser.y
	$(YACC) $(YFLAGS) $<

lex.yy.c: lexer.l
	$(LEX) $<

depend: y.tab.c y.tab.h lex.yy.c
	$(DEP_COMAMND) > Makefile.depend

clean:
	rm -f $(OBJS) $(TARGET) y.tab.[hc] y.output lex.yy.c

-include Makefile.depend
