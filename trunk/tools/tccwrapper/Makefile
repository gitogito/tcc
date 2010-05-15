TARGET = tccwrapper
INCLUDE =
OBJS = syntax.cmo parser.cmo lexer.cmo tccwrapper.cmo
DEPEND = parser.mli parser.ml lexer.ml

all: $(DEPEND) $(TARGET)

$(TARGET): $(OBJS)
	ocamlc $(INCLUDE) -o $(TARGET) $(OBJS)

%.cmi: %.mli
	ocamlc $(INCLUDE) $<

%.cmo: %.ml
	ocamlc -c $(INCLUDE) $<

parser.cmo: parser.cmi

parser.cmi: parser.mli
	ocamlc $(INCLUDE) $<

parser.mli: parser.mly
	ocamlyacc $<

parser.ml: parser.mly
	ocamlyacc $<

lexer.ml: lexer.mll
	ocamllex $<

clean::
	-rm -rf *.cmi *.cmo $(TARGET) *.cmo *.o lexer.ml parser.ml parser.mli

depend:: $(DEPEND)
	ocamldep *.mli *.ml > .depend

-include .depend
