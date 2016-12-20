OBJS = distkm.o getarg.o errprint.o math.o vax2sun.o
SRCS = distkm.c getarg.c errprint.c math.c vax2sun.c
PROGS = distkm getarg errprint math vax2sun

all: $(OBJS)
	@echo "all made now"
	@echo

$(OBJS): $(SRCS)
	gcc -c $(SRCS)
	ar rcv libth.a $(OBJS)
	ranlib libth.a
	@echo $(SRCS) built
