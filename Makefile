#makefile for the tomograph programs
#

CFLAGS =  -traditional
CFLAGS =
LDFLAGS = -L. -lth -lnr -lm
PROGS = rdphase ctomo syphase readsol readsta readpha 
CC = gcc


all: $(PROGS)
	@echo "all made now"
	@echo 

rdphase: data14.h rdphase.c 
	$(CC) $(CFLAGS) -o rdphase rdphase.c $(LDFLAGS)
	cp rdphase ../bin
	@echo "made rdphase"

ctomo: data14.h ctomo.c
	$(CC) $(CFLAGS) -o ctomo ctomo.c $(LDFLAGS)
	cp ctomo ../bin
	@echo "made ctomo"

syphase: data14.h syphase.c
	$(CC) $(CFLAGS) -o syphase syphase.c $(LDFLAGS)
	cp syphase ../bin
	@echo "made syphase"

readsol: data14.h readsol.c
	$(CC) $(CFLAGS) -o readsol readsol.c $(LDFLAGS)
	cp readsol ../bin
	@echo "made readsol"

readsta: data14.h readsta.c
	$(CC) $(CFLAGS) -o readsta readsta.c $(LDFLAGS)
	cp readsta ../bin
	@echo "made readsta"

readpha: data14.h readpha.c
	$(CC) $(CFLAGS) -o readpha readpha.c $(LDFLAGS)
	cp readpha ../bin
	@echo "made readpha"

data14.h:
	@echo "data14.h was updated"

badeve.h:
	@echo "badeve.h was updated"

badsta.h:
	@echo "badsta.h was updated"
	


clean:
	rm $(PROGS)
	@echo "directory cleaned"
