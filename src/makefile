# Makefile for scaff10x 
CC= gcc
CFLAGS= -O2
LFLAGS= -lm

SOURCES=scanPAV_seqout.c scanPAV_fasta.c scanPAV_shred.c scanPAV_confirm.c scanPAV_gapsize.c scanPAV.c scanPAV_bwasam.c scanPAV_bwatag.c scanPAV_bwagap.c scanPAV_clean.c scanPAV_cover.c scanPAV_process.c  

#SOURCES1=scanPAV_seqout.c scanPAV_fasta.c scanPAV_shred.c scanPAV_confirm.c scanPAV_gapsize.c #these are only codes that need fast.o

OBJS = $(patsubst %.c,%.o,$(SOURCES)) fast.o
EXECS = $(patsubst %.c,%,$(SOURCES))
EXECS_BIN = $(patsubst %.c,scanPAV-bin/%,$(SOURCES))
COMPILE = $(CC) $(CFLAGS)


all:  cleanall iprint $(OBJS) executables clean oprint

executables:
	for exe in $(EXECS);  do $(COMPILE) -o $$exe $$exe.o fast.o; cp $$exe scanPAV-bin/.; done

%.o: %.c fasta.h
	$(CC) $(CFLAGS)  -c $< -o $@

iprint:
	@echo '+++ Compiling All ... '

oprint:
	@echo 'All Done '


clean: 
	@echo '+++ Cleaning Up ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS)
	@cp scanPAV-bin/scanPAV .

cleanall: 
	@echo '+++ Cleaning All ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS) fast.o
	@rm -f $(EXECS_BIN)
