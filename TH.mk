PROGRAM = th1
FILES.c = PHS_thesis.c mincost.c 
FILES.h = PHS.h 
FILES.o = ${FILES.c:.c=.o}
	
CC      = gcc
SFLAGS  = -std=c11
OFLAGS  = -O
CFLAGS  = ${SFLAGS} ${OFLAGS} 
LDFLAGS = -lm

all: ${PROGRAM} 

${PROGRAM}: ${FILES.o}
	${CC} -o $@ ${CFLAGS} ${FILES.o} ${LDFLAGS}

PHS_thesis.o: ${FILES.h}
mincost.o:    ${FILES.h}

SeqID         = B29
LTSize        = 5
id            = b29p2
beta          = 0.8
maxiter       = 10000
iterstride    = 50
stoperr       = 0.0001
epsilon       = 0.0001
seed          = 36
numrun        = 1

run: ${PROGRAM}
	$(addprefix ./,${PROGRAM}) ${SeqID} ${LTSize} ${id} ${beta} ${maxiter} ${iterstride} ${stoperr} ${epsilon} ${seed} ${numrun}

.PHONY: run 
