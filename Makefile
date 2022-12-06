IDIR =include
CC=gcc
CFLAGS=-I$(IDIR)

ODIR=obj
LDIR=lib
SDIR=src

LIBS=-lm

_DEPS = pso_b.h pso_bi.h pso_c.h pso_ci.h pso_d.h pso_di.h pso_m.h pso_mi.h pso_n.h pso_ni.h  
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = PSO_v5.o pso_b.o pso_bi.o pso_c.o pso_ci.o pso_d.o pso_di.o pso_m.o pso_mi.o pso_n.o pso_ni.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

PSO_v5: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 