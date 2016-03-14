CC     = g++
CFLAGS = 
EFLAGS =  
EFILE  = Simera
LIBS   = -lm -lgsl -lgslcblas
OBJS   = Simera.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : Simera.h
	$(CC) $(CFLAGS) -c -w $*.c 

clean:
	rm -rf *.o Simera