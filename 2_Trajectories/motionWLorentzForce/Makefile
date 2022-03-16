CC=gcc-11
CFLAGS= -Wall -O2 -Wdouble-promotion
DEPS = 
INCDIR=-I. -I/usr/local/opt/lapack/include -I/usr/local/Cellar/openblas/0.3.20/include
LIBDIR=-L/usr/local/opt/lapack/lib -L/usr/local/Cellar/openblas/0.3.20/lib 
LIBS= -llapacke -lm -lopenblas -lpthread -lgfortran
OBJ = motionWLorentzForce.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) -c -o $@ $< 

motionWLorentzForce: $(OBJ)
	$(CC) $(CFLAGS) $(INCDIR) $(LIBDIR) $(LIBS) -o $@ $^

run:
	./motionWLorentzForce 1>logTime.tsv 2>&1 

plotsOnScreen:
	tail -500 logTime > file2Plot.dat
	gnuplot -p plotTrajectory.plt
	gnuplot -p plotPositions.plt

plotMatplotLib:
	./plotTrajectory.py -i logTime

cleanSim:
	rm -f logTime file2Plot.dat

clean:
	rm -f *.o motionWLorentzForce logTime file2Plot.dat

