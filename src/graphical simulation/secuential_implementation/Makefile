default: secuential_nbodies.exe
final: nbodies_final.exe

nbodies_secuencial_v0.o: nbodies_secuencial_v0.c
	gcc -c nbodies_secuencial_v0.c -o nbodies_secuencial_v0.o

secuential_nbodies.exe: nbodies_secuencial_v0.o
	gcc nbodies_secuencial_v0.o -o secuential_nbodies.exe -lm

nbodies_final.o: nbodies_final.c
	gcc -c nbodies_final.c -o nbodies_final.o

nbodies_final.exe: nbodies_final.o
	gcc nbodies_final.o -o nbodies_final.exe -lm

clean:
	-rm -f *.o
	-rm -f *.exe
