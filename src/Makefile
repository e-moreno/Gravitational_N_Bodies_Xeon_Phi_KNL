CC=icc
OMPFLAG=-qopenmp
#FLOATFLAG=-DDATATYPE=float -DPOW=powf -DSQRT=sqrtf
BLOCKFLAG=-DBLOCKSIZE=16
AVX256=-xAVX2
AVX512=-xMIC-AVX512
FPMODEL=-fp-model fast=2
DOUBLEFLAGS = -DDATATYPE=double -DPOW=pow -DSQRT=sqrt -DBLOCKSIZE=8

default: opt0.exe
opt0: opt0.exe
opt1:	opt1.exe
opt2:	opt2.exe
opt3:	opt3.exe
opt4:	opt4.exe
opt5:	opt5.exe
opt6:	opt6.exe
opt7:	opt7.exe
opt8:	opt8.exe
opt9: opt9.exe
opt10:	opt10.exe
opt11:	opt11.exe
opt12:	opt12.exe
opt13:	opt13.exe
move:
	mkdir bin
	mv *.exe bin
all: opt0 opt1 opt2 opt3 opt4 opt5 opt6 opt7 opt8 opt9 opt10 opt11 opt12 opt13 move

opt0.exe: nbodies_opt0.c
	$(CC) -o opt0.exe nbodies_opt0.c -lm

opt1.exe: nbodies_opt1.c
	$(CC) -o opt1.exe nbodies_opt1.c -lm $(OMPFLAG)

opt2.exe: nbodies_opt2.c
	$(CC) -o opt2.exe nbodies_opt2.c -lm $(OMPFLAG)

opt3.exe: nbodies_opt3.c
	$(CC) -o opt3.exe nbodies_opt3.c -lm $(OMPFLAG)

opt4.exe: nbodies_opt4.c
	$(CC) -o opt4.exe nbodies_opt4.c -lm $(OMPFLAG)

opt5.exe: nbodies_opt5_6_7.c
	$(CC) -o opt5.exe nbodies_opt5_6_7.c -lm $(OMPFLAG)

opt6.exe: nbodies_opt5_6_7.c
	$(CC) -o opt6.exe nbodies_opt5_6_7.c -lm $(OMPFLAG) $(AVX256)

opt7.exe: nbodies_opt5_6_7.c
	$(CC) -o opt7.exe nbodies_opt5_6_7.c -lm $(OMPFLAG) $(AVX512)

opt8.exe: nbodies_opt8.c
	$(CC) -o opt8.exe nbodies_opt8.c -lm $(OMPFLAG) $(AVX512)

opt9.exe: nbodies_opt9.c
	$(CC) -o opt9.exe nbodies_opt9.c -lm $(OMPFLAG) $(AVX512) $(BLOCKFLAG)

opt10.exe: nbodies_opt10_11.c
	$(CC) -o opt10.exe nbodies_opt10_11.c -lm $(OMPFLAG) $(AVX512) $(BLOCKFLAG)

opt11.exe: nbodies_opt10_11.c
	$(CC) -o opt11.exe nbodies_opt10_11.c -lm $(OMPFLAG) $(AVX512) $(BLOCKFLAG) $(FPMODEL)

opt12.exe: nbodies_opt12_13.c
	$(CC) -o opt12.exe nbodies_opt12_13.c -lm $(OMPFLAG) $(AVX512) $(DOUBLEFLAGS)

opt13.exe: nbodies_opt12_13.c
	$(CC) -o opt13.exe nbodies_opt12_13.c -lm $(OMPFLAG) $(AVX512) $(DOUBLEFLAGS) $(FPMODEL)

clean:
	-rm -f *.o
	-rm -f *.exe
	-rm -rf bin
