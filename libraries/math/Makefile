CC=g++
FLAGS = -O2 -funroll-loops -fomit-frame-pointer -std=c++11

OBJECTS=matrix.o
OUT=../../libs/libmath.a

build: $(OBJECTS)
	ar rcs $(OUT) $(OBJECTS)

clean:
	rm -rf *.o $(OUT)

.cpp.o:
	$(CC) $(FLAGS) -o $@ $< -c
