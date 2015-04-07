CC=g++
FLAGS=-O2 -funroll-loops -fomit-frame-pointer -std=c++11
LDFLAGS=-framework accelerate
OBJECTS= model.o particle.o patch.o tools.o

math_dir="../../../libraries/math"
STATIC_LIBS=-L../../../libs -lmath

clean:
	rm -rf ./*.o ./*.out

main: main.cpp $(OBJECTS)
	$(CC) $(FLAGS) $(LDFLAGS) -o $@.out $@.cpp $(OBJECTS) $(STATIC_LIBS)

.cpp.o: .cpp .h
	$(CC) $(FLAGS) -o $@ -c $<
