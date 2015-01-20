CC=g++
FLAGS=-O2 -funroll-loops -fomit-frame-pointer -std=c++11
LDFLAGS=-framework accelerate
OBJECTS=orbit.o polynomial.o basis.o vec.o

clean:
	rm -rf ./*.o ./*.out

main: main.cpp $(OBJECTS)
	$(CC) $(FLAGS) $(LDFLAGS) -o $@.out $@.cpp $(OBJECTS)

.cpp.o: .cpp .h
	$(CC) $(FLAGS) -o $@ -c $<
