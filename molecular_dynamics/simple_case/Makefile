CC=g++
FLAGS=-O2 -funroll-loops -fomit-frame-pointer -std=c++11
OBJECTS= model.o vec.o quat.o matrix.o

clean:
	rm -rf ./*.o ./*.out

main: main.cpp $(OBJECTS)
	$(CC) $(FLAGS) -o $@.out $@.cpp $(OBJECTS) $(STATIC_LIBS)

.cpp.o: .cpp .h
	$(CC) $(FLAGS) -o $@ -c $<
