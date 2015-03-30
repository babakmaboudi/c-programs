CC=g++
FLAGS=-O2 -funroll-loops -fomit-frame-pointer -std=c++11
OBJECTS= common.o

clean:
	rm -rf ./*.o ./*.out

main: main.cpp $(OBJECTS)
	$(CC) $(FLAGS) -o $@.out $@.cpp $(OBJECTS)

.cpp.o: .cpp .h
	$(CC) $(FLAGS) -o $@ -c $<
