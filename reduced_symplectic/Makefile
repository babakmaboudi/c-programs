CC=g++
FLAGS=-O2 -funroll-loops -fomit-frame-pointer -std=c++11
#FLAGS=-g -std=c++11
LDFLAGS=-framework accelerate
OBJECTS= model.o parameters.o tools.o grid.o node.o

math_dir="../libraries/math"
STATIC_LIBS=-L../libs -lmath

clean:
	rm -rf ./*.o ./*.out

clean_libs:
	@cd $(math_dir); make clean

build_libs:
	@cd $(math_dir); make build

main: main.cpp $(OBJECTS)
	$(CC) $(FLAGS) $(LDFLAGS) -o $@.out $@.cpp $(OBJECTS) $(STATIC_LIBS)

.cpp.o: .cpp .h
	$(CC) $(FLAGS) -o $@ -c $<
