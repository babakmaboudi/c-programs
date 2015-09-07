CC=g++
FLAGS=-O2 -funroll-loops -fomit-frame-pointer -std=c++11
LDFLAGS=-framework accelerate
OBJECTS= model.o vec.o quat.o matrix.o

VTK_INCLUDE=-I/usr/local/Cellar/vtk/6.2.0/include/vtk-6.2
VTK_LIBS=-L/usr/local/Cellar/vtk/6.2.0/lib

clean:
	rm -rf ./*.o ./*.out

main: main.cpp $(OBJECTS)
	$(CC) $(VTK_INCLUDE) $(FLAGS) $(LDFLAGS) -o $@.out $@.cpp $(OBJECTS) $(STATIC_LIBS) $(VTK_LIBS)

.cpp.o: .cpp .h
	$(CC) $(VTK_INCLUDE) $(FLAGS) -o $@ -c $<
