CC=nvcc
FLAGS=-arch=sm_20
OBJECTS = main.o hostfunctions.o devicefunctions.o

gpuCode: $(OBJECTS)
	$(CC) $(FLAGS) $(OBJECTS) -o $@.out

%.o: %.cpp
	nvcc -x cu -arch=sm_20 -I. -dc $< -o $@

clean:
	rm -rf ./*.o ./*.out
