NVCC = nvcc
NVCCFLAGS = -O2 -arch=native

TARGET = vecadd_cuda

all:
	$(NVCC) $(TARGET).cu $(NVCCFLAGS) -o $(TARGET)

clean:
	rm -f $(TARGET)

