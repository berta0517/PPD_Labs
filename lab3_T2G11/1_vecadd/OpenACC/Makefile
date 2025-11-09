CC = nvc
CFLAGS = -acc=gpu -gpu=mem:managed -Minfo=all
OBJ = vecadd_oacc
all:
	$(CC) $(CFLAGS) $(OBJ).c -mp -o $(OBJ)
clean:
	rm -f $(OBJ)