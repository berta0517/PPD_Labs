CC = nvc
CFLAGS = -acc=gpu -gpu=managed -Minfo=all
OBJ = partis_oacc_uni_mem
all:
	$(CC) $(CFLAGS) $(OBJ).c -mp -o $(OBJ) -lm
clean:
	rm -f $(OBJ)