CC = nvc
CFLAGS = -acc=gpu -gpu=cc90 -Minfo=all
OBJ = partis_oacc_async
all:
	$(CC) $(CFLAGS) $(OBJ).c -mp -o $(OBJ) -lm
clean:
	rm -f $(OBJ)