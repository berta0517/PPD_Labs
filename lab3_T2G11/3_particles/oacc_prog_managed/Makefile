CC = nvc
CFLAGS = -acc=gpu -gpu=cc90 -Minfo=all
OBJ = partis_oacc_prog_managed
all:
	$(CC) $(CFLAGS) $(OBJ).c -mp -o $(OBJ) -lm
clean:
	rm -f $(OBJ)