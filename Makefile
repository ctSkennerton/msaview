EXECUTABLE := msaview
OBJS := msaview.o

$(EXECUTABLE): $(OBJS)
	$(CC) -o $@ $^ easel/lib/libeasel.a -lm -lncurses
