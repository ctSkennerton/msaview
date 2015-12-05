EXECUTABLE := msaview
OBJS := msaview.o
CFLAGS := -g -O2

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ easel/lib/libeasel.a termbox/lib/libtermbox.a -lm 
