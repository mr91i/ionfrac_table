###
CC:=gcc
CFLAGS := -g -Wall  -Wextra -O3
OBJS := main_make_table_2.o const.o NLOhm.o
LIBS := -lm
TARGET := mktbl.out

###
all: clean $(TARGET)

###
$(TARGET): $(OBJS)
	$(CC) $(LIBS) -o $@ $(OBJS) 

clean:
	$(RM) $(TARGET) $(OBJS)
###
.c.o:
	$(CC) $(CFLAGS) -c $<

###
main_make_table_2.o NLOhm.o OHP.o: Headers.h

