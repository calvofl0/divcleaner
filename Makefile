SOLIB   := divcleaner_idl.so
SRCS    := bicgstab.c divcleaner_idl.c

CC      := cc
CFLAGS  := -fPIC -Wall -Wextra -Werror
LDFLAGS := -shared

OBJS     = $(SRCS:.c=.o)
DEPS     = $(SRCS:.c=.d)

all: $(SOLIB)

$(SOLIB): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -MMD -c $< -o $@

-include $(DEPS)

.PHONY: all clean

clean:
	rm -f $(OBJS) $(DEPS)

cleaner: clean
	rm -rf $(SOLIB)
