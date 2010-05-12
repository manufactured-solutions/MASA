CC=icpc #intel

SOURCES = masa_map.cpp masa_shell.cpp
OBJECTS = $(addsuffix .o,$(basename $(SOURCES)))
LDFLAGS += #-i-dynamic
CFLAGS  +=
INCLUDE=
LIBS=

## build the executable ##
EXEC=MASA_shell

$(EXEC): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(EXEC) $(INCLUDE) $@ $(LIBS) 

# build cpp object files
%.o: %.cpp
	@echo building $< 
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

# build c object files
%.o: %.c
	@echo building $<
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

# clean directive 'make clean'
clean:
	- /bin/rm $(EXEC) *.o *.mod *~ \#* 
	@echo 'files cleaned'