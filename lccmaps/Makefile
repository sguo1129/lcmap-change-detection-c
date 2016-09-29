# Set up compile options
CC    = gcc
RM    = rm -f
MV    = mv
CFLAGS = -Wall -g

# Define the include files
INC = utilities.h lccmaps.h 2d_array.h 
INCDIR  = -I. -I$(MATIO_INC) -I$(HDF5INC)
NCFLAGS = $(CFLAGS) $(INCDIR)

# Define the source code and object files
SRC = utilities.c lccmaps.c input.c output.c misc.c 2d_array.c
OBJ = $(SRC:.c=.o)

# Define the object libraries
LIB = -lm -lz -lpthread -lrt -L$(MATIO_LIB) -L$(HDF5LIB) -lmatio -lhdf5 -lhdf5_hl

# Define the executable
EXE = lccmaps

# Target for the executable
all: $(EXE)

lccmaps: $(OBJ) $(INC)
	$(CC) $(NCFLAGS) -o lccmaps $(OBJ) $(LIB)

install:
	cp $(EXE) $(BIN)

clean:
	$(RM) *.o $(EXE)

$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) $(INCDIR) -c $<

