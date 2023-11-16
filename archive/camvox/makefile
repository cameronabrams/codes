# general makefile for md series
# (c) 1997, 1998 Cam Abrams

C_SRCS_R3D   =  dblmat.c point.c colormap.c r3d_utils.c

INC_DIR = .

BIN_DIR = .

OBJ_DIR = .

OBJS_R3D      = $(C_SRCS_R3D:.c=.o)

CC = cc

RM = rm -f


CFLAGS = -O3 -I$(INC_DIR) 

LFLAGS = 
LFLAGS_GENERAL = 


camvox: $(OBJS_R3D) camvox.o
	$(CC) $(LFLAGS) -o $(BIN_DIR)/camvox $(OBJS_R3D) camvox.o -lm
