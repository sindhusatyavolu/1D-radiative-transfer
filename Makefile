CC ?= cc
CFLAGS ?= -Wall -Wextra -std=c11
LDFLAGS ?= -lm

TARGET := rte1d
OPT_TARGET := rte1d_opt
OPT_CFLAGS ?= -O3 -march=native -std=c11
SRC := main.c \
	src/evolve_cell.c \
	src/lya_spectra.c \
	src/params.c \
	src/photoion.c \
	src/proximity_zone.c \
	src/rates.c \
	src/save_outputs.c \
	src/thermochemistry.c \
	src/utils.c

.PHONY: all optimized clean

all: $(TARGET)

optimized: $(OPT_TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) -o $@

$(OPT_TARGET): $(SRC)
	$(CC) $(OPT_CFLAGS) $(SRC) $(LDFLAGS) -o $@

clean:
	rm -f $(TARGET) $(OPT_TARGET)
