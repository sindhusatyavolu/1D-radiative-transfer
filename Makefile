CC ?= cc
CFLAGS ?= -Wall -Wextra -std=c11
LDFLAGS ?= -lm

TARGET := rte1d
OPT_TARGET := rte1d_opt
FAST_TARGET := rte1d_fast
FASTER_TARGET := rte1d_faster
OPT_CFLAGS ?= -O3 -march=native -std=c11
FAST_CFLAGS ?= $(OPT_CFLAGS) -DRTE_CACHE_CROSS_SECTIONS
FASTER_CFLAGS ?= $(FAST_CFLAGS) -DRTE_CACHE_FREQUENCY_GRID
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

.PHONY: all optimized fast faster clean

all: $(TARGET)

optimized: $(OPT_TARGET)

fast: $(FAST_TARGET)

faster: $(FASTER_TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) -o $@

$(OPT_TARGET): $(SRC)
	$(CC) $(OPT_CFLAGS) $(SRC) $(LDFLAGS) -o $@

$(FAST_TARGET): $(SRC)
	$(CC) $(FAST_CFLAGS) $(SRC) $(LDFLAGS) -o $@

$(FASTER_TARGET): $(SRC)
	$(CC) $(FASTER_CFLAGS) $(SRC) $(LDFLAGS) -o $@

clean:
	rm -f $(TARGET) $(OPT_TARGET) $(FAST_TARGET) $(FASTER_TARGET)
