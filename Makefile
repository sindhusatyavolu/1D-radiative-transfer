CC ?= cc
CFLAGS ?= -Wall -Wextra -std=c11
LDFLAGS ?= -lm

TARGET := rte1d
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

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) -o $@

clean:
	rm -f $(TARGET)
