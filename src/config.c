#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/params.h"

void read_params_from_config(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening config file");
        exit(1);
    }

    char line[256], key[64];
    double val;

    while (fgets(line, sizeof(line), fp)) {
        // skip comments and blank lines
        if (line[0] == '#' || line[0] == ';' || line[0] == '\n')
            continue;

        // remove inline comments (// or ;)
        char *comment = strpbrk(line, ";#");
        if (comment) *comment = '\0';

        if (sscanf(line, "%63[^=]=%lf", key, &val) == 2) {
            // match keys (no whitespace around key)
            if (strcmp(key, "X") == 0) X = val;
            else if (strcmp(key, "Y") == 0) Y = val;
            else if (strcmp(key, "omega_m") == 0) omega_m = val;
            else if (strcmp(key, "H0") == 0) H0 = val;

            else if (strcmp(key, "z") == 0) z = val;
            else if (strcmp(key, "z_qso") == 0) z_qso = val;
            else if (strcmp(key, "Ng") == 0) Ng = val;
            else if (strcmp(key, "epsilon") == 0) epsilon = val;
            else if (strcmp(key, "t_on") == 0) tq = val;
            else if (strcmp(key, "fduty") == 0) fduty = val;
            else if (strcmp(key, "tcycle") == 0) t_on = val;
            else if (strcmp(key, "NLOS") == 0) NLOS = (int)val;
            else if (strcmp(key, "M1450") == 0) M1450 = val;
            else if (strcmp(key, "alpha") == 0) alpha = val;
            else if (strcmp(key, "alpha2") == 0) alpha2 = val;
            else if (strcmp(key, "fend") == 0) fend = val;
            else if (strcmp(key, "NFREQ") == 0) NFREQ = (int)val;

            else if (strcmp(key, "Ncells") == 0) Ncells = (int)val;
            else if (strcmp(key, "cellsize") == 0) cellsize = val;
            else if (strcmp(key, "boxsize") == 0) boxsize = val;
            else if (strcmp(key, "hubble_inv") == 0) hubble_inv = val;
            else if (strcmp(key, "Roffset") == 0) Roffset = val;

            else if (strcmp(key, "background") == 0) background = val;
        }
    }

    fclose(fp);
}
