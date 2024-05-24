
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>

    // Function to calculate the associated Legendre functions
    double** LegendrePoly(double latitude, int lmax) {
        latitude = latitude * M_PI / 180.0; // Convert latitude to radians

        // Allocate memory for a 2D array plm
        double** plm = (double**)malloc(lmax * sizeof(double*));
        for (int i = 0; i < lmax; i++) {
            plm[i] = (double*)malloc(lmax * sizeof(double));
            for (int j = 0; j < lmax; j++) {
                plm[i][j] = 0.0; // Initialize all elements to 0.0
            }
        }

        plm[0][0] = 1.0;
        if (lmax > 1) {
            plm[1][0] = sqrt(3) * sin(latitude);
            plm[1][1] = sqrt(3) * cos(latitude);
        }
        if (lmax > 2) {
            plm[2][1] = sqrt(5) * plm[1][1] * sin(latitude);
        }

        for (int i = 2; i < lmax; i++) {
            plm[i][i] = sqrt((i + 0.5) / i) * plm[i - 1][i - 1] * cos(latitude);
        }

        for (int j = 2; j < lmax - 1; j++) {
            plm[j + 1][j] = sqrt(2.0 * j + 3) * plm[j][j] * sin(latitude);
        }

        for (int k = 0; k < lmax - 2; k++) {
            for (int l = k + 1; l < lmax - 1; l++) {
                int n = l + 1;
                int m = k + 1;
                double pnm_next = (sqrt((2.0 * n + 1) / (n - (m - 1)) / (n + (m - 1)) * (2.0 * n - 1)) * sin(latitude) * plm[l][k]) -
                                (sqrt((2.0 * n + 1) / (n - (m - 1)) / (n + (m - 1)) * ((n + (m - 1) - 1) / (2.0 * n - 3)) * (n - (m - 1) - 1)) * plm[l - 1][k]);
                plm[l + 1][k] = pnm_next;
            }
        }

        return plm;
    }

    // Function to free the allocated memory for plm
    void freeLegendrePoly(double** plm, int lmax) {
        for (int i = 0; i < lmax; i++) {
            free(plm[i]);
        }
        free(plm);
    }

    int main(int argc, char *argv[]) {
        if (argc != 3) {
            printf("Usage: %s <latitude> <lmax>\n", argv[0]);
            return 1;
        }

        double latitude = atof(argv[1]); // Convert the first argument to a double
        int lmax = atoi(argv[2]); // Convert the second argument to an int

        // Calculate the associated Legendre functions
        double** plm = LegendrePoly(latitude, lmax);

        // Open file to write the results
        FILE *f = fopen("output.txt", "w");
        if (f == NULL) {
            printf("Error opening file!\n");
            return 1;
        }

        // Write the results to the file
        for (int i = 0; i < lmax; i++) {
            for (int j = 0; j < lmax; j++) {
                fprintf(f, "%f ", plm[i][j]);
            }
            fprintf(f, "\n");
        }

        // Close the file
        fclose(f);

        // Free the allocated memory
        freeLegendrePoly(plm, lmax);

        return 0;
    }
    