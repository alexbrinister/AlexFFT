#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Generic DFT function using _Generic
#define dft(in, out, N) _Generic((in), \
    float *: dft_float, \
    double *: dft_double  \
)(in, out, N) 

// DFT implementation for float
void dft_float(float *in, float complex *out, size_t N) {
    for (size_t k = 0; k < N; k++) {
        out[k] = 0;
        for (size_t n = 0; n < N; n++) {
            float angle = 2 * M_PI * k * n / N;
            out[k] += in[n] * cexpf(-I * angle);
        }
    }
}

// DFT implementation for double
void dft_double(double *in, double complex *out, size_t N) {
    for (size_t k = 0; k < N; k++) {
        out[k] = 0;
        for (size_t n = 0; n < N; n++) {
            double angle = 2 * M_PI * k * n / N;
            out[k] += in[n] * cexp(-I * angle);
        }
    }
}

int main() {
    // Example usage with float data
    float wave_data[] = {1.0, 2.0, 3.0, 4.0};
    size_t N = sizeof(wave_data) / sizeof(wave_data[0]);
    float complex *dft_result = malloc(N * sizeof(float complex));

    // Call the generic dft function
    dft(wave_data, dft_result, N);

    // Print the DFT result
    for (size_t i = 0; i < N; i++) {
        printf("X(%zu) = %f + %fi\n", i, crealf(dft_result[i]),
               cimagf(dft_result[i]));
    }

    free(dft_result);
    return 0;
}
