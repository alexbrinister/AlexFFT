#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Generic FFT function using _Generic
#define fft(in, out, N) _Generic((in), \
    float *: fft_float, \
    double *: fft_double \
)(in, out, N) 

// FFT implementation for float
void fft_float(float *in, float complex *out, size_t N) {
    // Check if N is a power of 2
    if (N <= 1) {
        out[0] = in[0];
        return;
    }
    if (N & (N - 1)) {
        fprintf(stderr, "Error: N must be a power of 2 for FFT.\n");
        exit(1);
    }

    // Split the input into even and odd parts
    float *even = malloc(N / 2 * sizeof(float));
    float *odd = malloc(N / 2 * sizeof(float));
    for (size_t i = 0; i < N / 2; i++) {
        even[i] = in[2 * i];
        odd[i] = in[2 * i + 1];
    }

    // Recursively compute the FFT of the even and odd parts
    float complex *even_fft = malloc(N / 2 * sizeof(float complex));
    float complex *odd_fft = malloc(N / 2 * sizeof(float complex));
    fft_float(even, even_fft, N / 2);
    fft_float(odd, odd_fft, N / 2);

    // Combine the results
    for (size_t k = 0; k < N / 2; k++) {
        float angle = -2 * M_PI * k / N;
        float complex wk = cexpf(angle * I);
        out[k] = even_fft[k] + wk * odd_fft[k];
        out[k + N / 2] = even_fft[k] - wk * odd_fft[k];
    }

    // Free allocated memory
    free(even);
    free(odd);
    free(even_fft);
    free(odd_fft);
}

// FFT implementation for double (similar structure to fft_float)
void fft_double(double *in, double complex *out, size_t N) {
    if (N <= 1) {
        out[0] = in[0];
        return;
    }
    if (N & (N - 1)) {
        fprintf(stderr, "Error: N must be a power of 2 for FFT.\n");
        exit(1);
    }

    double *even = malloc(N / 2 * sizeof(double));
    double *odd = malloc(N / 2 * sizeof(double));
    for (size_t i = 0; i < N / 2; i++) {
        even[i] = in[2 * i];
        odd[i] = in[2 * i + 1];
    }

    double complex *even_fft = malloc(N / 2 * sizeof(double complex));
    double complex *odd_fft = malloc(N / 2 * sizeof(double complex));
    fft_double(even, even_fft, N / 2);
    fft_double(odd, odd_fft, N / 2);

    for (size_t k = 0; k < N / 2; k++) {
        double angle = -2 * M_PI * k / N;
        double complex wk = cexp(angle * I);
        out[k] = even_fft[k] + wk * odd_fft[k];
        out[k + N / 2] = even_fft[k] - wk * odd_fft[k];
    }

    free(even);
    free(odd);
    free(even_fft);
    free(odd_fft);
}

int main() {
    // Example usage with float data
    float wave_data[] = {1.0, 2.0, 3.0, 4.0};
    size_t N = sizeof(wave_data) / sizeof(wave_data[0]);
    float complex *fft_result = malloc(N * sizeof(float complex));

    // Call the generic fft function
    fft(wave_data, fft_result, N);

    // Print the FFT result
    for (size_t i = 0; i < N; i++) {
        printf("X(%zu) = %f + %fi\n", i, crealf(fft_result[i]),
               cimagf(fft_result[i]));
    }

    free(fft_result);
    return 0;
}
