
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"

// Define output file name
#define OUTPUT_FILE "stencil.pgm"

void stencil(const short nx, const short ny, float * restrict image, float * restrict  tmp_image);
void init_image(const short nx, const short ny, float * restrict  image, float * restrict  tmp_image);
void output_image(const char * file_name, const short nx, const short ny, float *image);
double wtime(void);

int main(int argc, char *argv[]) {

  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // Initiliase problem dimensions from command line arguments
  short nx = atoi(argv[1]);
  short ny = atoi(argv[2]);
  short niters = atoi(argv[3]);

  // Allocate the image
  float * restrict image = malloc(sizeof(float)*nx*ny);
  float * restrict tmp_image = malloc(sizeof(float)*nx*ny);

  // Set the input image
  init_image(nx, ny, image, tmp_image);

  //initialise mpi
  MPI_Init( &argc, &argv );

  // Call the stencil kernel
  double tic = wtime();
  for (short t = 0; t < niters; ++t) {
    stencil(nx, ny, image, tmp_image);
    stencil(nx, ny, tmp_image, image);
  }
  double toc = wtime();


  // Output
  printf("------------------------------------\n");
  printf(" runtime: %lf s\n", toc-tic);
  printf("------------------------------------\n");

  output_image(OUTPUT_FILE, nx, ny, image);
  free(image);
}

void stencil(const short nx, const short ny, float * restrict image, float * restrict tmp_image) {
  //when i=0
  //when j=0
  tmp_image[0] = image[0] * 0.6f;
  tmp_image[0] += image[nx] * 0.1f;
  tmp_image[0] += image[0] * 0.1f;
  #pragma vector always
  for (int j = 1; j < ny-1; ++j) {
    tmp_image[j] = image[j] * 0.6f;
    tmp_image[j] += image[j  +nx] * 0.1f;
    tmp_image[j] += image[j-1] * 0.1f;
    tmp_image[j] += image[j+1] * 0.1f;
  }
  //when j=ny-1
  tmp_image[ny-1] = image[ny-1] * 0.6f;
  tmp_image[ny-1] += image[(ny-1) +nx] * 0.1f;
  tmp_image[ny-1] += image[(ny-1)-1] * 0.1f;

  #pragma vector always
  for (int i = 1; i < nx-1; ++i) {
    //when j=0
    tmp_image[i*nx] = image[i*nx] * 0.6f;
    tmp_image[i*nx] += image[(i-1)*nx] * 0.1f;
    tmp_image[i*nx] += image[(i+1)*nx] * 0.1f;
    tmp_image[i*nx] += image[1+i*nx] * 0.1f;
    #pragma vector always
    for (int j = 1; j < ny-1; ++j) {
      tmp_image[j+i*nx] = image[j+i*nx] * 0.6f;
      tmp_image[j+i*nx] += image[j  +(i-1)*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j  +(i+1)*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j-1+i*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j+1+i*nx] * 0.1f;
    }
    //when j=ny-1
    tmp_image[(ny-1)+i*nx] = image[(ny-1)+i*nx] * 0.6f;
    tmp_image[(ny-1)+i*nx] += image[(ny-1)  +(i-1)*nx] * 0.1f;
    tmp_image[(ny-1)+i*nx] += image[(ny-1)  +(i+1)*nx] * 0.1f;
    tmp_image[(ny-1)+i*nx] += image[(ny-1)-1+i*nx] * 0.1f;
  }
  //when i = nx-1
  //when j=0
  tmp_image[(nx-1)*nx] = image[(nx-1)*nx] * 0.6f;
  tmp_image[(nx-1)*nx] += image[((nx-1)-1)*nx] * 0.1f;
  tmp_image[(nx-1)*nx] += image[1+(nx-1)*nx] * 0.1f;
  #pragma vector always
  for (int j = 1; j < ny-1; ++j) {
    tmp_image[j+(nx-1)*nx] = image[j+(nx-1)*nx] * 0.6f;
    tmp_image[j+(nx-1)*nx] += image[j  +((nx-1)-1)*nx] * 0.1f;
    tmp_image[j+(nx-1)*nx] += image[j-1+(nx-1)*nx] * 0.1f;
    tmp_image[j+(nx-1)*nx] += image[j+1+(nx-1)*nx] * 0.1f;
  }
  //when j=ny-1
  tmp_image[(ny-1)+(nx-1)*nx] = image[(ny-1)+(nx-1)*nx] * 0.6f;
  tmp_image[(ny-1)+(nx-1)*nx] += image[(ny-1)  +((nx-1)-1)*nx] * 0.1f;
  tmp_image[(ny-1)+(nx-1)*nx] += image[(ny-1)-1+(nx-1)*nx] * 0.1f;
}

// Create the input image
void init_image(const short nx, const short ny, float * restrict image, float * restrict tmp_image) {
  // Zero everything
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      image[j+i*ny] = 0.0f;
      tmp_image[j+i*ny] = 0.0f;
    }
  }

  // Checkerboard
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      for (int ii = j*nx/8; ii < (j+1)*nx/8; ++ii) {
        for (int jj = i*ny/8; jj < (i+1)*ny/8; ++jj) {
          if ((i+j)%2)
          image[jj+ii*ny] = 100.0f;
        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char * file_name, const short nx, const short ny, float *image) {

  // Open output file
  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  float maximum = 0.0f;
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (image[j+i*ny] > maximum)
        maximum = image[j+i*ny];
    }
  }

  // Output image, converting to numbers 0-255
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      fputc((char)(255.0f*image[j+i*ny]/maximum), fp);
    }
  }

  // Close the file
  fclose(fp);

}

// Get the current time in seconds since the Epoch
double wtime(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec*1e-6;
}
