/*
* Adapted from:
* https://github.com/walterdejong/glwhyz
*
* Which in turn used code from:
* see http://tfcduke.developpez.com/tutoriel/format/tga/
*
* Usage:
*    TGA *image = load_tga("img01.tga");
*    ... do something cool ...
*    free(image);
*/
#include "mpi.h"
#include "tga.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/******************************************************************************/

TGA *load_tga(const char *filename, TGA_Header *header) {
  TGA *image;
  FILE *f;

  if ((filename == NULL) || (!*filename)) {
    return (NULL);
  }

  if ((f = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "load_tga(): failed to open file '%s'\n", filename);
    return (NULL);
  }

  /* read header */
  fread(header, sizeof(TGA_Header), 1, f);

  if (feof(f)) {
    fprintf(stderr, "load_tga(): failed to load '%s'\n", filename);
    fclose(f);
    return (NULL);
  }

  if (fseek(f, header->id_lenght, SEEK_CUR) == -1) {
    fprintf(stderr, "load_tga(): failed to load '%s'\n", filename);
    fclose(f);
    return (NULL);
  }

  /* allocate pixel data right after the struct TGA */
  if ((image = (TGA *)malloc(sizeof(TGA) +
                             header->width * header->height *
                                 header->pixel_depth / 8)) == NULL) {
    fprintf(stderr, "load_tga(): out of memory while loading '%s'\n", filename);
    fclose(f);
    return (NULL);
  }
  /*else {
  printf("%u",header->pixel_depth / 8);
}*/

  /* read image data */
  if ((header->image_type == 3) && (header->pixel_depth == 8)) {
    /*
    * After these lines the image is contained in the matrix 'pixels'.
    */
    image->w = header->width;
    image->h = header->height;
    fread(image->pixels, image->w, image->h, f);
  } else {
    fprintf(stderr, "load_tga(): '%s' is not an 8-bit grayscale TGA image\n",
            filename);
    fclose(f);
    free(image);
    return (NULL);
  }

  fclose(f);
  return (image);
}

/******************************************************************************/

void save_tga(char *filename, TGA_Header header, TGA *image) {
  FILE *f;

  if ((f = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "save_tga(): failed to open file '%s'\n", filename);
    return;
  }

  /* write header */
  fwrite(&header, sizeof(TGA_Header), 1, f);

  fwrite(image->pixels, image->w, image->h, f);

  fclose(f);
  // printf("\n**Output Image Saved Successfully\n");
  printf("\n**Blured Image Saved Successfully , name : %s\n", filename);
}

/******************************************************************************/

void free_tga(TGA *image) {
  if (image != NULL) {
    free(image);
  }
}

/******************************************************************************/
allocator(int height, int width) {
  int *data = (int *)malloc(height * width * sizeof(int));
  int **array = (int **)malloc(height * sizeof(int *));
  for (int i = 0; i < height; i++) {
    array[i] = &(data[width * i]);
  }
  return array;
}

int **blur_filter(int ht, int width, int extra, int firstRow, int lastRow,
                  int array[][width], int blur[][width], int my_rank) {

  int i = 0, j = 0, count = 0, sum = 0, c = 0, d = 0, avg = 0, counter = 0,
      height = 0;
  height = ht + extra;
  // Gauss blur filter implementation

  for (i = firstRow; i < height - lastRow; i++) {
    for (j = 0; j < width; j++) {
      count = 0;
      sum = 0;
      // Iterate to find Neighbors
      for (c = i - 1; c <= i + 1; c++) {
        for (d = j - 1; d <= j + 1; d++) {
          // Look for elements only within bounds!
          if ((c >= 0 && c <= height - 1) && (d >= 0 && d <= width - 1)) {
            // Rule out Self element for each calculation
            count++;
            sum = sum + array[c][d];
          }
        }
      }
      avg = sum / count;
      blur[counter][j] = avg;
    }

    counter++;
  }
  return blur;
}

int *transform_array(int height, int width, int array[][width], int *result) {
  int upper = 0, i = 0, j = 0, count = 0;

  if (result == NULL) {
    printf("ALLOCATION FAILED!");
  }

  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      // image->pixels[upper] = result[i][j];

      result[upper] = array[i][j];
      // printf("%d,%d\t", upper, result[upper]);
      upper++;
      count++;
      if (count == width) {
        // printf("\n");
        count = 0;
      }
    }
  }

  // free(result);
  return result;
}

histogram(int *array, int height, int width, int recursion) {
  int histogram[256], position, i = 0, j = 0;
  char histogram_name[100];
  for (i = 0; i < 256; i++) {
    histogram[i] = 0;
  }

  for (i = 0; i < height * width; i++) {
    position = array[i];
    histogram[position]++;
  }
  /* strcpy(histogram_name, "histogram");
   strcat(histogram_name, r);*/
  sprintf(histogram_name, "histogram %d", recursion);

  FILE *f = fopen(histogram_name, "w");
  if (f == NULL) {
    printf("Error opening file!\n");
    exit(1);
  }
  for (i = 0; i < 256; i++) {
    fprintf(f, "%d\n", histogram[i]);
  }
  printf("\n**Histogram Saved Successfully , name : %s\n", histogram_name);
  free(f);
}

void p() { printf("\n"); }

print_array_1d(int *A, int height, int width) {
  int count = 0;
  for (int i = 0; i < height * width; i++) {
    printf("%d\t", A[i]);
    count++;
    if (count == width) {
      printf("\n");
      count = 0;
    }
  }
  p();
}

print_array(int *A, size_t height, size_t width, int process) {
  int count = 0;
  for (size_t i = 0; i < height; ++i) {
    for (size_t j = 0; j < width; ++j) {

      // printf("A[%zu][%zu] = %d\t", i, j, A[i * width + j]);
      // printf("%d:%d\t", process, A[i * width + j]);
      count++;
      if (count == width) {
        printf(",");
        count = 0;
      }
    }
  }
  printf("\n");
}
/******************************************************************/

int main(int argc, char *argv[]) {
  TGA_Header header;
  char output_filename[100];

  if (argc != 3) {
    printf("\nUsage:\n");
    printf("\t./tga <filename>.tga numberOfRepeats\n\n");
    printf(
        "The input file should be a greyscale image with an 8-bit depth.\n\n");
    exit(0);
  }

  TGA *image = load_tga(argv[1], &header);

  /*
  * Image processing starts here
  */
  int my_rank = 0, size = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Variable Initialization
  int i = 0, j = 0, c = 0, d = 0, avg = 0, count = 0, target = 0, sum = 0,
      sliceSizeExt = 0, upper = 0, start = 0, end = 0, extrapixels = 0,
      sliceSize = 0, height = 0, width = 0, RECURSION, r, **array, **result,
      *RECVBUF, ierr, errclass, resultlen;
  time_t t0, t1;
  if (my_rank == 0) {
    gettimeofday(&t0, 0);
  }

  MPI_Status status;
  // Set Height and width!
  // width = image->w;
  width = image->w;
  height = image->h;
  RECURSION = atoi(argv[2]);

  // Dynamic allocation for both arrays

  /*
  First and lastRow arrays should be allocated with +1 line for the respected
ghost
borders
  Arrays for processes with rank 0 > my_rank > size should keep the same
allocation

  Exception:
  If the height/size modulo is (>) Greater than (0) Zeroprintf("\n Computing
bounds for process ")
}
Then it should be processed by the rank 0
Just to keep the communication dynamic.
*/
  sliceSize = height / size;
  extrapixels = height % size;
  sliceSizeExt = sliceSize + extrapixels;

  // Contiguous allocation for MPI Distribution

  int(*subArrayLast)[width] = malloc(sizeof(int[(sliceSize + 1)][width]));
  int(*subArray)[width] = malloc(sizeof(int[(sliceSize + 2)][width]));
  int(*subResult)[width] = malloc(sizeof(int[sliceSize][width]));
  int(*subArrayFext)[width] = malloc(sizeof(int[sliceSizeExt + 1][width]));
  int(*subResultFext)[width] = malloc(sizeof(int[sliceSizeExt][width]));

  if (subArrayLast == NULL) {
    printf("Wrong1");
    return 0;
  }

  if (subArray == NULL) {
    printf("Wrong2");
    return 0;
  }
  if (subArrayFext == NULL) {
    printf("Wrong3");
    return 0;
  }
  if (subResultFext == NULL) {
    printf("Wrong4");
    return 0;
  }
  if (subResult == NULL) {
    printf("Wrong5");
    return 0;
  }
  int *SDresult = malloc(sliceSize * width * sizeof(int));
  if (SDresult == NULL) {
    printf("Wrong7");
    return 0;
  }

  int *SDresultFext = malloc(sliceSizeExt * width * sizeof(int));

  if (SDresultFext == NULL) {
    printf("Wrong6");
    return 0;
  }
  /**/

  RECVBUF = (int *)malloc(sizeof(int) * height * width);
  if (RECVBUF == NULL) {
    printf("Wrong8");
    return 0;
  }

  result = array = malloc(sizeof(int *) * height);
  for (i = 0; i < height; i++) {
    array[i] = malloc(sizeof(int) * width);
  }

  if (array == NULL) {
    printf("Wrong9");
    return 0;
  }

  /*Some useful printouts!*/
  if (my_rank == 0) {
    printf("\nThis will run for %d times\n", RECURSION);
    printf("\n**************************************************************"
           "****"
           "\nRANK:%d \nSome statistics: \nSlave "
           "Processes:%d\nImage dimensions %d x %d\nArray "
           "Slice Size(height): %d\nThese are the extra pixels! : %d\nThis "
           "is the "
           "Ext "
           "dimension:%d\n\nMoving on to Compute slice "
           "bounds...\n******************************************************"
           "*****"
           "********\n\n\n\n",
           my_rank, size - 1, height, width, sliceSize, extrapixels,
           sliceSizeExt);
  }

  /*space allocated!*/
  /*master process*/

  for (r = 0; r < RECURSION; r++) {
    if (my_rank == 0) {
      upper = 0;

      // Parse image pixel array to 2d array
      for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
          // use incremental counter for 1d array(: image->pixels[width*height])
          // to 2d transformation
          array[i][j] = image->pixels[upper];
          upper++;
        }
      }
      // print_array(*array, height, width, my_rank);

      // Compute the bounds for processes 1-2-3

      for (target = 0; target < size; target++) {
        if (target == 0) {
          start = target * sliceSize;
          end = start + (sliceSize + extrapixels);
        } else if (target < size - 1) {
          start = (target * sliceSize) + extrapixels - 1;
          end = start + sliceSize + 1;
        } else {
          start = (target * sliceSize) + extrapixels - 1;
          end = start + sliceSize;
        }
        if (size == 1) {
          end--;
        }
        for (i = start; i <= end; i++) {
          for (j = 0; j < width; j++) {
            if (target == 0) {
              subArrayFext[i - start][j] = array[i][j];
            } else if (target < size - 1) {
              subArray[i - start][j] = array[i][j];
            } else {
              subArrayLast[i - start][j] = array[i][j];
            }
            // printf("array: %d ",subArray[i-start][j]);
          }
        }

        if (target > 0 && target < size - 1) {

          MPI_Send(subArray, (sliceSize + 2) * width, MPI_INT, target, 10,
                   MPI_COMM_WORLD);
        } else if ((target == size - 1) && (target != 0)) {

          MPI_Send(subArrayLast, (sliceSize + 1) * width, MPI_INT, target, 10,
                   MPI_COMM_WORLD);
        }
      }
      /*EVERYTHING IS OK */

      int(*blur)[width] = malloc(sizeof(int[sliceSizeExt][width]));
      subResultFext =
          blur_filter(sliceSizeExt, width, 1, 0, 1, subArrayFext, blur, 0);
      // printf("bluring");
      int *trans = malloc(sliceSizeExt * width * sizeof(int));
      SDresultFext = transform_array(sliceSizeExt, width, subResultFext, trans);

      for (target = 0; target < size; target++) {

        start = target * (sliceSize * width) + (extrapixels * width);
        end = start + sliceSize * width + 1;
        if (target == size - 1) {
          end = (height * width) - 1;
        }
        // printf("\n %d %d\n", start, end);

        if (target == 0) {
          for (int i = 0; i < sliceSizeExt * width; i++) {
            RECVBUF[i] = SDresultFext[i];
          }
        } else if (target > 0) {
          MPI_Recv(SDresult, sliceSize * width, MPI_INT, target, 20,
                   MPI_COMM_WORLD, &status);

          // print_array_1d(SDresult, sliceSize, width);
          for (i = start; i <= end; i++) {
            RECVBUF[i] = SDresult[count];
            count++;
          }
          count = 0;
        }
      }
      for (i = 0; i < height * width; i++) {
        image->pixels[i] = RECVBUF[i];
      }

      /*histogram(RECVBUF, height, width, r);
      strcpy(output_filename, "output_");
      strcat(output_filename, argv[1]);
      sprintf(output_filename, "%d", r);
      save_tga(output_filename, header, image);*/
      /*************************************************************************************************************************************************/
      free(blur);
      free(trans);
    } else {

      if (my_rank < size - 1) {
        MPI_Recv(subArray, (sliceSize + 2) * width, MPI_INT, 0, 10,
                 MPI_COMM_WORLD, &status);
        int(*blur)[width] = malloc(sizeof(int[sliceSize][width]));
        subResult = blur_filter(sliceSize, width, 2, 1, 1, subArray, blur, 1);

        /*int **blur_filter(int ht, int width, int extra, int firstRow, int
           lastRow,
                          int array[][width], int blur[][width], int my_rank) */

        // print_array(blur, sliceSize, width, my_rank);
        int *trans = malloc(sliceSize * width * sizeof(int));
        SDresult = transform_array(sliceSize, width, subResult, trans);
        MPI_Send(SDresult, sliceSize * width, MPI_INT, 0, 20, MPI_COMM_WORLD);
        free(blur);
        free(trans);
      } else if (my_rank == size - 1) {
        MPI_Recv(subArrayLast, (sliceSize + 1) * width, MPI_INT, 0, 10,
                 MPI_COMM_WORLD, &status);
        int(*blur)[width] = malloc(sizeof(int[sliceSize][width]));

        subResult =
            blur_filter(sliceSize, width, 1, 1, 0, subArrayLast, blur, 1);

        // printf("bluring");
        int *trans = malloc(sliceSize * width * sizeof(int));
        SDresult = transform_array(sliceSize, width, subResult, trans);
        // print_array_1d(SDresult, sliceSize, width);
        MPI_Send(SDresult, sliceSize * width, MPI_INT, 0, 20, MPI_COMM_WORLD);
        free(blur);
        free(trans);
      }
    }
  }

  /*
  * Save processed image to a new file.
  * The TGA header for the new file is the same as for the input file.
  */
  if (my_rank == 0) {

    strcpy(output_filename, "output_");
    strcat(output_filename, argv[1]);

    save_tga(output_filename, header, image);
    gettimeofday(&t1, 0);
    long elapsed = (t1 - t0) * 1000000 + t1 - t0;
    printf("\n\n\n\n**********************************************************"
           "********* \n=REPORTING\n=time elapsed : %ld or approx. "
           "%ldseconds\n=histograms are saved in the code DIR. histogram 0 / "
           "output 0\n"
           "**********************************************************"
           "*********\n",
           elapsed, elapsed / 1000000);
  }
  // printf("\n...OK\n");
  free(image);

  MPI_Finalize();
  return (0);
}
