#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cv.h>
#include <highgui.h>
#include "kmeans.h"
#include "kmeans_omp.h"

int output_assignments(IplImage *img, uchar *assignments, double *centers, int cluster_count, int particle_count, int dimensions);

int main(int argc, char **argv) {
    if (argc<2) {
        printf("Usage: main <image-file-name>\n\7");
        exit(0);
    }

    // load an image  
    IplImage *img=cvLoadImage(argv[1]);
    if (!img){
        printf("Could not load image file: %s\n",argv[1]);
        exit(0);
    }

    // get the image data
    int height = img->height;
    int width = img->width;
    int particle_count = width*height;
    int channels = img->nChannels;
    IplImage *hsvImg = cvCreateImage( cvGetSize(img), 8, 3 );
    cvCvtColor( img, hsvImg, CV_BGR2HSV );
    uchar *data = (uchar *)hsvImg->imageData;
    
    // call to kmeans framework
    int cluster_count = 32;
    double *centers = new double[cluster_count*channels];
    uchar *assignments = new uchar[particle_count];
    for( int particle_iter = 0; particle_iter < particle_count; particle_iter++ ) {
        assignments[particle_iter] = UCHAR_MAX;
    }
    //kmeans_serial(data, centers, particle_count, channels, cluster_count, assignments);
    int grainSize = 1024;
    int thread_count = 4;
    kmeans_omp(data, centers, particle_count, channels, cluster_count, assignments, grainSize, thread_count);

    IplImage *outImg = cvCreateImage( cvGetSize(img), 8, 3 );
    output_assignments(outImg, assignments, centers, cluster_count, particle_count, channels);

    // release memory
    delete [] assignments; 
    delete [] centers; 
    cvReleaseImage(&outImg );
    cvReleaseImage(&hsvImg );
    cvReleaseImage(&img );
    return 0;
}

int output_assignments(IplImage *img, uchar *assignments, double *centers, int cluster_count, int particle_count, int dimensions) {
    IplImage *tempImg = cvCreateImage( cvSize(cluster_count, 1), 8, 3);
    for( int centers_iter = 0; centers_iter < cluster_count*3; centers_iter++ ) {
        tempImg->imageData[centers_iter] = (char)floor(centers[centers_iter]);
    }
    cvCvtColor( tempImg, tempImg, CV_HSV2BGR );

    for( int particle_iter = 0; particle_iter < particle_count; particle_iter++ ) {
        int assn = (int)assignments[particle_iter];
        for( int dim_iter = 0; dim_iter < dimensions; dim_iter++ ) {
//            printf("here %f", array_load<double>(centers, assn, dim_iter, dimensions));
            img->imageData[particle_iter*dimensions+dim_iter] = array_load<char>(tempImg->imageData, assn, dim_iter, dimensions);
        }
    }
    
    char outFileName[] = "seg.png";
    if (!cvSaveImage(outFileName,img)) {
        printf("Could not save: %s\n",outFileName);
    }
    cvReleaseImage(&tempImg);
    return 0;
}
