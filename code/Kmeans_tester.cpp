#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cv.h>
#include <highgui.h>
#include "kmeans.h"
#include "kmeans_tbb.h"

int main(int argc, char **argv) {
    if (argc<2) {
        printf("Usage: main grainsize\n\7");
        exit(0);
    }
    int grainsize = atoi(argv[1]);

    // load an image  
    IplImage *img=cvLoadImage("portomoniz.jpg");
    if (!img){
        printf("Could not load image file portomoniz.jpg\n");
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
    int cluster_count = 50;
    double *centers = new double[cluster_count*channels];
    uchar *assignments = new uchar[particle_count];
    for( int particle_iter = 0; particle_iter < particle_count; particle_iter++ ) {
        assignments[particle_iter] = UCHAR_MAX;
    }
//    kmeans_serial(data, centers, particle_count, channels, cluster_count, assignments);
    kmeans_tbb( data, centers, particle_count, channels, cluster_count, assignments, grainsize );
 
    // release memory
    delete [] assignments; 
    delete [] centers; 
    cvReleaseImage(&hsvImg );
    cvReleaseImage(&img );
    return 0;
}
