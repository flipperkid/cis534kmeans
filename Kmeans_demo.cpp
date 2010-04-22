#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cv.h>
#include <highgui.h>
#include "kmeans.h"

int output_assignments(IplImage *img, uchar *assignments, int cluster_count, int particle_count);

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
    int cluster_count = 25;
    uchar *assignments = new uchar[particle_count];
    kmeans_serial(data, particle_count, channels, cluster_count, assignments);

    IplImage *outImg = cvCreateImage( cvGetSize(img), 8, 1 );
    output_assignments(outImg, assignments, cluster_count, particle_count);

    // release memory
    delete [] assignments; 
    cvReleaseImage(&hsvImg );
    cvReleaseImage(&img );
    return 0;
}

int output_assignments(IplImage *img, uchar *assignments, int cluster_count, int particle_count) {
    for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) {
        assignments[particle_iter] *= floor(255/cluster_count);
    }
    
    img->imageData = (char *)assignments;
    char outFileName[] = "seg.png";
    if (!cvSaveImage(outFileName,img)) {
        printf("Could not save: %s\n",outFileName);
    }
}
