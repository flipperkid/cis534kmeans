#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <cv.h>
#include <highgui.h>
#include "kmeans.h"
#include "kmeans_tri.h"
#include "kmeans_tbb.h"

int main(int argc, char **argv) {
    if (argc<6) {
        printf("Usage: main cluster_count grainsize thread_count scale mode\n\7");
        exit(0);
    }
    int cluster_count = atoi(argv[1]);
    int grainsize = atoi(argv[2]);
    int thread_count = atoi(argv[3]);
    std::string mode = std::string(argv[5]);

    // load an image  
//    IplImage *img=cvLoadImage("test.jpg");
    IplImage *img=cvLoadImage("portomoniz.jpg");
    if (!img){
        printf("Could not load image file portomoniz.jpg\n");
        exit(0);
    }

    // get the image data
	int imgW = floor(img->width*atof(argv[4]));
	int imgH = floor(img->height*atof(argv[4]));
    IplImage *img2 = cvCreateImage( cvSize(imgW, imgH), 8, img->nChannels );
    cvResize(img, img2);
    int height = img2->height;
    int width = img2->width;
    int particle_count = width*height;
    int channels = img2->nChannels;
    IplImage *hsvImg = cvCreateImage( cvGetSize(img2), 8, img2->nChannels );
    cvCvtColor( img2, hsvImg, CV_BGR2HSV );
    uchar *data = (uchar *)hsvImg->imageData;
    
    // call to kmeans framework
    float *centers = new float[cluster_count*channels];
    uchar *assignments = new uchar[particle_count];
    for( int particle_iter = 0; particle_iter < particle_count; particle_iter++ ) {
        assignments[particle_iter] = UCHAR_MAX;
    }

    if(mode.compare("triangle")) {
        kmeans_tri(data, centers, particle_count, channels, cluster_count, assignments);
    }
    else if(mode.compare("tbb")) {
        kmeans_tbb( data, centers, particle_count, channels, cluster_count, assignments, grainsize, thread_count );
    }
    else {
        kmeans_serial(data, centers, particle_count, channels, cluster_count, assignments);
    }
 
    // release memory
    delete [] assignments; 
    delete [] centers; 
    cvReleaseImage(&hsvImg );
    cvReleaseImage(&img2 );
    cvReleaseImage(&img );
    return 0;
}

