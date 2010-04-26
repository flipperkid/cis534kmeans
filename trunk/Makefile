all:	Kmeans_demo Kmeans_tester

Kmeans_demo: 	Makefile code/*
	g++ -o Kmeans_demo code/kmeans.cpp code/Kmeans_demo.cpp code/timer.cpp code/kmeans_tbb.cpp -I ~cis534/public/tbb/include/ -L ~cis534/public/tbb/build/linux_intel64_gcc_cc4.3.2_libc2.9_kernel2.6.27.29_release -I opencv/include/opencv -L opencv/lib64 -lcv -lcvaux -lcxcore -lhighgui -lml -ltbb -O3

Kmeans_tester: 	Makefile code/*
	g++ -o Kmeans_tester code/kmeans.cpp code/Kmeans_tester.cpp code/timer.cpp code/kmeans_tbb.cpp -I ~cis534/public/tbb/include/ -L ~cis534/public/tbb/build/linux_intel64_gcc_cc4.3.2_libc2.9_kernel2.6.27.29_release -I /project/cis/cis534-10a/users/wjc/cis534kmeans/opencv/include/opencv -L /project/cis/cis534-10a/users/wjc/cis534kmeans/opencv/lib64 -lcv -lcvaux -lcxcore -lhighgui -lml -ltbb -O3
