#Overview 
I have used complete c code for this project.The following presents you the way in which you can compile files and run code to get the image.

##C code
It contains my algorithm to find first k singular values.
There are codes for Greyscale(svd_final_greyscale).

##Compiling those  c files
To compile the c code use the following command
gcc new.c -02 -lm -o new
To run the c code to get the output image use the following command
./new inputfilename outputnameyouwant

After compiling the code it creates a shared object(.so files) which are used by python to call c functions

##Conclusion
You can view the image in your file manager
