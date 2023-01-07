# Software Project - Final Project: K-means Clustering and Matrix Calculations
This is the final project for the course "Project Software" at Tel-Aviv University.

## Compiling
To compile the program, use the following command:  
gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c -lm -o spkmeans

## Usage
The program has the following usage:  
./spkmeans [goal] [input_file]  
Where [goal] is one of the following strings:  
wam: Calculate the within-class average matrix (WAM)  
ddg: Calculate the difference of within-class and between-class scatter matrices (DDG)  
lnorm: Calculate the normalized DDG matrix (Lnorm)  
jacobi: Calculate the eigenvalues and eigenvectors of a symmetric matrix using the Jacobi method  
kmeanspp: Perform K-means clustering using the K-means++ initialization method  
[input_file] is the path to the input file containing the matrix or vectors to be used for the calculations.  

The input file should be in the following format:  
For the wam, ddg, and lnorm goals, the input file should contain a square matrix of size NxN with N being the number of rows (and columns). The elements of the matrix should be separated by commas.  
For the jacobi goal, the input file should contain a symmetric matrix of size NxN with N being the number of rows (and columns). The elements of the matrix should be separated by commas.  
For the kmeanspp goal, the input file should contain a matrix of size Nxd, where N is the number of rows and d is the number of dimensions. The elements of the matrix should be separated by commas.  
The program will then perform the desired calculation and print the result to the standard output.

## Grading
This program was graded 100.

## Notes
The program includes additional functions that were not specified in the original prompt, such as the calculation of the normalized Laplacian matrix and the Jacobi method for calculating eigenvalues and eigenvectors. These functions were implemented for additional functionality and as a learning exercise.  
The program also includes a K-means implementation, which was not a requirement for the original prompt. This implementation was included for additional functionality and as a learning exercise.
