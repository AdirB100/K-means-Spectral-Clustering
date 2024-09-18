# K-means Spectral Clustering

This repository contains the implementation of the **Normalized Spectral Clustering Algorithm** as part of the final project for the Software Project course (0368-2161) at Tel Aviv University. The project includes both Python and C components, with a focus on applying spectral clustering using K-means++ initialization.

## Project Overview

The normalized spectral clustering algorithm is an advanced method for clustering a set of data points by representing them as a graph. This repository implements the algorithm in two stages:

1. **Graph Representation and Spectral Decomposition (C)**: 
   - Construct the weighted adjacency matrix of the data points.
   - Compute the normalized graph Laplacian.
   - Find the eigenvalues and eigenvectors using the **Jacobi algorithm**.
   - Apply the **eigengap heuristic** to determine the optimal number of clusters.
  
2. **K-means Clustering (Python)**: 
   - Perform K-means++ clustering on the resulting eigenvectors to classify the data points into clusters.

### Key Features:
- Implements the **Jacobi algorithm** to compute eigenvalues and eigenvectors of a symmetric matrix.
- Uses **K-means++ initialization** for clustering, as in the course's homework assignments.
- Modular design, with the core spectral decomposition written in C and interfaced with Python using a C API.

## Algorithm

The normalized spectral clustering algorithm proceeds as follows:

1. **Graph Construction**:  
   Form the weighted adjacency matrix \(W\) based on the pairwise Euclidean distances between data points.

2. **Graph Laplacian**:  
   Compute the normalized graph Laplacian \(L_{norm}\).

3. **Eigenvalue Decomposition**:  
   Use the Jacobi algorithm to compute the eigenvalues and eigenvectors of \(L_{norm}\).

4. **Clustering**:  
   Apply the eigengap heuristic to determine the number of clusters \(k\), then cluster the eigenvectors using K-means++.

## Repository Structure

- **spkmeans.py**: Python interface for handling command-line arguments and invoking the spectral clustering process.
- **spkmeans.h**: C header file with function declarations.
- **spkmeans.c**: C code that performs matrix computations and spectral decomposition.
- **spkmeansmodule.c**: Python C API wrapper that links the C and Python components.
- **setup.py**: Setup file to build the C extension.
- **comp.sh**: Script for compiling the C program.
- **data/**: Directory containing input data files for testing.

## Installation

To build the project, first compile the C extension by running the following command:

\`\`\`bash
python setup.py build_ext --inplace
\`\`\`

To compile and run the C program independently:

\`\`\`bash
bash comp.sh
\`\`\`

## Usage

You can run the program by providing the following command-line arguments:

\`\`\`bash
python3 spkmeans.py <k> <goal> <input_file>
\`\`\`

- \`k\`: Number of clusters. If 0, the eigengap heuristic will determine the number of clusters.
- \`goal\`: Task to perform. Options are:
  - \`spk\`: Perform full spectral clustering.
  - \`wam\`: Output the weighted adjacency matrix.
  - \`ddg\`: Output the diagonal degree matrix.
  - \`lnorm\`: Output the normalized graph Laplacian.
  - \`jacobi\`: Output the eigenvalues and eigenvectors.
- \`input_file\`: Path to the input file containing the data points.

Example:

\`\`\`bash
python3 spkmeans.py 3 spk input.txt
\`\`\`

## Output

- When performing spectral clustering (\`spk\`), the program outputs the initial centroids and final centroids obtained from the K-means algorithm.
- For other goals, the corresponding matrix (e.g., the weighted adjacency matrix, diagonal degree matrix) will be printed.

## Testing

You can test the program with custom input files in \`.txt\` or \`.csv\` format. Input files should contain rows of data points, where each row represents a point in \(R^d\).

## Assumptions

- Input data points are distinct and have at most 10 features and 1000 points.
- Output matrices are formatted to four decimal places.

## References

1. Andrew Ng, Michael Jordan, and Yair Weiss. On spectral clustering: Analysis and an algorithm.
Advances in neural information processing systems, 14:849–856, 2001.
2. Ulrike Von Luxburg. A tutorial on spectral clustering. Statistics and computing, 17(4):395–416, 2007.

## Thanks

Special thanks to my collaborator **Offir Dassa**, who contributed equally to this project. It was a great experience working together on this challenging and rewarding project.

## Additional Resources

For reference, this repository also includes a directory titled Preliminary_Exercises, which contains the exercises completed throughout the "Software Project" course at Tel Aviv University. These exercises contributed to building the foundational skills and concepts that influenced the development of the final project.
