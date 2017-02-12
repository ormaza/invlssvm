# invlssvm

The Least Squares Support Vector Machine (LSSVM) requires, as part of its cross-validation algorithm, the calculation of the inverse of the kernel matrix. Although the whole inverse is computed only blocks of elements close to its main diagonal are used in the cross-validation. We present a solution for the computation of the needed inverse blocks as a way to decrease the memory footprint of the algorithm.

### Requires

##### LS-SVMlab
##### Armadillo
##### Matlab and C-Mex compilers

### Documents

##### TCC FINAL (pt-br)
##### ijcnn article (en-us)
