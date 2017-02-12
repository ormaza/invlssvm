# invlssvm

Least Square Support Vector Machine (LS-SVM) is an optimal margin classifier that requires a set of input parameters that, if correctly selected, will cause the machine to achieve good results in its predictions. Cross-Validation (CV) is a technique that assists in this selection, calculating the performance of each set of parameters tested. In spite of its importance, its application requires a high computational cost, making it impossible to execute in more modest hardware. The CV requires the inversion of large data matrices, called kernel matrix. As a solution to this problem, this work will propose the technique of Inversion of Matrices in Blocks. This technique will aid in the optimization of CV, allocating less memory for this calculation, and obtaining similar results. All necessary steps will be presented to adapt the original algorithm, as well as tests in time and memory consumption with the purpose of validating the proposal.

### Requires

##### LS-SVMlab
##### Armadillo
##### Matlab and C-Mex compilers
