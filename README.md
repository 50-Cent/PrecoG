# PrecoG

# PrecoG

PrecoG is a framework to compute preconditioning matrix for input data in online fashion in order to speed up the convergence of LMS filter. PrecoG is extended to solving linear systems of equations and unsupervised Hebb-LMS problems. Infact, learning paradigms that use the gradient descent approach as an intermediate step to update model weights as a function of online data will be substantially benefited from our work. 

mainProcess.m : a matlab file where generators of several processes (autoregressive (first order and second order), moving average) are enlisted and it yields errors, filter weights over time as outputs after applying transforms, such as DCT, DST and PrecoG.

computeLMS.m : is called from mainProcess.m and it is the canonical LMS filter algorithm implemented for size 3 and 4 of filter taps.

computePrecoGLMS.m : is called from mainProcess.m and calls the subroutine findPrecogTransform.m. It computes the LMS filter weights and error over time.

computeHebbLMS1D.m : is called from mainprocess.m. It computes the canonical Hebb-LMS filter and weights over time.

computePrecoGHebbLMS1D.m : is called from mainProcess.m and calls the subroutine findPrecogTransform.m. It computes the Hebb-LMS filter weights and error over time.


findPrecogTransform.m : This file finds the unitary preconditioning matrix, PrecoG, by updating graph weights using the gradient descent approach. It accepts arguments : R (estimated or asymtotic autocorrelation matrix or any square matrix), A (initialization of graph adjacency matrix), epsn1 (upper bound of maximum eigenvalue), epsn2 (lower bound of maximum eigenvalue), l2coeff (L2 regularization coefficient), learning_rate, maxIter, log_penalty (stoping edge weights from becoming negative). 




Relevant papers:

>> Widrow, Bernard, and Marcian E. Hoff. Adaptive switching circuits. Stanford Univ Ca Stanford Electronics Labs, 1960.

>> Widrow, Bernard, et al. "Nature's learning rule: The Hebbian-LMS algorithm." Artificial Intelligence in the Age of Neural Networks and Brain Computing. Academic Press, 2019. 1-30.

>> Beaufays, Francoise. "Transform-domain adaptive filters: An analytical approach." IEEE Transactions on Signal processing 43.2 (1995): 422-431.


