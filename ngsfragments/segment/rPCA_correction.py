import numpy as np
from scipy.linalg import pinv, svd
from sklearn.utils.extmath import randomized_svd
import warnings


class RRPCAResult:
    """
    Result object for Randomized Robust PCA decomposition.
    
    Attributes:
    -----------
    L : numpy.ndarray
        Low-rank component matrix
    S : numpy.ndarray  
        Sparse component matrix
    k : int
        Final estimated rank
    err : list
        List of Frobenius norm errors at each iteration
    """
    def __init__(self):
        self.L = None
        self.S = None
        self.k = None
        self.err = []

def rrpca(A, lambda_param=None, maxiter=100, tol=1.0e-4, p=10, q=2, trace=False, rand=True):
    """
    Randomized robust principal component analysis (rrpca).
    
    Robust principal components analysis separates a matrix into a low-rank plus sparse component.
    
    Robust principal component analysis (RPCA) is a method for the robust separation of a
    rectangular (m,n) matrix A into a low-rank component L and a sparse component S:
    
    A = L + S
    
    To decompose the matrix, we use the inexact augmented Lagrange multiplier
    method (IALM). The algorithm can be used in combination with either the randomized or deterministic SVD.
    
    Parameters:
    -----------
    A : numpy.ndarray
        A real (m, n) input matrix to be decomposed.
    lambda_param : float, optional
        Tuning parameter (default lambda = max(m,n)^-0.5).
    maxiter : int, optional
        Maximum number of iterations (default maxiter = 100).
    tol : float, optional
        Precision parameter (default tol = 1.0e-4).
    p : int, optional
        Oversampling parameter for randomized SVD (default p=10).
    q : int, optional
        Number of additional power iterations for randomized SVD (default q=2).
    trace : bool, optional
        Print progress.
    rand : bool, optional
        If True, the randomized SVD routine is used, otherwise standard SVD is used.
    
    Returns:
    --------
    RRPCAResult
        Object containing:
        - L: low-rank component; (m, n) dimensional array
        - S: sparse component; (m, n) dimensional array
        - k: final estimated rank
        - err: list of errors at each iteration
        
    References:
    -----------
    [1] Lin, Zhouchen, Minming Chen, and Yi Ma.
        "The augmented lagrange multiplier method for exact
        recovery of corrupted low-rank matrices." (2010).
        (available at arXiv http://arxiv.org/abs/1009.5055).
        
    [2] N. B. Erichson, S. Voronin, S. L. Brunton and J. N. Kutz. 2019.
        Randomized Matrix Decompositions Using R.
        Journal of Statistical Software, 89(11), 1-48.
    """
    if trace:
        print("This is Python version of rrpca")
    
    # Convert to numpy array and get dimensions
    A = np.asarray(A, dtype=float)
    m, n = A.shape
    
    # Initialize result object
    rrpca_obj = RRPCAResult()
    
    # Set target rank
    k = 1
    if k > min(m, n):
        k = min(m, n)
    
    # Handle missing values (set NaN to 0)
    A = np.nan_to_num(A, nan=0.0)
    
    # Set lambda, gamma, rho
    if lambda_param is None:
        lambda_param = max(m, n) ** -0.5
    gamma = 1.25
    rho = 1.5
    
    # Set SVD algorithm
    svdalg = 'rsvd' if rand else 'svd'
    
    # Compute matrix norms
    if svdalg == 'svd':
        spectral_norm = np.linalg.norm(A, ord=2)
    elif svdalg == 'rsvd':
        # Use randomized SVD to estimate spectral norm
        try:
            _, s, _ = randomized_svd(A, n_components=1, n_oversamples=10, n_iter=1)
            spectral_norm = s[0]
        except:
            # Fallback to standard norm if randomized SVD fails
            spectral_norm = np.linalg.norm(A, ord=2)
    else:
        raise ValueError("Selected SVD algorithm is not supported!")
    
    inf_norm = np.linalg.norm(A, ord=np.inf) / lambda_param
    dual_norm = max(spectral_norm, inf_norm)
    fro_norm = np.linalg.norm(A, ord='fro')
    
    # Initialize Lagrange multiplier
    Z = A / dual_norm
    
    # Initialize tuning parameter
    mu = gamma / spectral_norm
    mubar = mu * 1e7
    mu = min(mu * rho, mubar)
    muinv = 1.0 / mu
    
    # Initialize low-rank and sparse matrix
    L = np.zeros((m, n))
    S = np.zeros((m, n))
    
    niter = 1
    err = 1
    
    while err > tol and niter <= maxiter:
        
        # Update S using soft-threshold
        epsi = lambda_param / mu
        temp_S = A - L + Z / mu
        
        S = np.zeros((m, n))
        
        # Soft thresholding
        idx_L = temp_S < -epsi
        idx_H = temp_S > epsi
        S[idx_L] = temp_S[idx_L] + epsi
        S[idx_H] = temp_S[idx_H] - epsi
        
        # Singular Value Decomposition
        R = A - S + Z / mu
        
        if svdalg == 'svd':
            U, s, Vt = svd(R, full_matrices=False)
        elif svdalg == 'rsvd':
            # Use randomized SVD with adaptive strategy
            if k > min(m, n) / 5:
                auto_svd = 'svd'
            else:
                auto_svd = 'rsvd'
            
            if auto_svd == 'svd':
                U, s, Vt = svd(R, full_matrices=False)
            else:
                try:
                    n_components = min(k + 10, min(m, n))
                    U, s, Vt = randomized_svd(R, n_components=n_components, 
                                            n_oversamples=p, n_iter=q)
                except:
                    # Fallback to standard SVD
                    U, s, Vt = svd(R, full_matrices=False)
        
        # Predict optimal rank and update
        svp = np.sum(s > 1.0 / mu)
        
        if svp <= k:
            k = min(svp + 1, n)
        else:
            k = min(svp + round(0.05 * n), n)
        
        # Truncate SVD and update L
        if svp > 0:
            # Ensure we don't exceed available singular values
            svp = min(svp, len(s))
            # Soft threshold the singular values
            s_thresh = s[:svp] - 1.0 / mu
            # Keep only positive values
            pos_idx = s_thresh > 0
            if np.any(pos_idx):
                s_thresh = s_thresh[pos_idx]
                U_thresh = U[:, :svp][:, pos_idx]
                Vt_thresh = Vt[:svp, :][pos_idx, :]
                L = U_thresh @ np.diag(s_thresh) @ Vt_thresh
            else:
                L = np.zeros((m, n))
        else:
            L = np.zeros((m, n))
        
        # Compute error
        Astar = A - L - S
        Z = Z + Astar * mu
        
        err = np.linalg.norm(Astar, ord='fro') / fro_norm
        rrpca_obj.err.append(err)
        
        if trace:
            print(f'\nIteration: {niter}, predicted rank = {svp}, target rank k = {k}, Fro. error = {err:.6f}')
        
        # Update mu
        mu = min(mu * rho, mubar)
        muinv = 1.0 / mu
        
        niter += 1
    
    # Set final results
    rrpca_obj.L = L
    rrpca_obj.S = S
    rrpca_obj.k = k
    
    return rrpca_obj


def thresh(x, mu):
    """
    Soft-thresholding function for sparse regularization.
    
    This is the second phase in preparing the Panel of Normals.
    Internal shrinkage function for values in S vector:
    y = sgn(x)max(|x| - mu, 0)
    
    Parameters:
    -----------
    x : numpy.ndarray
        Numeric vector of length == number of markers in the sample. 
        Vector to which shrinkage is to be applied.
    mu : float
        Shrinkage parameter to use (default = lambda2 = 1/(sqrt(no. of markers in input vector)))
    
    Returns:
    --------
    numpy.ndarray
        Regularized vector y
    
    Author: Aditya Deshpande (translated to Python)
    """
    x = np.asarray(x)
    y = np.maximum(x - mu, 0)
    y = y + np.minimum(x + mu, 0)
    return y

def apg_project(m_vec, U, lambda1, lambda2):
    """
    Accelerated Proximal Gradient based projection method.
    
    Project new sample into the burnin space solving projection 
    by Accelerated Proximal Gradient.
    
    Parameters:
    -----------
    m_vec : numpy.ndarray
        Numeric vector of length m. Vector of GC corrected coverage 
        data of sample in question
    U : numpy.ndarray
        (m markers x n samples) numeric matrix. The basis of low rank subspace. 
        The dimensions are same as burnin matrix
    lambda1 : float
        Tuning parameter (default = 1/(sqrt(no. of markers in input vector))
    lambda2 : float  
        Tuning parameter (default = 1/(sqrt(no. of markers in input vector))
    
    Returns:
    --------
    list
        List containing [v, s] vectors where:
        - v: coefficient vector in the subspace
        - s: sparse residual vector
    
    Author: Aditya Deshpande (translated to Python)
    """
    # Convert inputs to numpy arrays
    m_vec = np.asarray(m_vec).reshape(-1, 1)
    U = np.asarray(U)
    
    q, p = U.shape
    v = np.zeros((p, 1))
    s = np.zeros((q, 1))
    I = np.eye(p)
    
    converged = False
    k = 0
    maxiter = 200
    
    # Pre-compute the projection matrix
    UUt = pinv(U.T @ U + lambda1 * I) @ U.T
    
    while not converged:
        k += 1
        v_old = v.copy()
        v = UUt @ (m_vec - s)
        
        s_old = s.copy()
        s = thresh(m_vec - (U @ v), lambda2).reshape(-1, 1)
        
        # Compute convergence criterion
        v_diff = np.linalg.norm(v - v_old, 'fro')
        s_diff = np.linalg.norm(s - s_old, 'fro')
        e = max(v_diff, s_diff) / q
        
        if e < 1e-6 or k > maxiter:
            converged = True
    
    return [v, s]

def wash_cycle(m_vec, L_burnin, S_burnin, r, U_hat, V_hat, sigma_hat, 
               lambda1=None, lambda2=None, verbose=True):
    """
    Function to perform online rPCA decomposition.
    
    Function that calls the online rPCA on sample under question and updates subspace
    basis and does decomposition.
    
    Parameters:
    -----------
    m_vec : numpy.ndarray
        Numeric vector of length m. Vector of GC corrected coverage data of sample in question
    L_burnin : numpy.ndarray
        (m markers x n samples) numeric matrix. L matrix of panel of normals after batch rPCA decomposition
    S_burnin : numpy.ndarray  
        (m markers x n samples) numeric matrix. S matrix of panel of normals after batch rPCA decomposition
    r : int
        Estimated rank of panel of normals after batch rPCA decomposition
    U_hat : numpy.ndarray
        (m markers x n samples) numeric matrix. Right singular matrix of L_burnin
    V_hat : numpy.ndarray
        (m markers x n samples) numeric matrix. Left singular matrix of L_burnin
    sigma_hat : numpy.ndarray
        Numeric vector of length r. Singular values of L_burnin
    lambda1 : float, optional
        Tuning parameter (default = 1/(sqrt(no. of markers in input vector))
    lambda2 : float, optional
        Tuning parameter (default = 1/(sqrt(no. of markers in input vector))
    verbose : bool, optional
        Outputs progress (default=True)
    
    Returns:
    --------
    list
        List containing [L_vec, si] where:
        - L_vec: L vector for sample in question
        - si: S vector for sample in question
    
    Author: Aditya Deshpande (translated to Python)
    """
    if verbose:
        print("Using the detergent provided to start washing")
    
    # Convert inputs to numpy arrays
    m_vec = np.asarray(m_vec)
    L_burnin = np.asarray(L_burnin)
    S_burnin = np.asarray(S_burnin)
    U_hat = np.asarray(U_hat)
    V_hat = np.asarray(V_hat.copy())  # Make a copy since we'll modify it
    sigma_hat = np.asarray(sigma_hat)
    
    m, n = L_burnin.shape
    
    # Set default lambda values if not provided
    if lambda1 is None:
        lambda1 = 1.0 / np.sqrt(m)
        lambda2 = 1.0 / np.sqrt(m)
        if verbose:
            print("lambdas calculated")
    
    # Compute U from U_hat and sigma_hat
    U = U_hat @ np.diag(np.sqrt(sigma_hat))
    
    # Commented out matrix calculations that aren't used in final implementation
    # A = np.zeros((r, r))
    # B = np.zeros((m, r))
    
    if verbose:
        print("calculating A and B")
    
    # Original R code had a loop to calculate A and B matrices but it's commented out
    # for i in range(V_hat.shape[1]):
    #     A += np.outer(V_hat[:, i], V_hat[:, i])
    #     B += np.outer((L_burnin[:, i] - S_burnin[:, i]), V_hat[:, i])
    
    if verbose:
        print("calculating v and s")
    
    # Project the sample using APG
    projection = apg_project(m_vec, U, lambda1, lambda2)
    vi = projection[0].flatten()  # Convert to 1D array
    si = projection[1].flatten()  # Convert to 1D array
    
    # Update V_hat matrix
    vi_subtract = V_hat[:, 0]  # First column to be replaced
    V_hat = np.column_stack([V_hat, vi])  # Add new column
    
    # Commented out online update calculations
    # A = A + np.outer(vi, vi) - np.outer(vi_subtract, vi_subtract)  
    # B = B + np.outer((m_vec - si), vi) - np.outer((L_burnin[:, 0] - S_burnin[:, 0]), vi_subtract)
    
    if verbose:
        print("Calculating b")
    
    # Commented out U update step
    # U = update_cols(U, A, B, lambda1)  # This function isn't defined
    
    # Compute final L vector
    L_vec = U @ vi
    
    return [L_vec, si]