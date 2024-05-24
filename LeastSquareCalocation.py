import numpy as np 
import matplotlib.pyplot as plt

class LSCalocation:
    
    def spatial_cov(self,X, Y, G):
        smax = np.sqrt((np.max(X) - np.min(X))**2 + (np.max(Y) - np.min(Y))**2)
        ds = np.sqrt((2 * np.pi * (smax / 2)**2) / len(X))
        max_index = int(np.round(smax / ds)) + 2
        covariance = np.zeros(max_index)
        ncov = np.zeros(max_index)
        
        for i in range(len(G)):
            for j in range(len(G)):
                if j != i:
                    xx = X[i] - X[j]
                    yy = Y[i] - Y[j]
                    r = np.sqrt(xx**2 + yy**2)
                    ir = int(np.round(r / ds)) + 1
                    if r < smax:
                        covariance[ir] += G[i] * G[j]
                        ncov[ir] += 1
        
        for i in range(len(covariance)):
            if ncov[i] != 0:
                covariance[i] /= ncov[i]
        
        covdist = np.arange(len(covariance)) * ds
        return covariance, covdist


    def expcovarfit(self,X, Y, G, covariance, covdist, figplease):
        Dmax = np.sqrt((np.max(X) - np.min(X))**2 + (np.max(Y) - np.min(Y))**2)
        Step = np.sqrt((2 * np.pi * (Dmax / 2)**2) / len(X))
        C0 = np.std(G)**2
        s = covdist
        fiterr = float('inf')
        Dbest = None

        for D in np.arange(Step, Dmax, Step):
            covp = C0 * np.exp(-s / D)
            current_fiterr = np.sum((covp - covariance)**2)
            if current_fiterr < fiterr:
                fiterr = current_fiterr
                Dbest = D
                if figplease == 'covfigure':
                    plt.figure(10)
                    plt.clf()
                    plt.plot(covdist, covariance, 'o', label='Empirical covariance')
                    plt.plot(covdist, covp, 'g', label='Fitted model')
                    plt.xlabel('Distance')
                    plt.ylabel('Covariance')
                    plt.title('Empirical covariance in blue, fitted model in green.')
                    plt.legend()
                    plt.show()

        return C0, Dbest



    def LSCexponential(self,Xi, Yi, X, Y, C0, D, N, G):
        # Initialize the covariance matrices
        Czz = np.zeros((len(X), len(X)))
        Csz = np.zeros((len(Xi), len(X)))
        
        # Compute the distances and covariances
        s2 = (X[:, None] - X[None, :])**2 + (Y[:, None] - Y[None, :])**2
        r = np.sqrt(s2)
        Czz = np.exp(-r / D)
        
        s2i = (Xi[:, None] - X[None, :])**2 + (Yi[:, None] - Y[None, :])**2
        ri = np.sqrt(s2i)
        Csz = np.exp(-ri / D)
        
        # Scale the covariance matrices by C0
        Czz *= C0
        Csz *= C0
        
        # Compute the solution using the provided formulas
        LF = np.linalg.inv(Czz + np.diag(N))
        SolG = Csz @ (LF @ G)
        
        return SolG

