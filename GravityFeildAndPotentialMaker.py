import numpy as np

class GravityCalculator:
    def __init__(self,Do):

        self.Do = Do
    def Tn(self,lmax, anm, bnm, landa, rg, a, plm):
        """
        Calculate the sum of anomaly Potential up to degree lmax.
        
        Parameters:
        lmax (int): Maximum degree
        Clm (numpy.ndarray): Cosine coefficients
        Slm (numpy.ndarray): Sine coefficients
        lambda_ (float): Longitude in degrees
        rg (float): Geoidal radius
        a (float): Semi-major axis
        Plm (numpy.ndarray): Associated Legendre functions
        
        Returns:
        float: Sum of Potential
        """
        Cosnm = np.zeros((lmax,lmax))
        Sinnm = np.zeros((lmax,lmax))
        Rr = np.zeros((lmax,lmax))
        for i in range(3,lmax+1):
            m = np.linspace(1,i,i)-1
            Cosnm[i-1,0:len(m)] = np.cos(np.radians(m*landa))
            Sinnm[i-1,0:len(m)] = np.sin(np.radians(m*landa))
            Rr[i-1,0:len(m)] = ((a/rg)**(i))*np.ones((1,len(m)))
        Ylm = anm*Cosnm + bnm*Sinnm
        Tn = np.sum(Rr*Ylm*plm)

        return Tn

    def Dg_n(self,lmax, anm, bnm, landa, rg, a, plm):
        """
        Calculate the sum of Gravity anomaly up to degree lmax.
        
        Parameters:
        lmax (int): Maximum degree
        Clm (numpy.ndarray): Cosine coefficients
        Slm (numpy.ndarray): Sine coefficients
        lambda_ (float): Longitude in degrees
        rg (float): Geoidal radius
        a (float): Semi-major axis
        Plm (numpy.ndarray): Associated Legendre functions
        
        Returns:
        float: Sum of Delta g_n 
        """
        Cosnm = np.zeros((lmax,lmax))
        Sinnm = np.zeros((lmax,lmax))
        Rr = np.zeros((lmax,lmax))
        for i in range(3,lmax+1):
            m = np.linspace(1,i,i)-1
            Cosnm[i-1,0:len(m)] = np.cos(np.radians(m*landa))
            Sinnm[i-1,0:len(m)] = np.sin(np.radians(m*landa))
            Rr[i-1,0:len(m)] = ((a/rg)**(i))*((i-2)/rg)*np.ones((1,len(m)))
        Ylm = anm*Cosnm + bnm*Sinnm
        Tn = np.sum(Rr*Ylm*plm)

        return Tn

    def Dg_d(self,lmax, anm, bnm, landa, rp, a, plm):
        """
        Calculate the sum of Gravity Disturbance up to degree lmax.
        
        Parameters:
        lmax (int): Maximum degree
        Clm (numpy.ndarray): Cosine coefficients
        Slm (numpy.ndarray): Sine coefficients
        lambda_ (float): Longitude in degrees
        rg (float): Geoidal radius
        a (float): Semi-major axis
        Plm (numpy.ndarray): Associated Legendre functions
        
        Returns:
        float: Sum of Delta g_d 
        """
        Cosnm = np.zeros((lmax,lmax))
        Sinnm = np.zeros((lmax,lmax))
        Rr = np.zeros((lmax,lmax))
        for i in range(3,lmax+1):
            m = np.linspace(1,i,i)-1
            Cosnm[i-1,0:len(m)] = np.cos(np.radians(m*landa))
            Sinnm[i-1,0:len(m)] = np.sin(np.radians(m*landa))
            Rr[i-1,0:len(m)] = ((a/rp)**(i))*((i)/rp)*np.ones((1,len(m)))
        Ylm = anm*Cosnm + bnm*Sinnm
        Tn = np.sum(Rr*Ylm*plm)

        return Tn


    def T_rr(self,lmax, anm, bnm, landa, rp, a, plm):
        """
        Calculate the sum of Second deritive of potential up to degree lmax.
        
        Parameters:
        lmax (int): Maximum degree
        Clm (numpy.ndarray): Cosine coefficients
        Slm (numpy.ndarray): Sine coefficients
        lambda_ (float): Longitude in degrees
        rg (float): Geoidal radius
        a (float): Semi-major axis
        Plm (numpy.ndarray): Associated Legendre functions
        
        Returns:
        float: Sum of Trr
        """
        Cosnm = np.zeros((lmax,lmax))
        Sinnm = np.zeros((lmax,lmax))
        Rr = np.zeros((lmax,lmax))
        for i in range(3,lmax+1):
            m = np.linspace(1,i,i)-1
            Cosnm[i-1,0:len(m)] = np.cos(np.radians(m*landa))
            Sinnm[i-1,0:len(m)] = np.sin(np.radians(m*landa))
            Rr[i-1,0:len(m)] = ((a/rp)**(i-1))*((i)*(i-1))*np.ones((1,len(m)))
        Ylm = anm*Cosnm + bnm*Sinnm
        Tn = np.sum(Rr*Ylm*plm)

        return Tn