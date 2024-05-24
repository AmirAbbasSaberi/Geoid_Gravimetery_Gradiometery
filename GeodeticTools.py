import numpy as np

class GeodeticCalculations:
    def __init__(self, a, e2, f, mm):
        """
        Initialize the GeodeticCalculations class with geodetic parameters.
        
        Parameters:
        a (float): Semi-major axis
        e2 (float): Square of the first eccentricity
        f (float): Flattening
        mm (float): Mass model
        """
        self.a = a
        self.e2 = e2
        self.f = f
        self.mm = mm
    
    def radius_physical(self, phi, landa, h):
        """
        Calculate the physical radius at given geodetic coordinates.
        
        Parameters:
        phi (float): Latitude in degrees
        landa (float): Longitude in degrees
        h (float): Ellipsoidal height
        
        Returns:
        float: Physical radius
        """
        phi = np.radians(phi)
        landa = np.radians(landa)
        N = self.a / np.sqrt(1 - self.e2 * np.sin(phi)**2)
        X = (N + h) * np.cos(phi) * np.cos(landa)
        Y = (N + h) * np.cos(phi) * np.sin(landa)
        Z = (N * (1 - self.e2) + h) * np.sin(phi)
        rp = np.sqrt(X**2 + Y**2 + Z**2)
        return rp
    
    def radius_geoid(self, phi, landa):
        """
        Calculate the geoidal radius at given geodetic coordinates.
        
        Parameters:
        phi (float): Latitude in degrees
        landa (float): Longitude in degrees
        
        Returns:
        float: Geoidal radius
        """
        phi = np.radians(phi)
        landa = np.radians(landa)
        N = self.a / np.sqrt(1 - self.e2 * np.sin(phi)**2)
        X = N * np.cos(phi) * np.cos(landa)
        Y = N * np.cos(phi) * np.sin(landa)
        Z = N * (1 - self.e2) * np.sin(phi)
        rg = np.sqrt(X**2 + Y**2 + Z**2)
        return rg
    
    def gammas_grs80(self, phi, h):
        """
        Calculate the normal gravity potential at given geodetic coordinates.
        
        Parameters:
        phi (float): Latitude in degrees
        h (float): Ellipsoidal height
        
        Returns:
        float: Normal gravity potential
        """
        # Calculate normal gravity potential
        gamma = 9.780327 * (1 + 0.0053024 * np.sin(np.radians(phi))**2 - 0.0000058 * np.sin(2 * np.radians(phi))**2)
        return gamma


