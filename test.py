import numpy as np 
import random as rd
import matplotlib.pyplot as plt

class volume:
    '''
    Définit une classe d'objets 3D servant de volumes unitaires dans le calcul
    d'effet gravimétrique de terrain.
    '''
    def __init__(self,density):
        
        self.density = density
        
        
    
class prism(volume):
    ''' 
    Type de volume simple ppermettant de calculer des effets gravimétriques de
    terrain.
    
    -shape : taille planimétrique du prisme
    -height : hauteur du prisme récupérée dans le MNT
    '''
    
    def __init__(self,shape,height,density):
        
        super().__init__(density)
        self.shape = shape
        self.size = self.shape[0]*self.shape[1]
        self.height = height
        
    def attraction_gravi(self,r,dz):
        ''' 
        Calcul l'attarction gravitationnelle d'un prisme en fonction d'une 
        distance et d'une dénivelée
        '''
        
        volume = self.size * self.height
        mass = volume * self.density
        
        g = G * mass /r**2 * dz/r # dz/r --> Direction component 
        
        return g
        
        
class Grid:
    
    ''' 
    Définit la grille représentant l'espace dans lequel on observe le champ
    de pesanteur.
    
    - X : largeur de la grille
    - Y : longueur de la grille
    - stations : liste des stations gravimétriques à traiter
    '''
    
    def __init__(self,length,width,stations):
        
        self.length = length
        self.width = width
        self.stations = stations
    
       
def masque_adap(grid,seuil):
    
    n = grid.shape[0]
    masque = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            # distance par rapport à la station
            d = np.sqrt((i - S[1])**2 + (j - S[0])**2)
            if  d <= seuil :
                masque[i,j] = 1
                
            else:
                masque[i, j] = 1/d
                
    
    
    return masque


            
            
    

if __name__ == "__main__":
    
    G = 6.67430e-11   # Constante gravitationnelle
    
    n = 20

    X_grid, Y_grid = np.meshgrid(
    np.arange(0,n,1),  
    np.arange(0,n,1)  
    )

    ### Génération MNT aléatoire
    random_mnt = np.zeros((n,n))
    for i in range(len(random_mnt)):
        for j in range(len(random_mnt)): 
        
            random_mnt[i,j] = rd.random()*100
            
    ### Position de la station


    S = (10,10)

    ### Masque autour de la station
    rayon = 10

    MNT = masque_adap(random_mnt,7)      
                
    
     
    plt.contourf(X_grid,Y_grid,MNT, cmap = 'terrain')
    # plt.contourf(masque, cmap='Reds', alpha=0.4)
    plt.scatter(S[0],S[1],s = 60,color = 'red')
    
    plt.colorbar( label='Alti')
    plt.xticks(np.arange(0, n, 0.5))
    plt.yticks(np.arange(0, n, 0.5))
    plt.grid(True, color='black', linewidth=0.5)
    plt.show()