import numpy as np 
import random as rd
import matplotlib.pyplot as plt

       
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
                
    mnt_filtre = grid*masque
    
    return mnt_filtre

if __name__ == "__main__":
    
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