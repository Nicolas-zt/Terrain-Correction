import flopy
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from scipy.ndimage import zoom, gaussian_filter
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
import pickle 



# Grid
nlay = 3 # Number of sub-soil layers
nrow = 600 # Number of rows
ncol = 600 # Nuber of columns
delr = delc = 10.0 # Size of rows / columns in meters


# Time periods
nper = 8 # Number of time periods
perlen = [30,60,45,45,60,30,60,35] # Duration of each period in days
nstp = [1]*8 # Number ofsteps per period, determines the temporal precision



# Opening the dtm file as .vrt
try :
    with rasterio.open("../4.Data/MNT/MNT_TOULON.vrt") as src:
        dtm = src.read(1)  
        
except : 
    with rasterio.open("../../4.Data/MNT/MNT_TOULON.vrt") as src:
        dtm = src.read(1)
       
       
def DTM():
    '''
    Creates a 2D plot of a Digital Terrain Model

    '''
    
     
        
        
    plt.figure(figsize=(8, 6))
    plt.contourf(dtm,levels = 15, cmap="terrain")
    plt.colorbar(label="Altitude (m)")
    plt.title("Digital Terrain Model (DTM)")
    plt.gca().invert_yaxis()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()
    

 

def layer_creator(pos,size,shape,thickness,depth = 50):
    '''
    Creates a surface to be used with modflow as a bottom layer. the grid is manually adapted 
    to the one used with the modflow model.
    
    Parameters:
        
        - pos: Center(s) of the layer.
        - size: size of the aquifere.
        - shape: shape of the aquifere (ellipse, gaussian, rectangle or flat).
                 
        - thickness: thickness of the aquifere whithin the layer. Different from the the layer's depth
        
    Out:
        
        - layer: matrix of depths at each location of a grid
    '''
    
    layer = np.zeros((nrow, ncol))
    
    x = np.arange(0.5, ncol) * delr
    y = np.arange(0.5, nrow) * delc
    
    X, Y = np.meshgrid(x, y)
    
    if shape == "ellipse":
        
        cx, cy = pos[0], pos[1]  # centers
        rx, ry = size[0], size[1]  # radius
        
        layer = -thickness*np.exp(-(((X - cx)**2 / rx**2 + (Y - cy)**2 / ry**2))) - depth
        
        
    elif shape == "gaussian":
        
        cx, cy = pos[0], pos[1] 
        sx, sy = size[0], size[1]
        
        layer = thickness * np.exp(-(((X - cx)**2)/(2 * sx**2) + ((Y - cy)**2)/(2 * sy**2))) - depth
          
    elif shape == "rectangle":
            
        
        cx, cy = pos[0], pos[1]
        sx, sy = size[0], size[1]
        
        mask = (np.abs(X - cx) <= sx / 2) & (np.abs(Y - cy) <= sy / 2)
        layer[mask] += thickness
        
        layer -= depth
        
    elif shape == "flat":
        
        layer -= thickness
    
    
    elif shape == "dual":
        
        cx, cy = pos[0], pos[1]
        sx, sy = size[0], size[1]
        
        mask1 = (np.abs(X - cx + 440) <= sx / 2) & (np.abs(Y - cy) <= sy / 2)
        layer[mask1] += thickness
        
        mask2 = (np.abs(X - (cx + 440)) <= sx / 2) & (np.abs(Y - cy) <= sy / 2)
        layer[mask2] += thickness
        
        layer -= depth
    
    return layer 

def layer_plotter(layer):
    '''
    Creates a 2D plot of a layer's surface with matplotlib

    '''
    x = np.arange(0, ncol * delr, delr)  # 0 à 6000 m
    y = np.arange(0, nrow * delc, delc)  # 0 à 6000 m
    X, Y = np.meshgrid(x, y)
    
    plt.figure(figsize=(8,6))
    plt.contourf(X,Y,layer,levels = 100, cmap="terrain")
    plt.colorbar(label="Altitude du fond de nappe (m)")
    plt.gca().invert_yaxis()
    plt.title("Surface basse de la nappe (modélisée)")
    plt.show()
    
def mf_layers():
    ######################################### Layers ###################################################

    # Calcul des facteurs de réduction
    zoom_factors = (nrow / dtm.shape[0], ncol / dtm.shape[1])
    
    # Redimensionnement bilinéaire
    top = zoom(dtm, zoom_factors, order=1)  # order=1 = interpolation bilinéaire
    
    top = layer_creator([0,0],[0,0],'flat',50)
    
    
    ### Ground thickness
    botm1 = top + layer_creator([1000,1000],[2000,1000],'flat',50)  # sol d'épaisseur 50m
    
     
    ### 1st layer - ground
    
    
    botm2 = botm1 + layer_creator([5000,5000],[1000,4000],"flat",50)
    
    botm2b = botm1 + layer_creator([3000,3000],[400,2000],"dual",1)
    
    ### 2nd layer - water 
    
    botm3 = botm2 + layer_creator([1000,4000],[3000,1500],'flat',50)
    
    botm3b = botm2b + layer_creator([3000,3000],[400,2000],"dual",-1)
    
    
    ######################################### Layers ###################################################


    return top,botm1,botm2,botm2b, botm3, botm3b


def modflow(botm):
    '''
    Runs the hydrology model Modflow 6 with different created sub-soils and a real Digital Terrain Model.
    The Modflow 6 version is different from Modflow 2005 which is the previous version.
    The documentation for this model can be found on flopy.readthedocs.io, or by typing the name of a precise function.
    Different parameters can be set up to complexify the simulations.

    '''  
    
    
    
    
    ######################################### Modflow setup ############################################
    # Calcul des facteurs de réduction
    zoom_factors = (nrow / dtm.shape[0], ncol / dtm.shape[1])
    
    # Redimensionnement bilinéaire
    top = zoom(dtm, zoom_factors, order=1)  # order=1 = interpolation bilinéaire
    '''
    The part above is only needed to use a digital terrain mode, otherwise a matrix representing
    a surface is enough
    '''
    
    # top = layer_creator([0,0],[0,0],'flat',50)

    
    
    ### Création du modèle
    modelname = "nappe_perte"
    mf = flopy.mf6.MFSimulation(sim_name=modelname, version="mf6", exe_name="mf6", sim_ws="./modflow6_model")

    

    tsmult = [1.0] * nper
    
    tdis = flopy.mf6.ModflowTdis(
    mf,
    nper=nper,
    time_units="DAYS",
    perioddata=[(perlen[i], nstp[i], tsmult[i]) for i in range(nper)],
    )
    

    
    gwf = flopy.mf6.ModflowGwf(mf, modelname=modelname)
    
    dis = flopy.mf6.ModflowGwfdis(
    gwf,
    nlay=nlay,
    nrow=nrow,
    ncol=ncol,
    delr=delr,
    delc=delc,
    top=top,
    botm=botm,
    )
    
    
    
    
    ### Start conditions
    
    ibound = np.ones((nlay, nrow, ncol), dtype=int)

    
    # Start levels
    strt = 20 * np.ones((nlay, nrow, ncol))

    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)
    
    
    ### Hydrogeological parameters
    hk_values = [ 
       1E-4, # Ground m/s
       1E-1, # 1st layer
       5E-1, # 2nd layer
        ]
        
    vk_values = [val for val in hk_values]
    sy_values = [
        0.15, # Ground
        0.25, # 1st layer
        0.30, # 2nd layer
        ]
    
    hk = np.empty((nlay,nrow,ncol))
    vk = np.empty((nlay,nrow,ncol))
    sy = np.empty((nlay,nrow,ncol))
    
    for i in range(nlay):
        hk[i, :, :] = hk_values[i]
        vk[i, :, :] = vk_values[i]
        sy[i, :, :] = vk_values[i]
    
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=hk, k33=vk)
    
    sto = flopy.mf6.ModflowGwfsto(
    gwf,
    sy=sy,
    ss=1e-5,
    iconvert=1,
    steady_state={0: False},
    transient={0: True}
    )
    
    ### Water variations
    
    # Data
    
    
    perioddata = {
        0 : [(0, "RATE", -0.0),
              # (1, "RATE", 0.0),
              # (2, "RATE", 0.0)
            ],
        1 : [(0, "RATE", -50.0),
               (1, "RATE", -50.0),
               (2, "RATE", -20.0)
            ],
        2 : [(0, "RATE", -30.0),
               (1, "RATE", -30.0),
               (2, "RATE", -10.0)
            ],
        3 : [(0, "RATE", -100.0),
               (1, "RATE", -100.0),
               (2, "RATE", -50.0)
            ],
        4 : [(0, "RATE", -100.0),
               (1, "RATE", -100.0),
               (2, "RATE", -50.0)
            ],
        5 : [(0, "RATE", -30.0),
               (1, "RATE", -30.0),
               (2, "RATE", -10.0) 
            ],
        6 : [(0, "RATE", -50.0),
               (1, "RATE", -50.0),
               (2, "RATE", -20.0)
            ],
        7 : [(0, "RATE", -0.0),
               (1, "RATE", -0.0),
               (2, "RATE", -0.0)
            ],
        }
    
    
    
    '''
    The following section is to use precise points with a loss of water. In the model it is called maws
    but it can represent pumping positions for example. It is very local, the effect is only 
    visible a few meters around.
    
    Packagedata is the list of maws to use. Each have an id, a radius, a bottom and other 
    more specific parameters.
    
    connectiondata is a list with an element per layer crossed by the maw. For example if the maw
    with the id 0 goes through 7 layers, then 7 elements are needed with the positions where 
    the maw crosses the layer. It is equivalent to having 7 maws that need to be connected 
    to each other.
    
    perioddata is a list that contains, for each maw, the volume of water that varies.
    '''
    # maw_pkg = flopy.mf6.ModflowGwfmaw(
    # gwf,
    # pname="MAW",
    # boundnames=True,
    # print_input=True,
    # print_flows=True,
    # save_flows=True,
    # mover=False,
    # auxiliary = ["rate_supp"],
    # packagedata=[
    #     (0, 0.1, botm7[500,200], top[500,200], "SKIN", 7,0.0),  # (maw_id, radius, bottom, type, pump elev, head)
    #     (1, 0.1, botm3[400,500], top[400,500], "SKIN", 2,0.0),
    #     (2, 0.1, botm5[300,300], top[300,300], "SKIN", 4,0.0)
    # ],
    # connectiondata=[
    #     (0, 0, (0, 500, 200), top[500,200]  , botm1[500,200], hk_values[0] - 1E-7, 0.2),  # 1st layer
    #     (0, 1, (1, 500, 200), botm1[500,200], botm2[500,200], hk_values[1] - 1E-7, 0.2),  # 2nd layer
    #     (0, 2, (2, 500, 200), botm2[500,200], botm3[500,200], hk_values[2] - 1E-7, 0.2),  # 3rd layer
    #     (0, 3, (3, 500, 200), botm3[500,200], botm4[500,200], hk_values[3] - 1E-7, 0.2),  # 4th layer
    #     (0, 4, (4, 500, 200), botm4[500,200], botm5[500,200], hk_values[4] - 1E-7, 0.2),  # 5th layer
    #     (0, 5, (5, 500, 200), botm5[500,200], botm6[500,200], hk_values[5] - 1E-7, 0.2),  # 6th layer
    #     (0, 6, (6, 500, 200), botm6[500,200], botm7[500,200], hk_values[6] - 1E-7, 0.2),  # 7th layer
    #     (1, 0, (0, 400,500), top[400,500]  , botm1[400,500], hk_values[0] - 1E-7, 0.2),
    #     (1, 1, (1, 400,500), botm1[400,500], botm2[400,500], hk_values[1] - 1E-7, 0.2),
    #     (2, 0, (0, 300, 300), top[300,300]  , botm1[300,300], hk_values[0] - 1E-7, 0.2),
    #     (2, 1, (1, 300, 300), botm1[300,300], botm2[300,300], hk_values[1] - 1E-7, 0.2),
    #     (2, 2, (2, 300, 300), botm2[300,300], botm3[300,300], hk_values[2] - 1E-7, 0.2),
    #     (2, 3, (3, 300, 300), botm3[300,300], botm4[300,300], hk_values[3] - 1E-7, 0.2)    
    # ],
    # perioddata=perioddata
    # )
    
   
    # Reload
    
    '''
    The reload is a volume of water that is added at each period (it is a matrix) it can represent rain or drought periods
    The way it is now, it is added on the whole area
    '''
    
    rch = flopy.mf6.modflow.ModflowGwfrch(
    gwf,
    stress_period_data={
    0: [(0, i, j, 0.0001) for i in range(nrow) for j in range(ncol)],
    1: [(0, i, j, -0.0007) for i in range(nrow) for j in range(ncol)],
    2: [(0, i, j, -0.0001) for i in range(nrow) for j in range(ncol)],
    3: [(0, i, j, 0.0007) for i in range(nrow) for j in range(ncol)],
    4: [(0, i, j, 0.0007) for i in range(nrow) for j in range(ncol)],
    5: [(0, i, j, -0.0007) for i in range(nrow) for j in range(ncol)],
    6: [(0, i, j, -0.0007) for i in range(nrow) for j in range(ncol)],
    7: [(0, i, j, -0.00001) for i in range(nrow) for j in range(ncol)],
},  
    pname="RCH-1",
    filename="model.rch"
    )




    ### Simulation
    
    ims = flopy.mf6.ModflowIms(
    mf,
    print_option='SUMMARY',
    complexity='SIMPLE',
    outer_dvclose=1e-2,
    outer_maximum=100,
    under_relaxation='NONE',
    inner_maximum=100,
    inner_dvclose=1e-4,
    rcloserecord=1e-3,
    linear_acceleration='BICGSTAB'
)
    
    oc = flopy.mf6.ModflowGwfoc(
    gwf,
    head_filerecord=f"{modelname}.hds",
    budget_filerecord=f"{modelname}.cbc",
    saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]
    )
    
    
    ######################################### Modflow setup ############################################
    
    
    
    ### Exécution
    mf.write_simulation()
    success, buff = mf.run_simulation()
    
    ### Lecture des résultats
    # headobj = MF6HeadFile("nappe_perte.gwf.nappe_perte.hds")  # ou adapte le chemin si tu as changé le nom du modèle
    # headobj = flopy.utils.HeadFile('./modflow6_model/nappe_perte.hds') 
    # times = headobj.get_times()
    # head = headobj.get_data(totim=times[-1])
    
    
    
    # fig = plt.figure(figsize=(10, 8))
    # ax = fig.add_subplot(111, projection='3d')
    
    # X, Y = np.meshgrid(np.arange(ncol), np.arange(nrow))
    
    # # Exemple : top et botm[1] (première couche)
    # ax.plot_surface(X, Y, top, cmap='terrain', alpha=0.7)
    # ax.plot_surface(X, Y, botm[0], cmap='Greens', alpha=0.5)
    # ax.plot_surface(X, Y, botm[1], cmap='Blues', alpha=0.4)
    
    # ax.set_title("Représentation 3D des couches du modèle")
    # ax.set_xlabel("Colonne")
    # ax.set_ylabel("Ligne")
    # ax.set_zlabel("Altitude (m)")
    # plt.show()
    
    
    # Ouvrir le fichier de têtes
    # Lire les flux du fichier de budget
    # cbc = gwf.output.budget()
    
    # Chercher les flux du package MAW
    # maw_flows = cbc.get_data(text='MAW')
    
    # # Afficher les débits à chaque pas de temps
    # for i, data in enumerate(maw_flows):
    #     print(f"--- Période {i} ---")
    #     for record in data:
    #         print(record)
            
    # head = gwf.output.head().get_data()
    # print("Niveau d'eau dans les cellules connectées en 500,200 :", head[:, 500, 200])
    # print("Niveau d'eau dans les cellules connectées en 400,500 :", head[:, 400,500])
    # print("Niveau d'eau dans les cellules connectées en 300,300 :", head[:, 300, 300])

        
        
    idomain = gwf.modelgrid.idomain
    print(np.unique(idomain))
    
    head = gwf.output.head().get_data()
    print("Min head:", np.nanmin(head), "Max head:", np.nanmax(head))
    print("Nombre de cellules sèches :", np.isnan(head).sum())

    
    headfile = flopy.utils.HeadFile('./modflow6_model/nappe_perte.hds')  # adapter chemin si besoin
    times = gwf.output.head().get_times()  # Liste des temps totaux
    
    # Gets the average head at each time step
    periods_head = []
    avg_heads = []
    anomaly_heads = []
    for per in range(nper):
        for kstp in range(nstp[0]):
            head = headfile.get_data(kstpkper=(kstp, per))
            periods_head.append(head)
            head[head < -1e+20] = np.nan  # Cellules sèches
            anomaly_heads.append(head[2,300,300])
            avg_heads.append(np.nanmean(head))
    
    
            ### Piezometric levels visualisation
            plt.figure(figsize=(8, 6))
            plt.imshow(head[2], cmap="Blues", extent=(0, ncol * delr, 0, nrow * delc), origin ='lower')
            plt.colorbar(label="Niveau piézométrique (m)")
            plt.title("Gravi simulation")
            plt.xlabel("Distance Est (m)")
            plt.ylabel("Distance Nord (m)")
            plt.gca().invert_yaxis()
            plt.grid(False)
            plt.show()
    
    # Tracer
    plt.figure(figsize=(8, 4))
    plt.plot(times, avg_heads, marker='o')
    plt.xlabel("Temps (jours)")
    plt.ylabel("Niveau d'eau moyen (m)")
    plt.title("Average water level variations")
    plt.grid(True)
    plt.ylim(min(avg_heads) - 0.01, max(avg_heads) + 0.01)
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(8, 4))
    plt.plot(times, anomaly_heads, marker='o')
    plt.xlabel("Temps (jours)")
    plt.ylabel("Niveau d'eau moyen (m)")
    plt.title("Average water level variations in 3000 3000")
    plt.grid(True)
    plt.ylim(min(avg_heads) - 0.01, max(avg_heads) + 0.01)
    plt.tight_layout()
    plt.show()
    
    print( times)
        
    
    '''
    ################################### Area 3D visualisation ##########################################
    
    # Define the grid's size to be the same as the modflow's one
    x = np.arange(0, ncol * delr, delr)
    y = np.arange(0, nrow * delc, delc)
    X, Y = np.meshgrid(x, y)
    
    plotter = pv.Plotter()
    colors = ['lightgreen', 'blue', 'cyan', 'orange', 'blue','orange','blue','orange']
    
    for i in range(len(botm)):
        # Determine the top and bottom surfaces
        z_top = top if i == 0 else botm[i - 1]
        z_botm = botm[i]
    
        # Creates the 3D points for top and bottom surface of each layer
        z1 = z_top.ravel()
        z2 = z_botm.ravel()
        points = np.concatenate([
            np.c_[Y.ravel(), X.ravel(), z1],
            np.c_[Y.ravel(), X.ravel(), z2]
        ])
    
        # 3D grid
        grid = pv.StructuredGrid()
        grid.points = points
        grid.dimensions = (ncol, nrow, 2)
    
        # Adds volume to the plot
        plotter.add_mesh(grid, color=colors[i % len(colors)], opacity=1, show_edges=False)
    
    plotter.add_axes()
    plotter.show()
    

    ################################### Area 3D visualisation ##########################################
    '''
    return periods_head , times




if __name__ == "__main__":
    
   # layer_plotter(layer_creator([2000,1000],[1800,1000],"ellipse",5))
   # layer_plotter(layer_creator([3000,3000],[1000,4000],'rectangle',5))
   # layer_plotter(layer_creator([0,0],[0,0],'flat',0,5))
   # layer_plotter(layer_creator([3000,3000],[1000,750],'gaussian',-5))
   # layer_plotter(layer_creator([0,0],[0,0],'flat',0,5))
   # layer_plotter(layer_creator([4000,3000],[3000,3000],'ellipse',5))
   # layer_plotter(layer_creator([0,0],[0,0],'flat',0,5))
   # top,botm1,botm2,botm3,botm4,botm5,botm6,botm7,botm8 = mf_layers()
   # modflow()    
   # test()  
   # DTM()
   pass