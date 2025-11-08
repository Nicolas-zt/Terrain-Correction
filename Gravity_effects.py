import numpy as np 
import matplotlib.pyplot as plt
from Test_area.Area import *
import os 
from scipy.interpolate import griddata , LinearNDInterpolator
from scipy.interpolate import RegularGridInterpolator



    

def signal(head,  cell_size, botm, mode):
    
    '''
    Calculates a grid of gravity signal. Each cell's gravity contribution is added 
    for each measurement position.The gravimetric signal of each cell is calculated and then added.
    It is not the 
    
    Parameters : 
        
        - head : grid of piezometric levels
        - cell_size : size of each cell
        - botm : list of every bottom layer
        
    Returns :
        
        - total_gravity : signal grid
        - X_obs , Y_obs : measurement positions

    '''
    
    
    G = 6.67430e-11   # Constante gravitationnelle
    rho = 1000        # Water density
    rho_ground = 2700 # Ground density
    spatial_step = 5


    X_full, Y_full = np.meshgrid(
    np.arange(0.5 * delc, (ncol + 0.5) * delc + 1, delc),  # The 0.5 is only there to have measurements in the center of the cells
    np.arange(0.5 * delr, (nrow + 0.5) * delr + 1, delr)   
    )
    
    
    
    X_obs = X_full[::spatial_step, ::spatial_step]      # Measuring positions under sampling -> 1 pt every spatial_step m
    Y_obs = Y_full[::spatial_step, ::spatial_step]
    
    
   
    
    
    
    total_gravity = np.zeros_like(X_obs, dtype=float) # Matrix initialization

    for row in range(nrow):         # Double loop on the measured cells 
        for col in range(ncol):
            
            
            for i in range(nlay) :  # Loop on the layers
                
                
                # Cell's center of mass
                x = (col + 0.5) * cell_size
                y = (row + 0.5) * cell_size
                  
    
                # Horizontal distance between cell and sensor
                dx = x - X_obs
                dy = y - Y_obs
                

                ### Mass calculation   
                if i == 0   :                                       # Ground
                    height = top[row,col] - botm[i][row,col]
                    volume = height * cell_size ** 2
                    mass = rho_ground * volume
                    
                elif i == 1 :                                       # Considered as ground continuity
                
                    height = botm[i-1][row,col] - botm[i][row,col]   # Distance from top layer
                    
                    volume = height * cell_size ** 2  
                    mass = rho_ground * volume  
                        
                elif i == 2 :
                    
                    if mode == 'spatial':
                        
                        height = botm[i-1][row,col] - botm[i][row,col]
                        
                    if mode == 'temporal':
                        
                        height = botm[i-1][row,col] - head[i][row,col]  # The head is the water level in the cavity 
                    
                    volume = height * cell_size ** 2
                    mass = rho * volume
                                           
                ### Distance calculation
                z = botm[i][row,col] + height / 2   # Vertical middle of the cell
                dz = z - top[row,col] # Vertical distance from sensor
                r = np.sqrt(dx**2 + dy**2 + dz**2) 

                ### Gravity calculation
                g =  abs( G * mass * dz  / r**3) # gravity in m/s², dz/r --> Direction component 
                total_gravity += g
                    
    ### Signal map
    dir_path = '../5.Results'
    os.makedirs(dir_path, exist_ok=True)  
    
    plt.figure(figsize=(8, 6))
    plt.contourf(X_obs,Y_obs,1e8 * total_gravity, levels = 50, cmap = 'viridis')
    plt.colorbar( label='Gravity (µgal)')
    plt.title("Map of gravity signal")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.gca().invert_yaxis()
    plt.legend()
    plt.axis('equal')
    plt.grid(False)
    
    fig_path = dir_path + '/signal.svg'
    plt.savefig(fig_path)
    
    plt.show()
    
    
    return total_gravity , X_obs , Y_obs

def mass(head,botm,mode):
    
    
    G = 6.67430e-11  # gravitational constant
    rho = 1000       # Densité de l'eau (kg/m³)
    rho_ground = 2700
    cell_size = 10 
    
    
    
    Mass_map = np.empty((nrow,ncol))
    
    for row in range(nrow):
        for col in range(ncol):
            # Masse d'eau dans la cellule
            
            mass = 0
            
            for i in range(nlay):
            ### Mass calculation   
                if i == 0   :                                       # Ground
                    height = top[row,col] - botm[i][row,col]
                    volume = height * cell_size ** 2
                    mass = rho_ground * volume
                    
                elif i == 1 :                                       # Considered as ground continuity
                
    
                    height = botm[i-1][row,col] - botm[i][row,col]   # Distance from top layer
                    
                    volume = height * cell_size ** 2  
                    mass = rho_ground * volume  
                        
                elif i == 2 :
                    
                    if mode == 'spatial':
                        
                        height = botm[i-1][row,col] - botm[i][row,col]
                        
                    if mode == 'temporal':
                        
                        height = botm[i-1][row,col] - head[i][row,col]  # The head is the water level in the cavity 
                    
                    volume = height * cell_size ** 2
                    mass = rho * volume
                           

                Mass_map[row,col] = mass
    '''        
    plt.figure(figsize=(8, 6))
    plt.contourf(Mass_map ,levels = 10, cmap="viridis")
    plt.colorbar(label = 'mass')
    plt.title("Mass map")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.gca().invert_yaxis()
    plt.axis('equal')
    plt.grid(False)
    plt.show()       
    '''
def generate_pink_noise(n,shape, amplitude=1e-8):
    """
    Génère un bruit rose (1/f) de longueur n en utilisant un filtre fréquentiel.
    """
    # Génère du bruit blanc
    white_noise = np.random.normal(size=(n,*shape))

    # Transformée de Fourier
    f = np.fft.rfft(white_noise,axis = 0)
    freqs = np.fft.rfftfreq(n)

    # Évite la division par zéro à f = 0
    freqs[0] = 1e-6

    # Atténuation 1/f (bruit rose)
    f /= np.sqrt(freqs[:,np.newaxis,np.newaxis])

    # Retour au domaine temporel
    pink = np.fft.irfft(f, n=n, axis = 0)

    # Normalisation à l'amplitude souhaitée
    pink *= amplitude / np.std(pink, axis=0, keepdims=True)

    return pink





def classical_gravimeter_noise(shape,n_days = 365,timestep = 8,time = True):
    
    '''
    Generates a grid of noise corresponding to a classical gravimeter.
    
    Parameters :
        
        - shape : shape of the grid on which noise is generated

    Returns :
        
        - noise : grid of noise

    '''
    
    if time == True:
        
        t = np.arange(n_days)
        
        # drift = 4.6E-8 / n_days * t
        drift = 0
        
        
        noise = generate_pink_noise(n_days,shape, amplitude = 4E-8)
        noise -= np.mean(noise, axis=0, keepdims=True)
        
        measured_signal = noise 
    
        # Sample one value per period
        # cumulative_days = np.cumsum([0] + perlen)
        # sampled_signal = np.array([
        #     measured_signal[cumulative_days[i], :, :] for i in range (len(perlen))
        # ])
        
        sampled_signal = measured_signal[::timestep, :, :]
        
        return sampled_signal
    
    else :
        
        noise = np.random.normal(0, 1E-8 * 3,size = shape) # White noise ~3 µGal
    
        return noise






def quantum_gravimeter_noise(shape):
    
    '''
    Generates a grid of noise corresponding to a quantum gravimeter.
    For this sensor, a white noise is sufficient.
    
    Parameters :
        
        - shape : shape of the grid on which noise is generated

    Returns :
        
        - noise : grid of noise

    '''
    
    
    noise = np.random.normal(0, 1E-8 * 0.3, size = shape)  # White noise ~0.3 µGal   
    return noise 





def gravity_map(heads,botm):
    
    ''' 
    Calculates a grid of gravity measurements.
    
    Parameters : 
        
        - botm : list of every bottom layer
        
    Returns : 
        
        - classical_measurements : grid of measurements simulated for a classical gravimeter
        - quantum_measurements : grid of measurements simulated for a quantum gravimeter
        - X_obs , Y_obs : measurement positions
    '''
   

    # periods_head, times = modflow(botm)
    
    # Initilization 
    classical_measurements = []
    quantum_measurements = []
 
    
    
    
    
    g , X_obs , Y_obs = signal(heads, delr, botm, mode = 'spatial') # Getting signal and grid's shape


    ### Noise
    
    classical_noise = classical_gravimeter_noise(X_obs.shape,time = False)
    quantum_noise = quantum_gravimeter_noise(X_obs.shape)


    ### Measurements
    
    classical_measurement = 1E8 * (g + classical_noise)
    quantum_measurement = 1E8 * (g + quantum_noise)
    
    
    ### Appends
    
    classical_measurements.append(classical_measurement)
    quantum_measurements.append(quantum_measurement)
    

    '''
    ### Classical visualisation 
    plt.figure(figsize=(8, 6))
    plt.contourf([X_obs,Y_obs],classical_measurements,levels = 50, cmap = 'viridis')
    # plt.scatter(positions[:,0] ,positions[:,1]  , c = classical_values , cmap='hot', label='Mesures')
    plt.colorbar( label='Gravité (m/s²)')
    plt.title("Map of classical gravimetry measurements")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.gca().invert_yaxis()
    plt.legend()
    plt.axis('equal')
    plt.grid(False)
    plt.show()
        
        
    ### Quantum visualisation
    plt.figure(figsize=(8, 6))
    plt.contourf([X_obs , Y_obs], quantum_measurements,levels = 50, cmap = 'viridis')
    plt.colorbar(label='Gravity (m/s²)')
    plt.title("Map of quantum gravimetry measurements")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.gca().invert_yaxis()
    plt.axis('equal')
    plt.grid(False)
    plt.show()
    '''
    
    return classical_measurement, quantum_measurement, X_obs, Y_obs

def time_series(botm):
    
    '''
    Calculates time series of gravity measurements at each measurement position.
    
    Parameters : 
        
        - botm : list of every bottom layer
    '''
    
    n_days = 365
    timestep_quantum = 28  # jours → plus fréquent
    timestep_classical = 28 # jours → moins fréquent
    
    period = 90 # période du signal en jours
    
    # --------------------------------------
    # High resolution signal (quantum)
    # --------------------------------------
    times_quantum = np.arange(0, n_days, timestep_quantum)
    coeff_quantum = 5 + 0.1 * np.cos(2 * np.pi * times_quantum / period)
    periods_head_quantum = [np.ones((3, 600, 600)) * c for c in coeff_quantum]
    
    # --------------------------------------
    # Low resolution signal (classical)
    # --------------------------------------
    times_classical = np.arange(0, n_days, timestep_classical)
    coeff_classical = 5 + 0.1 * np.cos(2 * np.pi * times_classical / period)
    periods_head_classical = [np.ones((3, 600, 600)) * c for c in coeff_classical]
    
    # --------------------------------------
    # Observation positions
    # --------------------------------------
    obs_step = 500
    positions = [(row, col) for row in range(0, nrow, obs_step)
                            for col in range(0, ncol, obs_step)]
    
    x = [col * delc for (row, col) in positions]
    y = [row * delr for (row, col) in positions]
    
    # --------------------------------------
    # Directory creation
    # --------------------------------------
    base_dir = '../5.Results/1.Time_Series/34.14pts_0.1cos_1year_2'
    os.makedirs(base_dir, exist_ok=True)
    
    # --------------------------------------
    # Principal loop
    # --------------------------------------
    for (row, col) in positions:
        i = row // obs_step
        j = col // obs_step
    
        # Classical noise
        X_obs = signal(periods_head_classical[0], cell_size=delr, botm=botm,mode = 'temporal')[1]
        full_classical_noise = classical_gravimeter_noise(X_obs.shape, n_days, timestep_classical)
    
        classical_time_serie = []
        quantum_time_serie = []
        true_g_classical = []
        true_g_quantum = []
        classical_time_noise = []
        quantum_time_noise = []
    
        # --- CLASSICAL ---
        for t_c, head in enumerate(periods_head_classical):
            g_grid = signal(head, cell_size=delr, botm=botm,mode = 'temporal')[0]
            g = g_grid[i, j]
            noise = full_classical_noise[t_c, i, j]
            classical_time_noise.append(noise)
            classical_time_serie.append(g + noise)
            true_g_classical.append(g)
    
        # --- QUANTUM ---
        for t_q, head in enumerate(periods_head_quantum):
            g_grid = signal(head, cell_size=delr, botm=botm,mode = 'temporal')[0]
            g = g_grid[i, j]
            noise = quantum_gravimeter_noise(X_obs.shape)[i, j]  # optionnel
            quantum_time_noise.append(noise)
            quantum_time_serie.append(g + noise)
            true_g_quantum.append(g)

        
        
        ### Arrays + Units conversion
        
        # Noise
        classical_time_noise = 1E8 * np.array(classical_time_noise)
        quantum_time_noise = 1E8 * np.array(quantum_time_noise)
        
        # Signal
        true_g_quantum = 1E8 * np.array(true_g_quantum)
        
        # Measurements
        classical_time_serie = 1E8 * np.array(classical_time_serie)
        quantum_time_serie = 1E8 * np.array(quantum_time_serie)
        
        # SnR
        SnR_classical_time = np.abs(true_g_classical) / np.abs(classical_time_noise)
        SnR_quantum_time   = np.abs(true_g_quantum) / np.abs(quantum_time_noise)
        
        
        
        # Point directory creation 
        
        pt_dir = base_dir + '/' + str(col) + '_' + str(row)
        os.makedirs(pt_dir, exist_ok = True)
        
        # Measurement positions display on the map
        plt.figure(figsize=(8, 6)) 
        plt.imshow(top, cmap="terrain", extent=(0, ncol * delr, 0, nrow * delc))
        plt.scatter(x, y, s=30, c='red')
        plt.xlabel("X (m)")
        plt.ylabel("Y (m)")
        plt.title("Measurement positions")
        plt.grid(False)
        plt.show()
        
        ### --- SIGNAL ---
        plt.figure(figsize=(8, 4))
        plt.plot(times_quantum, true_g_quantum, marker='o', label="Signal")
        plt.xlabel("Time (days)")
        plt.ylabel("Simulated gravity (µGal)")
        plt.title(f"Simulated gravity over time (signal) at ({col*10},{row*10})")
        plt.grid(True)
        plt.legend()
        fig_path = pt_dir + f'/{col*10}_{row*10}_Signal.svg'
        plt.savefig(fig_path)
        plt.show()
        
        
        ### --- NOISE ---
        plt.figure(figsize=(8, 4))
        plt.plot(times_classical, classical_time_noise, marker='o', label='Classical noise')
        plt.plot(times_quantum, quantum_time_noise, marker='o', label='Quantum noise')
        plt.xlabel("Time (days)")
        plt.ylabel("Noise (µGal)")
        plt.title(f"Noise over time at ({col*10},{row*10})")
        plt.grid(True)
        plt.legend()
        fig_path = pt_dir + f'/{col*10}_{row*10}_Noise.svg'
        plt.savefig(fig_path)
        plt.show()
        
        
        ### --- MEASUREMENTS ---
        plt.figure(figsize=(8, 4))
        plt.plot(times_quantum, true_g_quantum,marker = 'o', label='Signal', color='blue', linewidth=2)
        plt.plot(times_classical, classical_time_serie,marker =  'o', label='Classical measurement', color='orange')
        plt.plot(times_quantum, quantum_time_serie, marker = 'o', label='Quantum measurement', color='green')
        plt.xlabel("Time (days)")
        plt.ylabel("Simulated gravity (µGal)")
        plt.title(f"Simulated measurements over time at ({col*10},{row*10})")
        plt.grid(True)
        plt.legend()
        fig_path = pt_dir + f'/{col*10}_{row*10}_Measurements.svg'
        plt.savefig(fig_path)
        plt.show()
        
        
        ### --- SNR ---
        
        plt.figure(figsize=(8, 4))
        plt.plot(times_classical, SnR_classical_time, marker = 'o', label='Classical SNR')
        plt.plot(times_quantum, SnR_quantum_time,marker = 'o', label='Quantum SNR')
        plt.xlabel("Time (days)")
        plt.ylabel("Signal-to-Noise Ratio")
        plt.title(f"SNR over time at ({col*10},{row*10})")
        plt.grid(True)
        plt.legend()
        fig_path = pt_dir + f'/{col*10}_{row*10}_SNR.svg'
        plt.savefig(fig_path)
        plt.show()



 
    
    
def anomaly():
    '''
    Calculates a reference signal and a perturbated signal to make an anomaly map per period of time 
    and per sensor.
    The number of measurement points is defined in signal().
    
    '''
    
    
    head_ref = np.zeros((3,600,600))
    head = np.ones((3,600,600))
    
    x = np.array([i/2 for i in range(15)])
    heads = 0.5 + 0.5*np.cos(x) 
    
    thickness = [0.1,0.25,0.35,0.5,0.75,1]
     
    
    for thick in thickness:
        
        botm2b = botm1 + layer_creator([3000,3000],[400,2000],"dual",thick)
        botm3b = botm2b + layer_creator([3000,3000],[400,2000],"dual",-thick) 
        
        botm2 = botm1 + layer_creator([5000,5000],[1000,4000],"flat")
        botm3 = botm2 + layer_creator([5000,5000],[1000,4000],"flat",50)
        '''
        ref_classical, ref_quantum, X_obs, Y_obs = gravity_map(head_ref   ,[botm1,botm2,botm3]) # Measurements without any anomaly
        test_classical, test_quantum = gravity_map(head,[botm1,botm2b,botm3b])[0:2] # Measurements with anomaly
        '''
        ref_perfect,X_obs, Y_obs = signal(head_ref,10,[botm1,botm2,botm3],'spatial')
        test_perfect = signal(head,10,[botm1,botm2b,botm3b],'spatial')[0]
        
        base_dir = '../5.Results/2.Spatial_Detections/Perfect_solution'
        os.makedirs(base_dir, exist_ok = True)
        
        pt_dir = base_dir + '/Signal_variation'
        os.makedirs(pt_dir, exist_ok = True)
        
        '''
        anomaly_classical =  test_classical - ref_classical
        anomaly_quantum = test_quantum - ref_quantum
        '''
        anomaly_perfect = test_perfect - ref_perfect
         
        '''
        plt.figure(figsize=(8, 6))
        plt.contourf(X_obs, Y_obs, anomaly_classical , levels=100, cmap='seismic')  # µGal
        plt.colorbar(label='Gravimetric anomaly (µGal)')
        # plt.scatter(X_obs,Y_obs, c = 'g', s = 0.5)
        plt.title('Map of classical gravimetric anomaly')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.gca().invert_yaxis()
        plt.axis('equal')
        
        
        fig_path = pt_dir + '/classical_' + str(50 + depth) + '.svg'
        plt.savefig(fig_path)
        
        plt.show()
        
        plt.figure(figsize=(8, 6))
        plt.contourf(X_obs, Y_obs, anomaly_quantum , levels=100, cmap='seismic')  # µGal
        plt.colorbar(label='Gravimetric anomaly (µGal)')
        # plt.scatter(X_obs,Y_obs, c = 'g', s = 0.5)
        plt.title('Map of quantum gravimetric anomaly')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.gca().invert_yaxis()
        plt.axis('equal')
        
        
        fig_path = pt_dir + '/quantum_' + str(50 + depth) + '.svg'
        plt.savefig(fig_path)
        
        plt.show()
        '''
        
        plt.figure(figsize=(8, 6))
        plt.contourf(X_obs, Y_obs, anomaly_perfect , levels=100, cmap='seismic')  # µGal
        plt.colorbar(label='Gravimetric anomaly (µGal)')
        # plt.scatter(X_obs,Y_obs, c = 'g', s = 0.5)
        plt.title('Map of perfect gravimetric anomaly')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.gca().invert_yaxis()
        plt.axis('equal')
        
        
        fig_path = pt_dir + '/perfect' + str(50 + depth) + '.svg'
        plt.savefig(fig_path)
        
        plt.show()
        
if __name__ == "__main__":
    
    top,botm1, botm2,botm2b ,botm3 ,botm3b = mf_layers()
    # gravity_map(0,[botm1,botm2b])
    time_series([botm1,botm2,botm3])
    # anomaly()
    # mass(1,[botm1,botm2,botm3],'spatial')
    