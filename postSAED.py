from Constant import *
from Constant import radial_average_2d, chunk_average

class SingleSAED():
    
    def __init__(self, FILE, is_print=True):

        file = os.popen('grep DIMENSIONS '+FILE)
        tt = file.readlines()[0].split()
        Nx = int(tt[1])
        Ny = int(tt[2])
        Nz = int(tt[3])
        if is_print:
            print(Nx,Ny,Nz)   
        self.L = np.array([Nx,Ny,Nz])

        self.axis_ave = np.where(self.L==1)[0][0]
        self.x1_idx = np.where(self.L!=1)[0][0]
        self.x2_idx = np.where(self.L!=1)[0][1]
        
        file = os.popen('grep ASPECT_RATIO '+FILE)
        tt = file.readlines()[0].split()
        dx = float(tt[1])
        dy = float(tt[2])
        dz = float(tt[3])
        self.dL = np.array([ dx, dy, dz])

        file = os.popen('grep ORIGIN '+FILE)
        tt = file.readlines()[0].split()
        x0 = float(tt[1])
        y0 = float(tt[2])
        z0 = float(tt[3])
        self.L0 = np.array([x0, y0, z0])
        
        raw_data = np.loadtxt(FILE,skiprows=10)
        self.data = raw_data.reshape(self.L)
        
        
    def _get_2d_intensity(self):
    
        self.intensity = np.average(self.data, axis=self.axis_ave)

        grid_x1 = np.arange( self.L[self.x1_idx] )* self.dL[self.x1_idx] + self.L0[self.x1_idx]
        grid_x2 = np.arange( self.L[self.x2_idx] )* self.dL[self.x2_idx] + self.L0[self.x2_idx]
        self.X1, self.X2 = np.meshgrid(grid_x1, grid_x2)


    def _get_3d_intensity(self):

        grid_x = np.arange( self.L[0] )* self.dL[0] + self.L0[0]
        grid_y = np.arange( self.L[1] )* self.dL[1] + self.L0[1]
        grid_z = np.arange( self.L[2] )* self.dL[2] + self.L0[2]
             
        temp = np.meshgrid(grid_x,grid_y,grid_z)
        self.x = temp[0].reshape(-1,1)
        self.y = temp[1].reshape(-1,1)
        self.z = temp[2].reshape(-1,1)
        self.intensity = self.data.reshape(-1,1)
  
    def _get_radius_intensity(self, init_r = 0.2, N=200):
        
        self._get_2d_intensity()
        
        Z = np.array(self.intensity)
        R = ( np.sqrt(self.X1**2 + self.X2**2) )

        self.radius = np.linspace(init_r, np.max(R), N)
        self.int_r = np.zeros(N)
        
        for i in range(N-1):
    
            r = self.radius[i]
            r_dr = self.radius[i+1]
            dr = r_dr - r

            left= (R>r)
            right = (R<=r_dr)
            cri = (left & right)

            self.int_r[i] = np.sum(Z[cri])/(2*np.pi*r*dr)
    
#     def _get_radius_intensity(self, is_cut = False, DELTA = 10, is_thres = False, thres = 1e3):
    
#         self._get_2d_intensity()
        
#         if is_cut:
#             mid = int(self.intensity.shape[0]/2)+1
#             self.intensity[mid-DELTA:mid+DELTA,:] = 0
#             self.intensity[:, mid-DELTA:mid+DELTA] = 0
        
#         R = ( np.sqrt(self.X1**2 + self.X2**2) ).reshape(-1,1)
        
#         temp_int = self.intensity.reshape(-1,1)
#         temp_int[temp_int == -1] = 0
        
#         if is_thres:
#             temp_int[temp_int > thres] = thres
        
#         self.radius, self.int_r = radial_average_2d(R, temp_int )
        

# def plot_2d_ele_frac(self, ax, thres, kmax):

#     self._get_2d_intensity()

#     Z = np.array(self.intensity)
#     Z[Z>thres] = thres
#     Z = np.exp(np.log(Z))
#     cs = ax.contourf(self.X1, self.X2, Z, 100, vmin=0, vmax=thres,
#                         cmap=plt.cm.viridis)

#     plt.colorbar(cs,label='Intensity')

#     ax.set_xlim(-kmax,kmax)
#     ax.set_ylim(-kmax,kmax)

#     ax.set_xlabel('y ($\\rm{\mathring A^{-1}}$)')
#     ax.set_xlabel('z ($\\rm{\mathring A^{-1}}$)')