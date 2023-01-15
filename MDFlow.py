from Constant import *

class MDSys():
    
    def __init__(self, traj_dir):
        
        self.traj_dir = traj_dir
        self.thermo = os.path.join(self.traj_dir , 'thermo.dat')
        
        self.rdf = os.path.join(self.traj_dir , 'rdf.txt')
        if os.path.exists(self.rdf):
            try:
                self._read_rdf()
            except:
                pass
        
        self.log = os.path.join(self.traj_dir , 'log.run' )
        
        self.natoms = int(os.popen("grep ' atoms' "+ self.log).readlines()[0].split()[1])
        
    def _read_thermo(self, is_printf=True):
        
        data = np.loadtxt(self.thermo)
        
        ave_data = np.average(data, axis=0)
        
        self.temps = ave_data[1]
        self.press = ave_data[2] * bar2Pa/1e9
        self.vol = ave_data[3] / self.natoms
        self.density = ave_data[4]
        self.energy = ave_data[5] / self.natoms
    
        if is_printf:
            print('%.2f %.2f %.4f %.4f %.4f'%(self.temps, self.press, self.vol, self.density, self.energy))
    
    def _read_rdf(self):
        
        self.rdf_data = np.loadtxt(self.rdf, skiprows=4)
        
    def _get_coordination(self):
        
        delta_r = self.rdf_data[1,1] - self.rdf_data[0,1]
        
        try:
            self.coord = 4 * np.pi / self.vol *  np.cumsum( self.rdf_data[:,1] **2 * self.rdf_data[:,2] * delta_r)
        except:
            self._read_thermo(is_printf=False)
            
            self.coord = 4 * np.pi / self.vol *  np.cumsum( self.rdf_data[:,1] **2 * self.rdf_data[:,2] * delta_r)
            
        
def plot_RDF(data):
    
    ax.plot(data[:,1],data[:,2])
    ax.set_xlim(0,6)
    ax.set_ylim(0,)        
    ax.set_xlabel('Distance ($\\rm{\mathring A}$)')
    ax.set_ylabel('$g(r)$')
