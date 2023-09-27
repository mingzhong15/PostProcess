from Constant import *

import dpdata

def _get_EOS(raw_dir):
    
    natoms = np.loadtxt(os.path.join(raw_dir, 'type.raw')).shape[0]

    try:
        energy = np.loadtxt(os.path.join(raw_dir, 'energy.raw'))
        
        virials = np.loadtxt(os.path.join(raw_dir, 'virial.raw'))
        alg = (virials[:,0] + virials[:,4] + virials[:,-1])/3

        # Angstrom^3
        box = np.loadtxt(os.path.join(raw_dir, 'box.raw'))
        vol = np.linalg.det(box.reshape(-1,3,3))

        press =  alg/vol  /J2eV * (m2A**3) *1e-9
        try:
            fpar = np.loadtxt(os.path.join(raw_dir, 'fparam.raw'))
            therm_term = natoms/vol * kb * fpar * (m2A**3) *1e-9

            press += therm_term
        except:
            fpar = np.ones(energy.shape) * 0
        
        
    except:
        energy = np.load(os.path.join(raw_dir,'set.000', 'energy.npy'))
    
        virials = np.load(os.path.join(raw_dir, 'set.000', 'virial.npy'))
        alg = (virials[:,0] + virials[:,4] + virials[:,-1])/3

        # Angstrom^3
        box = np.load(os.path.join(raw_dir,'set.000',  'box.npy'))
        vol = np.linalg.det(box.reshape(-1,3,3))

        press =  alg/vol  /J2eV * (m2A**3) *1e-9
        try:
            fpar = np.loadtxt(os.path.join(raw_dir,'set.000', 'fparam.npy'))
            therm_term = natoms/vol * kb * fpar * (m2A**3) *1e-9

            press += therm_term
        except:
            fpar = np.ones(energy.shape) * 0
        
        
    return fpar, press, vol/natoms, energy/natoms, natoms

class SingleSys():

    def __init__(self, path):
        self.path = path
        self.temps, self.press, self.vol, self.energy, self.natoms = _get_EOS(self.path)
        
    def _read_force(self):
        
        self.force = np.loadtxt(os.path.join(self.path, 'force.raw'))
    
    def _show_plot(self, ax, INTERVAL, INIT, color, label):
        #print('Number of Frames: %.d'%self.press.shape[0])#([self.press<=thres][::INTERVAL].shape[0])
        
        ax.plot(self.press[INIT::INTERVAL], self.temps[INIT::INTERVAL], 
                'o', color=color,alpha=0.4, mew=0.5, ms=3, mfc='none',label=label)
   
    def _filter_plot(self, ax, thres, prefix, INTERVAL, color, label):
        
        if prefix == 'press':
            print('Number of Frames: %.d'%self.press[self.press<=thres][::INTERVAL].shape[0])
        
            ax.plot(self.press[self.press<=thres][::INTERVAL], self.temps[self.press<=thres][::INTERVAL],
                    'o', color=color,alpha=0.4, mew=0.5, mfc='none',ms=3, label=label)
        elif prefix == 'temps':
            
            print('Number of Frames: %.d'%self.press[self.temps<=thres][::INTERVAL].shape[0])
        
            ax.plot(self.press[self.temps<=thres][::INTERVAL], self.temps[self.temps<=thres][::INTERVAL],
                    'o', color=color,alpha=0.4, mew=0.5, mfc='none',ms=3, label=label)    
 
    def cri_filter_all(self, cri, INTERVAL, OUT_DIR):
        
        if not os.path.exists(OUT_DIR):
            os.mkdir(OUT_DIR)
            
        for prefix in ['box','coord','energy','force','fparam','virial']:
            CAT_FILE= os.path.join(self.path, prefix+'.raw') 
            OUT_FILE= os.path.join(OUT_DIR, prefix+ '.raw')

            data = np.loadtxt(CAT_FILE)

            filter_data = data[cri][::INTERVAL]

            print(data.shape,' -->',filter_data.shape)

            np.savetxt(OUT_FILE, filter_data)
            print(OUT_FILE + ' is done \n')
            
        os.chdir(self.path)
        os.system('cp type.raw type_map.raw '+ OUT_DIR)
        
        os.chdir(OUT_DIR)
        os.system('~/raw_to_set.sh 50000')
                    
    def _filter_all(self, thres, filter_prefix, INTERVAL, OUT_DIR):
        
        if not os.path.exists(OUT_DIR):
            os.mkdir(OUT_DIR)
            
        for prefix in ['box','coord','energy','force','fparam','virial']:
            CAT_FILE= os.path.join(self.path, prefix+'.raw') 
            OUT_FILE= os.path.join(OUT_DIR, prefix+ '.raw')

            data = np.loadtxt(CAT_FILE)

            if filter_prefix == 'press':
                filter_data = data[self.press<=thres][::INTERVAL]
        
            elif filter_prefix == 'temps':
                filter_data = data[self.temps<=thres][::INTERVAL]   

            print(data.shape,' -->',filter_data.shape)

            np.savetxt(OUT_FILE, filter_data)
            print(OUT_FILE + ' is done \n')
            
        os.chdir(self.path)
        os.system('cp type.raw type_map.raw '+ OUT_DIR)
        
        os.chdir(OUT_DIR)
        os.system('~/raw_to_set.sh 50000')
        
    def _get_isotherm(self, T):
        
        p = self.press[self.temps==T]
        V = self.vol[self.temps==T]
        E = self.energy[self.temps==T]
        
        return p,V,E        
