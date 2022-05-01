from Constant import *

import dpdata

def _get_EOS(raw_dir):

    natoms = np.loadtxt(os.path.join(raw_dir, 'type.raw')).shape[0]

    fpar = np.loadtxt(os.path.join(raw_dir, 'fparam.raw'))
    virials = np.loadtxt(os.path.join(raw_dir, 'virial.raw'))
    alg = (virials[:,0] + virials[:,4] + virials[:,-1])/3

    # Angstrom^3
    box = np.loadtxt(os.path.join(raw_dir, 'box.raw'))
    vol = np.linalg.det(box.reshape(-1,3,3))

    therm_term = natoms/vol * kb * fpar * (m2A**3) *1e-9

    press = therm_term + alg/vol  /J2eV * (m2A**3) *1e-9

    energy = np.loadtxt(os.path.join(raw_dir, 'energy.raw'))
    
    return fpar, press, vol/natoms, energy/natoms

class SingleSys():

    def __init__(self, path):
        self.path = path
        self.temps, self.press, self.vol, self.energy = _get_EOS(self.path)
        
    def _show_plot(self, ax, INTERVAL, INIT, color, label):
        #print('Number of Frames: %.d'%self.press.shape[0])#([self.press<=thres][::INTERVAL].shape[0])
        
        ax.plot(self.press[INIT::INTERVAL], self.temps[INIT::INTERVAL], 
                'o', color=color,alpha=0.4, mew=0.5, ms=3, label=label)
   
    def _filter_plot(self, ax, thres, prefix, INTERVAL, color, label):
        
        if prefix == 'press':
            print('Number of Frames: %.d'%self.press[self.press<=thres][::INTERVAL].shape[0])
        
            ax.plot(self.press[self.press<=thres][::INTERVAL], self.temps[self.press<=thres][::INTERVAL],
                    'x',alpha=0.4, mew=0.5, mfc='none',ms=3, label=label)
        elif prefix == 'temps':
            
            print('Number of Frames: %.d'%self.press[self.temps>=thres][::INTERVAL].shape[0])
        
            ax.plot(self.press[self.temps<=thres][::INTERVAL], self.temps[self.temps<=thres][::INTERVAL],
                    'x',alpha=0.4, mew=0.5, mfc='none',ms=3, label=label)    
            
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
        
    def _get_isotherm(self, T):
        
        p = self.press[self.temps==T]
        V = self.vol[self.temps==T]
        E = self.energy[self.temps==T]
        
        return p,V,E        
    
class DPDataSet():
    
    def __init__(self, prefix, SET_DIR, natoms):
        
        self.prefix = prefix
        try:
            self.SET_DIR = SET_DIR
            self.natoms, self.temps, self.press, self.vol, self.energy = self._get_EOS(SET_DIR)
        except:
            self.natoms = natoms
            
    def _get_EOS(self, raw_dir):

        natoms = np.loadtxt(raw_dir + 'type.raw').shape[0]

        fpar = np.loadtxt(raw_dir + 'fparam.raw')
        virials = np.loadtxt(raw_dir + 'virial.raw')
        alg = (virials[:,0] + virials[:,4] + virials[:,-1])/3

        # Angstrom^3
        box = np.loadtxt(raw_dir + 'box.raw')
        vol = np.linalg.det(box.reshape(-1,3,3))

        therm_term = natoms/vol * kb * fpar * (m2A**3) *1e-9

        press = therm_term + alg/vol  /J2eV * (m2A**3) *1e-9

        energy = np.loadtxt(raw_dir + 'energy.raw')

        return natoms, fpar, press, vol, energy

    def _read_energy_test(self):
        self.energy = np.loadtxt(self.prefix+'.e.out')/self.natoms
        
        E_rmse = np.sqrt((self.energy[:,1] - self.energy[:,0] )**2)
        
        self.E_rmse_ave = np.average(E_rmse)
        
        print('RMSE of energy per natom is %.5f eV'%self.E_rmse_ave)
        
        
    def _read_force_test(self):

        force = np.loadtxt(self.prefix+'.f.out')

        self.norm_dft = np.sqrt((force[:,0] )**2 + 
                         (force[:,1] )**2 +
                         (force[:,2] )**2 )
        self.norm_dp = np.sqrt((force[:,3] )**2 + 
                         (force[:,4] )**2 +
                         (force[:,5] )**2 )
        self.f_diff = (force[:,0] - force[:,3])**2 + (force[:,1] - force[:,4])**2 + (force[:,2] - force[:,5])**2 
                
        f_rmse = np.sqrt((force[:,0] - force[:,3])**2 + 
                         (force[:,1] - force[:,4])**2 +
                         (force[:,2] - force[:,5])**2 )

        self.f_rmse_ave = np.average(f_rmse)
        print('RMSE of force is %.5f eV/Angstrom'%self.f_rmse_ave)
        
        conf_ave_f_norm_std = np.std(self.norm_dft.reshape(-1,self.natoms),axis=1)
        conf_ave_f_rmse = np.average(f_rmse.reshape(-1,self.natoms),axis=1)

        self.rrmse = conf_ave_f_rmse/conf_ave_f_norm_std * 100
        print('RMSE/STD of force is %.2f '%np.average(self.rrmse)+'%')
    
    def _read_virial_test(self):
        virials = np.loadtxt(self.prefix+'.v.out') 
        
        self.alg_dft = (virials[:,0] + virials[:,4] + virials[:,8])/3
        self.alg_dp = (virials[:,9] + virials[:,13] + virials[:,17])/3
        try:
            self.p_dft = self.alg_dft/self.vol/J2eV * (m2A**3) *1e-9
            self.p_dp = self.alg_dp/self.vol/J2eV * (m2A**3) *1e-9
            
            self.p_rmse = np.sqrt((self.p_dft-self.p_dp)**2)
            print('RMSE of pressure is %.5f GPa'%np.average(self.p_rmse))
        except:         
            
            self.v_rmse =  np.sqrt((self.alg_dft-self.alg_dp)**2)
            print('RMSE of virial is %.5f eV'%np.average(self.v_rmse))
            
    def _plot_data(self, data_dft, data_dp, INTERVAL):
        
        ax.plot(data_dft[::INTERVAL], data_dp[::INTERVAL],'o',ms=1)
        
        vmin = np.min(data_dft)
        vmax = np.max(data_dft)
        
        x1 = np.linspace(vmin,vmax)
        ax.plot(x1,x1,'--k', lw=1)
        ax.set_xlim(vmin,vmax)
        ax.set_ylim(vmin,vmax)

    def _subplot_data(self, idx, data_dft, data_dp, INTERVAL):
        
        ax[idx].plot(data_dft[::INTERVAL], data_dp[::INTERVAL],'o',ms=1)
        
        vmin = np.min(data_dft)
        vmax = np.max(data_dft)
        
        x1 = np.linspace(vmin,vmax)
        ax[idx].plot(x1,x1,'--k', lw=1)
        ax[idx].set_xlim(vmin,vmax)
        ax[idx].set_ylim(vmin,vmax)
        