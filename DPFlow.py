from Constant import *
from Dataset import _get_EOS
from Constant import parity_plot
from Constant import _get_RMSE

class DPSys():
    
    def __init__(self, dp_dir):
        
        self.dir = dp_dir
        
        
        self.ener_dft = []
        self.ener_dp  = []
        self.norm_dft = []
        self.norm_dp  = []
        self.f_diff   = []
        
        self.alg_dft  = []
        self.alg_dp   = []        
        self.p_dft  = []
        self.p_dp   = []

        
    def _plot_lc(self, ax, label, is_vali = False):
        
        try:
            lc_file = os.path.join(self.dir, 'lc.train.out')
            data = np.loadtxt(lc_file)
        except:
            lc_file = os.path.join(self.dir, 'lcurve.out')
            data = np.loadtxt(lc_file)
            
        ll_list = ['Energy (eV)','Force (eV/$\\rm{\mathring A}$)','Virial (eV)']
        
        if not is_vali:

            for idx in range(3):
                if idx == 0:
                    ax[idx].plot(data[:,0],data[:,2+idx],'o',ms=1,mew=0.5, 
                                 label=label)
                else:
                    ax[idx].plot(data[:,0],data[:,2+idx],'o',ms=1,mew=0.5)

                ax[idx].set_xscale('log')
                ax[idx].set_yscale('log')
                ax[idx].set_title(ll_list[idx])
                ax[idx].set_xlabel('stopbatch')
    
    def _read_dp_test(self, prefix, SET_DIR,  is_print=False):

        self.prefix = prefix
        try:
            self.SET_DIR = SET_DIR
            self.temps, self.press, self.vol, self.dft_energy, self.natoms = _get_EOS(SET_DIR)
        except:
            pass
        
        self._update_natoms(is_print)
        
        self._read_energy_test(is_print)
        self._read_force_test(is_print)
        self._read_virial_test(is_print) 

    def _update_natoms(self, is_print):
        
        energy = np.loadtxt( os.path.join(self.dir, self.prefix+'.e.out') )
        
        force = np.loadtxt( os.path.join(self.dir, self.prefix+'.f.out'))
        
        nframes = energy.shape[0]
        self.natoms = int(force.shape[0]/nframes)
        
        if is_print:
            print('Total atoms in this phase is %.d '%self.natoms)
        
    def _read_energy_test(self, is_print):
        energy = np.loadtxt( os.path.join(self.dir, self.prefix+'.e.out') )/self.natoms
        
        self.ener_dft = np.append(self.ener_dft, energy[:,0])
        self.ener_dp  = np.append(self.ener_dp,  energy[:,1])
        
        E_rmse = np.sqrt((self.ener_dft - self.ener_dp )**2)
        
        self.E_rmse_ave = np.average(E_rmse)
        
        if is_print:
            print('RMSE of energy per natom is %.5f meV'%(self.E_rmse_ave*1e3))
        
    def _read_force_test(self, is_print):

        force = np.loadtxt( os.path.join(self.dir, self.prefix+'.f.out'))
        
        self.norm_dft = np.append(self.norm_dft, np.sqrt((force[:,0] )**2 + 
                         (force[:,1] )**2 +
                         (force[:,2] )**2 ))
        
        self.norm_dp = np.append(self.norm_dp, np.sqrt((force[:,3] )**2 + 
                         (force[:,4] )**2 +
                         (force[:,5] )**2 ))
        
        self.f_diff = np.append(self.f_diff,
        (force[:,0] - force[:,3])**2 + (force[:,1] - force[:,4])**2 + (force[:,2] - force[:,5])**2) 
                
        f_rmse = np.sqrt(self.f_diff )

        self.f_rmse_ave = np.average(f_rmse)
        
        if is_print:
            print('RMSE of force is %.5f eV/Angstrom'%self.f_rmse_ave)
        
        #conf_ave_f_norm_std = np.std(self.norm_dft.reshape(-1,self.natoms),axis=1)
        #conf_ave_f_rmse = np.average(f_rmse.reshape(-1,self.natoms),axis=1)

        #self.rrmse = conf_ave_f_rmse/conf_ave_f_norm_std * 100
        self.rrmse = f_rmse/self.norm_dft * 100
        if is_print:
            print('RMSE/NORM of force is %.2f '%np.average(self.rrmse)+'%')
    
    def _read_virial_test(self, is_print):
        virials = np.loadtxt(os.path.join(self.dir, self.prefix+'.v.out')) 
        
        alg_dft = (virials[:,0] + virials[:,4] + virials[:,8])/3
        alg_dp  = (virials[:,9] + virials[:,13] + virials[:,17])/3
        
        self.alg_dft = np.append(self.alg_dft, alg_dft)
        self.alg_dp = np.append(self.alg_dp,   alg_dp)
        
        try:
            self.p_dft = np.append(self.p_dft,
                                 alg_dft/(self.vol*self.natoms)/J2eV * (m2A**3) *1e-9)
            self.p_dp = np.append(self.p_dp,
                                 alg_dp/(self.vol*self.natoms)/J2eV * (m2A**3) *1e-9)

            self.p_rmse = np.sqrt((self.p_dft-self.p_dp)**2)
            self.p_rmse_ave = np.average(self.p_rmse)
            if is_print:
                print('RMSE of pressure is %.5f GPa'%self.p_rmse_ave)
        except:         
            
            self.v_rmse =  np.sqrt((self.alg_dft-self.alg_dp)**2)
            self.v_rmse_ave = np.average(self.v_rmse)
            if is_print:
                print('RMSE of virial is %.5f eV'%self.v_rmse_ave)

    def _plot_pred_parity(self, ax, INTERVAL = 1, cl='dodgerblue'):
        
        idx = 0
        parity_plot(ax[idx], self.ener_dft,self.ener_dp, INTERVAL = INTERVAL, cl=cl)
        ax[idx].set_xlabel('$E_{DFT}\ ({\\rm eV/atom})$',fontsize=10)  
        ax[idx].set_ylabel('$E_{DP}\ ({\\rm eV/atom})$',fontsize=10)  
        ax[idx].set_title('$\sigma_E$ = %.3f meV/atom'%(self.E_rmse_ave*1e3),fontsize=10)  

        idx = 1
        parity_plot(ax[idx], self.norm_dft,self.norm_dp, INTERVAL =  INTERVAL, cl=cl)
        ax[idx].set_xlabel('$\parallel f_{DFT}\parallel \ ({\\rm eV/\mathring A})$',fontsize=10)  
        ax[idx].set_ylabel('$\parallel f_{DP}\parallel \ ({\\rm eV/\mathring A})$',fontsize=10)  
        ax[idx].set_title('$\sigma_f$ = %.4f ${\\rm eV/\mathring A}$'%self.f_rmse_ave,fontsize=10)  

        idx = 2

        try:
            parity_plot(ax[idx], self.p_dft,self.p_dp, INTERVAL =  INTERVAL, cl=cl)
            ax[idx].set_xlabel('$p_{DFT}\ ({\\rm GPa})$',fontsize=10)  
            ax[idx].set_ylabel('$p_{DP}\ ({\\rm GPa})$',fontsize=10)  
            ax[idx].set_title('$\sigma_p$ = %.4f GPa'%np.average(self.p_rmse),fontsize=10)  
        except:
            parity_plot(ax[idx], self.alg_dft,self.alg_dp, INTERVAL =  INTERVAL, cl=cl)
            ax[idx].set_xlabel('$v_{DFT}\ ({\\rm eV})$',fontsize=10)  
            ax[idx].set_ylabel('$v_{DP}\ ({\\rm eV})$',fontsize=10)  
            ax[idx].set_title('$\sigma_v$ = %.4f eV'%np.average(self.v_rmse),fontsize=10)  


    def _plot_RMSE_PT(self, ax, INTERVAL):
        

        rmse = _get_RMSE(self.ener_dft,self.ener_dp,)

        sr = plt.scatter(self.press[::INTERVAL], self.temps[::INTERVAL], 
                         c=rmse[::INTERVAL]*1e3,
                         s=10, marker='|',  cmap='coolwarm')

        plt.colorbar(sr, label='RMSE of Energy (meV/atom)')
        ax.set_xlabel('Pressure (GPa)')
        ax.set_ylabel('Temperature (K)')
