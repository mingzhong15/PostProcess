from Constant import *

const_I = np.sqrt(np.pi/4 / np.log(2))

def _heat_capacity(X, a1, a2, a3, a4):
    
    return a1 * X + a2 * X**2 + a3 * X**3 + a4 * X**4 
#c0 (a0 + +  *np.exp(-(A*X)**2))

def _fit_it(data, init, end):
    
    p_est, err_est = curve_fit(_heat_capacity, data[init:end,0], data[init:end,1])
    
def _cut_profile(DIR):
    
    file = open(DIR+'/MD.profile','r')
    outfile = open(DIR+'/cut.profile','w')

    lines =file.readlines()

    for line in lines:

        if len(line.split()) == 9:
            outfile.writelines(line)

    file.close()
    outfile.close()

class TTMSys():
    
    def __init__(self, mass_mole, a, latt_type):
        
        self.mass = mass_mole / kg2gmol
        self.a = a
        
        if latt_type == 'bcc':
            cell_atom = 2
        elif latt_type == 'fcc':
            cell_atom = 4
                
        # unit : 1/m^3
        self.num_rho = cell_atom /(self.a/m2A)**3       
        
        # unit : kg/m^3
        self.rho = self.num_rho * self.mass         
        #self.rho_gcc = self.rho * kg2g/(m2cm**3)
        
        print('density %.2f g/cm^3'%(self.rho * kg2g/(m2cm**3)))
        
    def set_system_size(self, L):
        
        self.L = L
        
        r_verlet = 6.10 + 1.0     # unit: angstrom

        self.Nx = int(self.L * 10 /self.a ) + 1

        delta = r_verlet / 2

        self.N_grid = int(L * 10 / delta) +1

        print('Lattice :',self.Nx, delta, '\nGrid for FD :', self.N_grid, self.N_grid * 3)

        self.l_surface = int(self.N_grid)
        self.r_surface = int(self.N_grid*2-1) 

        print('left Surface: ',self.l_surface, self.r_surface)        
        
        self.coords = np.arange(self.N_grid*3) * delta /10 - L
        
    # ==================================== #
    #  Generate input script
    # ==================================== #     
    
    def _generate_init_TE_grid(self, savedir, T_0):
        
        grid_x = np.arange(self.N_grid * 3).astype(int)

        T_init = np.vstack( (grid_x, np.zeros(grid_x.shape), np.zeros(grid_x.shape), np.zeros(grid_x.shape) ))

        thick = self.N_grid

        TE_field = T_init.T.astype(int)

        TE_field[self.N_grid:self.N_grid*2,-1] += T_0
        
        print('\n from grid %.d to %.d, electronic temperature is initialized at %.1f kelvin \n'%(self.N_grid,self.N_grid*2-1,T_0))
        
        np.savetxt(savedir+'TTM_init.dat',TE_field, delimiter='   ', fmt='%.d')
        
    def _obtain_gamma(self, G):
        
        # unit : kg/s
        gamma_p_SI = self.mass  /(3 * self.num_rho * kb) * G      
        self.gamma_p_metal = gamma_p_SI * kg2gmol / s2ps
        
        print('e-p coupling: ',self.gamma_p_metal, 'g/mol / ps')


    def _obtain_laser(self, tau_L, epsilon=0, F_abs=0):
        # Gaussian width, unit: ps
        self.sigma = tau_L / np.sqrt(8 * np.log(2))  
        
        if epsilon !=0 :
            I_abs_SI = epsilon * self.rho * self.L * nm2A / m2A / np.sqrt(2*np.pi) / self.sigma
            self.I_abs = I_abs_SI * J2eV / (m2A**2)
        elif F_abs !=0 :
            self.I_abs = F_abs / const_I / tau_L * J2eV / (m2A**2)

        print('Intensity (eV/A^2 ps): ', self.I_abs)
        print('laser pulse (ps): ', self.sigma)
        
    def _plot_laser(self):
        t = np.linspace(0, self.sigma*10, 100)
        t_0 = 5 * self.sigma

        fig, ax = plt.subplots(figsize=(2,2),dpi=200)

        ax.plot(t * 1e3, self.I_abs * np.exp(-(t-t_0)**2/(2*self.sigma**2)) )#, label='Laser Pulse')

        ax.axhline(self.I_abs, linestyle=':',color='k')
        ax.axhline(self.I_abs/2, linestyle=':',color='k')

        ax.set_xlim(0,)
        ax.set_xlabel('Time (fs)')
        ax.set_ylabel('Intensity $\\rm{eV/(\mathring A^2 \cdot ps)}$)')
        
    def _fit_capacity(self, data, end, init):
        
        #compare_root = '/data/home/djy4/compare/tungsten/'
        #data = np.loadtxt(file)
        #data = np.loadtxt(compare_root+file)
        #data[:,0] /= 1e3 
        
        fig, ax = plt.subplots(figsize=(2,2),dpi=200)
        ax.plot(data[:,0],data[:,1],'-',color='r',mfc='none',mew=0.5, linewidth=1.0, label='Raw Data')
        ax.plot(data[:,0],data[:,0]* 1e3* 137.3,':',color='k',mfc='none',mew=0.5, linewidth=1.0, label='FEG')
        
        # ======== lower fit ========== #
        p_est, err_est = curve_fit(_heat_capacity, data[:end,0], data[:end,1])
        self.temp_thres = data[init,0]
        self.esheat = p_est * J2eV /(m2A**3)
        
        X = np.linspace(0,data[init,0],100)
        Ce_fit = _heat_capacity(X, *p_est)
        ax.plot(X, Ce_fit,':',color='b',mfc='none',mew=0.5, linewidth=2.5, label='PolyFit (Low)')
        print('X < %.2f '%data[init,0],'low temperature (a1-a4): ',p_est * J2eV /(m2A**3))

        # ======== higher fit ========== #
        p_est, err_est = curve_fit(_heat_capacity, data[init:,0], data[init:,1])
        self.h_esheat = p_est * J2eV /(m2A**3)
        X = np.linspace(data[init,0], data[-1,0],100)
        Ce_fit = _heat_capacity(X, *p_est)
        ax.plot(X, Ce_fit,':',color='cyan',mfc='none',mew=0.5, linewidth=2.5, label='PolyFit (High)')

        print('high temperature (a1-a4): ',p_est * J2eV /(m2A**3) )
        ax.legend(fontsize=6)
        #ax.set_xlim(0, data[-1,0])
        ax.set_xlim(0,20)
        ax.set_ylim(0,np.max(data[:,1]))
        ax.set_ylabel('$C_e$ ($J/(m^3\cdot K)$)',fontsize=8)
        ax.set_xlabel('Temperature ($10^3K$)',fontsize=8)
        
    def _generate_param(self, savedir):

        file = open(savedir+'parameter.txt','w')
        
        # == part 1. low temperautre C_e polyfit == #
        file.writelines(' =========== [*esheat_0] =========== # \n')
        temp = 0
        file.writelines('%.1f \n'%temp)
        file.writelines(' =========== [*esheat_1] =========== # \n')
        file.writelines('%.32f \n'%self.esheat[0])
        file.writelines(' =========== [*esheat_2] =========== # \n')
        file.writelines('%.32f \n'%self.esheat[1])
        file.writelines(' =========== [*esheat_3] =========== # \n')
        file.writelines('%.32f \n'%self.esheat[2])     
        file.writelines(' =========== [*esheat_4] =========== # \n')
        file.writelines('%.32f \n'%self.esheat[3]) 
        file.writelines(' =========== [C_limit] =========== # \n')
        temp = 0
        file.writelines('%.1f \n'%temp)
        file.writelines(' =========== [T_damp] =========== # \n')
        temp = 0
        file.writelines('%.1f \n'%temp)        
        # ======================================= #
        
        file.writelines(' =========== [rho_e] =========== # \n')
        temp = 1
        file.writelines('%.1f \n'%temp)      
        
        file.writelines(' =========== [el_th_diff] =========== # \n')
        temp = 0.0
        file.writelines('%.1f \n'%temp)
        
        file.writelines(' =========== [*gamma_p] =========== # \n')
        file.writelines('%.32f \n'%self.gamma_p_metal)             
        file.writelines(' =========== [gamma_s] =========== # \n')
        temp = 0.0
        file.writelines('%.1f \n'%temp)          
        file.writelines(' =========== [v_0] =========== # \n')
        temp = 0.0
        file.writelines('%.1f \n'%temp)        
        file.writelines(' =========== [*intensity] =========== # \n')
        file.writelines('%.32f \n'%self.I_abs)
        file.writelines(' =========== [*surface_l] =========== # \n')
        file.writelines('%.d \n'%self.l_surface)    
        file.writelines(' =========== [*surface_r] =========== # \n')
        file.writelines('%.d \n'%self.r_surface)       
        
        file.writelines(' =========== [*skin_layer] =========== # \n')
        file.writelines('%.d \n'%(self.r_surface-self.l_surface+1))          
        file.writelines(' =========== [*width] =========== # \n')
        file.writelines('%.32f \n'%self.sigma)          
        file.writelines(' =========== [pres_factor] =========== # \n')
        temp = 0.0
        file.writelines('%.1f \n'%temp)                
        file.writelines(' =========== [free_path] =========== # \n')
        temp = 0.0
        file.writelines('%.1f \n'%temp)                  
        file.writelines(' =========== [ionic_density] =========== # \n')
        file.writelines('%.32f \n'%(self.num_rho /(m2A)**3)  )     
        file.writelines(' =========== [*movsur]=========== # \n')
        temp = 1
        file.writelines('%.d \n'%temp)          
        file.writelines(' =========== [electron_temperature_min]=========== # \n')
        temp = 300
        file.writelines('%.2f \n'%temp)          
        
        # == extra setttings == #
        file.writelines(' =========== *extra_steps=========== # \n')
        temp = 1
        file.writelines('%.d \n'%temp)  
        file.writelines(' =========== *minimum to activate new grid=========== # \n')
        temp = 40
        file.writelines('%.d \n'%temp)      
        
        # == part 1. high temperautre C_e polyfit == #        
        file.writelines(' =========== [temp_threshold] =========== # \n')
        file.writelines('%.32f \n'%self.temp_thres  )             
        file.writelines(' =========== [*h_esheat_1] =========== # \n')
        file.writelines('%.32f \n'%self.h_esheat[0])
        file.writelines(' =========== [*h_esheat_2] =========== # \n')
        file.writelines('%.32f \n'%self.h_esheat[1])
        file.writelines(' =========== [*h_esheat_3] =========== # \n')
        file.writelines('%.32f \n'%self.h_esheat[2])     
        file.writelines(' =========== [*h_esheat_4] =========== # \n')
        file.writelines('%.32f \n'%self.h_esheat[3]) 
                    
        file.close()
        
        
    # ==================================== #
    #  Read Results from Simulations
    # ==================================== #
        
    def _read_TTM_out(self, DIR, delta_t):
        data = np.loadtxt(DIR+'/TTM_out.dat')  #  max_rows=4000
        
        self.TI = data[:,1:self.N_grid*3+1]
        self.TE = data[:,self.N_grid*3+1:]
        self.times = np.arange(self.TE.shape[0])*delta_t - 5*self.sigma
        
        print(np.max(self.TE[:,self.l_surface]))
        
    def _obtain_laser(self, tau_L, epsilon=0, F_abs=0):
        # Gaussian width, unit: ps
        self.sigma = tau_L / np.sqrt(8 * np.log(2))  
        
        if epsilon !=0 :
            I_abs_SI = epsilon * self.rho * self.L * nm2A / m2A / np.sqrt(2*np.pi) / self.sigma
            self.I_abs = I_abs_SI * J2eV / (m2A**2)
        elif F_abs !=0 :
            self.I_abs = F_abs / const_I / tau_L * J2eV / (m2A**2)

        print('Intensity (eV/A^2 ps): ', self.I_abs)
        print('laser pulse (ps): ', self.sigma)
        
    def _get_TE_evolution_centroid(self, DELTA):

        self.TE_AVE = np.average(self.TE[:, self.l_surface + int(self.N_grid/2) - DELTA : self.l_surface + int(self.N_grid/2) + DELTA], axis=1 )
        self.TI_AVE = np.average(self.TI[:, self.l_surface + int(self.N_grid/2) - DELTA : self.l_surface + int(self.N_grid/2) + DELTA], axis=1 )
        
        self.TI_STD = np.std(self.TI[:, self.l_surface + int(self.N_grid/2) - DELTA : self.l_surface + int(self.N_grid/2) + DELTA], axis=1 ) 
    

    
    def _plot_TE_evolution(self, ax, DELTA, show_error=False):

        self._get_TE_evolution_centroid(DELTA)

        ax.plot(self.times, self.TE_AVE,'-', mfc='none',mew=0.5, 
                linewidth=1.0, label='$T_e$ (bulk)')
        
        if show_error:
            ax.errorbar(self.times, self.TI_AVE, yerr=self.TI_STD, fmt='-o',mew=0.5, ms=2,
                linewidth=1.0, label='$T_i$ (bulk)')            
        else:
            ax.plot(self.times, self.TI_AVE,'-', mfc='none',mew=0.5, 
                linewidth=1.0, label='$T_i$ (bulk)')

    def _plot_spatial(self,):
        
        X,Y = np.meshgrid(self.times, self.coords)
        fig,ax = plt.subplots(figsize=(3,2),dpi=200) 
        delta_x = 15
        
        cs = ax.contourf(X, Y, self.TE.T, 20,  cmap=plt.cm.Spectral_r)
        plt.colorbar(cs,label='$T_e$ (K)')
        ax.set_ylabel('Depth (nm)')
        ax.set_xlabel('Time (ps)')
        ax.set_ylim(0-delta_x, 30+delta_x)
        
        fig,ax = plt.subplots(figsize=(3,2),dpi=200) 
        cs = ax.contourf(X, Y, self.TI.T, 20,  cmap=plt.cm.Spectral_r)
        plt.colorbar(cs,label='$T_i$ (K)')
        ax.set_ylabel('Depth (nm)')
        ax.set_xlabel('Time (ps)')
 
        ax.set_ylim(0-delta_x, 30+delta_x)
    
    def _read_profile(self, DIR, delta_t):

        data_profile = np.loadtxt(DIR+ '/cut.profile')

        N_chunk = int(np.max(data_profile[:,0]))

        self.coords = data_profile[:,1].reshape(-1,N_chunk)[0]/10    
        self.ncount = np.array( data_profile[:,2].reshape(-1,N_chunk) )
        self.rho = np.array( data_profile[:,3].reshape(-1,N_chunk) )
        self.temps = np.array( data_profile[:,4].reshape(-1,N_chunk)  )
        self.press = - np.array( data_profile[:,5].reshape(-1,N_chunk) /1e4 )
        self.Q4 = np.array( data_profile[:,6].reshape(-1,N_chunk) )
        self.Q6 = np.array( data_profile[:,7].reshape(-1,N_chunk) )
        self.Q8 = np.array( data_profile[:,8].reshape(-1,N_chunk) )

        self.times = np.arange(self.temps.shape[0])*delta_t

        self.X,self.Y = np.meshgrid(self.times, self.coords)
        
    def _get_centroid(self, data, DELTA):
        center = int((self.l_surface+self.r_surface)/2)
        
        return np.average(data[:, center - DELTA : center + DELTA], axis=1 )
    
    def _get_surface(self, data, DELTA):
        
        data_front = []
        data_rear = []

        for i in range(self.TE.shape[0]):
            
            idx_f = int(self.front[i])
            idx_r = int(self.rear[i])
            
            data_front = np.append(data_front, np.average(data[i, idx_f :idx_f+DELTA]) )
            data_rear = np.append(data_rear, np.average(data[i, idx_r-DELTA:idx_r ]) )
        
        return data_front, data_rear
        
    def _get_boundary(self):
        
        self.front = []
        self.rear = []
        
        for i in range(self.TE.shape[0]):
            
            idx_f = np.argwhere(self.TE[i]!=0)[0][0]
            idx_r = np.argwhere(self.TE[i]!=0)[-1][0]

            self.front = np.append(self.front, idx_f)
            self.rear = np.append(self.rear, idx_r)

        self.front.astype(int)
        self.rear.astype(int)

        
    def _plot_PT(self, ax,  INTERVAL, color, label):
        
        DELTA = 2

        self.temps_c = self._get_centroid(self.temps, DELTA)
        self.press_c = self._get_centroid(self.press, DELTA)

        ax.plot( self.press_c[::INTERVAL], self.temps_c[::INTERVAL],
            'o',color=color,mew=0.5,alpha=0.6, ms=4,label=label)   