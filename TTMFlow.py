from Constant import *

sigma2duration = 2 * np.sqrt(2*np.log(2))

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
                
        self.num_rho = cell_atom /(self.a/m2A)**3   # [1/m^3]
        
        self.rho = self.num_rho * self.mass         # [kg/m^3]
        #self.rho_gcc = self.rho * kg2g/(m2cm**3)
        
        print('density %.2f g/cm^3'%(self.rho * kg2g/(m2cm**3)))
        
    def set_system_size(self, L, mode='free', delta =3.55, min_act = 80):
        
        self.L = L
        #r_verlet = 6.10 + 1.0     # unit: angstrom
        #delta = r_verlet / 2
        
        self.Nx = int(self.L * 10 /self.a ) + 1
        self.N_grid = int(L * 10 / delta) +1
        self.delta = delta

        self.boundary_mode = mode
        if self.boundary_mode == 'free':
            self.N_expand = 3
        elif self.boundary_mode == 'substrate':
            self.N_expand = 2
            
        print('Lattice : ',self.Nx, self.delta, '\nGrid for FD :', 
              self.N_grid, self.N_grid * self.N_expand)

        self.l_surface = int(self.N_grid)
        self.r_surface = int(self.N_grid*2-1) 

        self.coords = np.arange(self.N_grid*self.N_expand) * self.delta /10 - L

                    
        print('left Surface: %.d \n right Surface: %.d \n'%(self.l_surface, self.r_surface))        
        
        self.min_act = min_act
    # ==================================== #
    #  Generate input script
    # ==================================== #     
    
    def _generate_init_TE_grid(self, savedir, T_0):
        
        grid_x = np.arange(self.N_grid * self.N_expand).astype(int)

        T_init = np.vstack( (grid_x, np.zeros(grid_x.shape), np.zeros(grid_x.shape), np.zeros(grid_x.shape) ))

        thick = self.N_grid

        TE_field = T_init.T.astype(int)

        TE_field[self.N_grid:self.N_grid*2,-1] += T_0
        
        print('\n from grid %.d to %.d, electronic temperature is initialized at %.1f kelvin \n'%(self.N_grid,self.N_grid*2-1,T_0))
        
        np.savetxt(savedir+'TTM_init.dat',TE_field, delimiter='   ', fmt='%.d')

    
    def _obtain_gamma(self, G0, data, init, const=False):
        
        print('=====================================')
        print('e-ph coupling setting starts')
        
        
        # unit : kg/s
        gamma_p_SI = self.mass  /(3 * self.num_rho * kb) * G0      
        self.gamma_p_metal = gamma_p_SI * kg2gmol / s2ps
        
        print('e-p coupling (G0): ',self.gamma_p_metal, 'g/mol / ps')
        
        
        fig, ax = plt.subplots(figsize=(2,2),dpi=200)
        ax.plot(data[:,0],data[:,1],'-',color='silver',mfc='none',mew=0.5, linewidth=1.0, label='Raw Data')
        
        gei_metal = self.mass  /(3 * self.num_rho * kb) * data[:,1] * kg2gmol / s2ps    
        # ======== higher fit ========== #
        
        
        if const:
            self.thres_G = 100
        else:
            self.thres_G = data[init,0]
        
        p_est, err_est = curve_fit(_heat_capacity, data[init:,0], gei_metal[init:]) #data[init:,1])
        self.gp = p_est #self.mass  /(3 * self.num_rho * kb) * p_est * kg2gmol / s2ps
        
        X = np.linspace(data[init,0], data[-1,0],100)
        gamma_fit = _heat_capacity(X, *p_est)
        gei_fit =  gamma_fit / (kg2gmol / s2ps) * (3 * self.num_rho * kb)  / self.mass  
        
        ax.plot(X, gei_fit,':',color='blue',mfc='none',mew=0.5, linewidth=1.5, label='PolyFit (High)')

        if const:
            print('constant G0 is set')
        else:
            print('above %.4f 10^3 K '%self.thres_G,'high temperature (a1-a4): ', self.gp )
        ax.legend(fontsize=6)
        ax.set_xlim(0, data[-1,0])
        #ax.set_xlim(0,30)
        ax.set_ylim(0,np.max(data[:,1]))
        ax.set_ylabel('$g_{ei}$ ($W/(m^3\cdot K)$)',fontsize=8)
        ax.set_xlabel('Temperature ($10^3K$)',fontsize=8)
  
        print('e-ph coupling setting ends')
        print('=====================================')
        
        
        return ax
               
    def _obtain_laser(self, tau_L, skin_layer=10, is_decay=False, epsilon=0, F_abs=0,  is_print=True):
    
    # ======================
    # tau_L [ps], skin_layer [nm], epsilon [J/kg], F [J/m^2]
    # L     [nm], rho     [kg/m^3],
    # ======================

        self.sigma = tau_L / sigma2duration                   # [ps]
        self.skin_layer = int(skin_layer*nm2A/self.delta)     # [int]
        
        self.is_decay = is_decay
        decay_factor = self.L / skin_layer / (1 - np.exp( - self.L / skin_layer))
        
        
        fluence2intensity = 1/np.sqrt(2*np.pi)/self.sigma     # [ps^-1]
        epsilon2fluence   = self.rho * (self.L* nm2A/m2A)           # [kg/m^2]

        if epsilon != 0 :
            F_abs = epsilon * epsilon2fluence                # [J/m^2] <-- [J/kg] * [kg/m^2]
     
        if self.is_decay:
            self.I_abs = F_abs * fluence2intensity * J2eV/m2A**2 * decay_factor
        else:
            self.I_abs = F_abs * fluence2intensity * J2eV/m2A**2 # [eV/A^2/ps] 
        
        if is_print:
            if self.is_decay:
                print('exponential decay with optical peneration depth : %.d grids'%(self.skin_layer))
            else:
                print('uniform heating condition')
            print('Intensity (eV/A^2 ps): ', self.I_abs)
            print('laser pulse (ps): ', self.sigma)

    def _set_ele_diff(self, v_fermi, A_temp, B_temp):
        
        self.v_fermi = v_fermi
        self.A_temp = A_temp
        self.B_temp = B_temp
            
    def _plot_laser(self):
        t = np.linspace(0, self.sigma*10, 100)
        t_0 = 5 * self.sigma

        fig, ax = plt.subplots(figsize=(2,2),dpi=200)

        pulse = self.I_abs * np.exp(-(t-t_0)**2/(2*self.sigma**2))
        ax.plot(t * 1e3,  pulse )#, label='Laser Pulse')

        ax.axhline(self.I_abs, linestyle=':',color='k')
        ax.axhline(self.I_abs/2, linestyle=':',color='k')
        
        #print('FWHM : %.7f ps'%)

        ax.set_xlim(0,)
        ax.set_xlabel('Time (fs)')
        ax.set_ylabel('Intensity $\\rm{eV/(\mathring A^2 \cdot ps)}$)')
        
    def _fit_capacity(self, data, end, init):

        self.Ce_data = data
        self.Ce_end = end
        self.Ce_init = init
        
        # ======== lower fit ========== #
        p_est, err_est = curve_fit(_heat_capacity, data[:end,0], data[:end,1])
        self.temp_thres = data[init,0]
        self.esheat = p_est * J2eV /(m2A**3)
        
        print('X < %.2f '%data[init,0],'low temperature (a1-a4): ',p_est * J2eV /(m2A**3))

        # ======== higher fit ========== #
        p_est, err_est = curve_fit(_heat_capacity, data[init:,0], data[init:,1])
        self.h_esheat = p_est * J2eV /(m2A**3)
        print('high temperature (a1-a4): ',p_est * J2eV /(m2A**3) )
    
    def _plot_capacity(self):
                
        fig, ax = plt.subplots(figsize=(2,2),dpi=200)
        ax.plot(self.Ce_data[:,0],self.Ce_data[:,1],'-',color='silver',
                mew=1.5, linewidth=1.0, label='Raw Data')
        
        X = np.linspace(0,self.Ce_data[self.Ce_init,0],100)
        Ce_fit = _heat_capacity(X, *self.esheat/J2eV*(m2A**3))
        ax.plot(X, Ce_fit,'-',color='navy',mfc='none',mew=0.5, linewidth=0.5, label='PolyFit (Low)')

        X = np.linspace(self.Ce_data[self.Ce_init,0], self.Ce_data[-1,0],100)
        Ce_fit = _heat_capacity(X, *self.h_esheat/J2eV*(m2A**3))
        
        ax.plot(X, Ce_fit,'-',color='crimson',mfc='none',mew=0.5, linewidth=0.5, label='PolyFit (High)')

        ax.legend(fontsize=6)
        ax.set_xlim(0, self.Ce_data[-1,0])
        #ax.set_xlim(0,30)
        ax.set_ylim(0,np.max(self.Ce_data[:,1]))
        ax.set_ylabel('$C_e$ ($J/(m^3\cdot K)$)',fontsize=8)
        ax.set_xlabel('Temperature ($10^3K$)',fontsize=8)
        
        return ax    
        
    def _generate_param(self, savedir, savefile):

        file = open( os.path.join(savedir,savefile),'w')
        
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
        
        file.writelines(' =========== [v_fermi] (ang/ps)=========== # \n')
        file.writelines('%.2f \n'%self.v_fermi)   
        file.writelines(' =========== [A_temp] (K^-2 ps^-1)=========== # \n')
        file.writelines('%.32f \n'%self.A_temp)  
        file.writelines(' =========== [B_temp] (K^-1 ps^-1) =========== # \n')
        file.writelines('%.32f \n'%self.B_temp)          
        
        file.writelines(' =========== [*gamma_0] =========== # \n')
        file.writelines('%.32f \n'%self.gamma_p_metal)             
        file.writelines(' =========== [*thres_G] =========== # \n')
        file.writelines('%.32f \n'%self.thres_G)                 
        file.writelines(' =========== [*gp_1] =========== # \n')
        file.writelines('%.32f \n'%self.gp[0])
        file.writelines(' =========== [*gp_2] =========== # \n')
        file.writelines('%.32f \n'%self.gp[1])
        file.writelines(' =========== [*gp_3] =========== # \n')
        file.writelines('%.32f \n'%self.gp[2])     
        file.writelines(' =========== [*gp_4] =========== # \n')
        file.writelines('%.32f \n'%self.gp[3])         
        
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
        #file.writelines('%.d \n'%(self.r_surface-self.l_surface+1))
        file.writelines('%.d \n'%(self.skin_layer))
        file.writelines(' =========== [*is_decay] =========== # \n')
        if self.is_decay:
            temp = 1
        else:
            temp = 0
        file.writelines('%.d \n'%temp)
        
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
        temp = self.min_act
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
        
        self.TI = data[:,1:self.N_grid*self.N_expand+1]
        self.TE = data[:,self.N_grid*self.N_expand+1:]
        self.times = np.arange(self.TE.shape[0])*delta_t - 5*self.sigma
        
        mask = np.ma.masked_equal(self.TI, 0)
        self.TI_min = np.ma.min(mask, axis=1).data
        self.TI_max = np.ma.max(mask, axis=1).data
        self.TI_ave = np.ma.mean(mask, axis=1).data
        
        mask = np.ma.masked_equal(self.TE, 0)
        self.TE_min = np.ma.min(mask, axis=1).data
        self.TE_max = np.ma.max(mask, axis=1).data
        self.TE_ave = np.ma.mean(mask, axis=1).data
        
        print(np.max(self.TE_max))
        

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
        fig,ax = plt.subplots(figsize=(2,2),dpi=200) 
        delta_x = 15
        
        cs = ax.contourf(X, Y, self.TE.T, 20,  cmap=plt.cm.Spectral_r)
        plt.colorbar(cs,label='$T_e$ (K)')
        ax.set_ylabel('Depth (nm)')
        ax.set_xlabel('Time (ps)')
        ax.set_ylim(0-delta_x, 30+delta_x)
        
        fig,ax = plt.subplots(figsize=(2,2),dpi=200) 
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
        #self.Q8 = np.array( data_profile[:,8].reshape(-1,N_chunk) )

        self.times = np.arange(self.temps.shape[0])*delta_t

        self.X,self.Y = np.meshgrid(self.times, self.coords)
        
    def _get_average(self, data):
        
        cri = self.ncount >= self.min_act
        data[cri] = 0
        
        mask = np.ma.masked_equal(data, 0)
        #cri = sys.ncount[:,1:] >= sys.min_act
        return np.ma.mean(mask, axis=1).data
        
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