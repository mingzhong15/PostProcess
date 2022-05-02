from Constant import *
import dpdata

class msstTraj():
    
    def __init__(self, traj_dir, printf = True):
        
        self.traj_dir = traj_dir
        self.thermo_list = glob.glob( os.path.join(traj_dir, 'output_*'))
        self.thermo_list.sort()
        
        if printf:
            print('%.d thermo file found :'%(len(self.thermo_list)))
            self._print_list(self.thermo_list)
                
        self.dump_list = glob.glob( os.path.join(traj_dir, 'dump.*'))
        self.dump_list.sort()
        
        if printf:
            print('%.d dump file found :'%(len(self.dump_list)))
            self._print_list(self.dump_list)

    def _print_list(self, file_list):
        for file in file_list:
            print(file.split('/')[-1])
            
    def _creat_label(self, prefix):
        self.label_list = []
        for output in self.thermo_list:
            #label = prefix+ output.split('/')[-1].split('.')[0].split('_')[-1]
            
            label = prefix + output.split('_')[-1][:-4]
            self.label_list.append(label)

    def _evolution_plot(self, INIT, END):
        
        fig, ax = plt.subplots(2,1,figsize=(3,4),dpi=200)

        for i in range(len(self.label_list)):

            times, temps, press = _get_shock(self.thermo_list[i])

            i = 0
            ax[i].plot(times[INIT:END], temps[INIT:END], label=self.label_list[i])
            ax[i].set_ylabel('Temperature (k)')

            i = 1
            ax[i].plot(times[INIT:END], press[INIT:END])
            ax[i].set_ylabel('Preussure (GPa)')

            ax[i].set_xlabel('Time (ps)')

            print(np.average(temps[INIT:END]), np.average(press[INIT:END]))

        ax[0].legend(fontsize=8)
        
    def _sampling_plot(self, ax, INIT, END, INTERVAL, prefix):
        
        for i in range(len(self.label_list)):
            times, temps, press = _get_shock(self.thermo_list[i])
            print(times.shape)
            
            ax.plot(press[INIT:END:INTERVAL], temps[INIT:END:INTERVAL], 'o', 
            ms=3, alpha=0.5, mew=0.5,label=prefix+self.label_list[i])
            print(press[INIT:END:INTERVAL].shape)
            
    def _generate_vasp(self):
        
        for i in range(len(self.label_list)):
            
            out_dir = os.path.join(self.traj_dir, self.label_list[i])
            print(out_dir, ' is working')
            
            os.system('mkdir '+out_dir)
            
            sys = dpdata.System(self.dump_list[i],fmt='lammps/dump',type_map=['C','H'])
            data = np.loadtxt(self.thermo_list[i], skiprows=1)
            
            temps = data[:,1]
            sigma = data[:,1]*kb*J2eV
            Nf = sys.get_nframes()
            Nt = temps.shape[0]
            
            print('%.d Configurations , %.d Thermodynamic results'%(Nf,Nt))
    
            for kk in range(Nt):
                sys.to('vasp/poscar', os.path.join(out_dir, 'POSCAR.%.d'%kk), frame_idx=kk)
                
                incar_kk = os.path.join(out_dir,'INCAR.%.d'%kk)
                os.system('cp '+os.path.join(self.traj_dir,'INCAR')+' '+ incar_kk)
                os.system("sed -i 's/SIGMA = xx/SIGMA = %.10f/g' "%sigma[kk]+' '+incar_kk)
            
            print(self.label_list[i],' is done, %.d frames are added'%Nt)     
            
def _get_shock(DIR):
    
    data = np.loadtxt(DIR,skiprows=1)
    
    temps = data[:,1]
    press = data[:,4]*bar2Pa/1e9
    
    return data[:,0],temps, press