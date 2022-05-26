from Constant import *

def dump_to_xyz(FILE, OUTFILE):

    file = open(FILE, 'r')
    outfile = open(OUTFILE, 'w')

    line = file.readline()

    while line:

        if 'NUMBER OF ATOMS' in line:
            line = file.readline()
            natoms = int(line)
            outfile.write(line)

        if 'BOX BOUNDS' in line:
            line = file.readline()

            #print(line.split())
            xl, xh = line.split()
            ret ='%.11f 0 0 \t'%(float(xh) - float(xl))
            ret += '0 %.11f 0 \t'%(float(xh) - float(xl))
            ret += '0 0 %.11f\n'%(float(xh) - float(xl))

            outfile.write(ret)

        if 'ATOMS' in line:

            for i in range(natoms):
                line = file.readline()

                output = line.split()[1:]

                idx = int(output[0])
                x = float(output[1])
                y = float(output[2])
                z = float(output[3])

                outfile.write('%.d %.11f %.11f %.11f \n'%(idx,x,y,z))

        line = file.readline()
    
    file.close()
    outfile.close()
    
class lammpsTraj():
    
    def __init__(self, traj_dir, output_prefix, dump_prefix, printf = True, mode='msst'):
        
        self.traj_dir = traj_dir
        self.thermo_list = glob.glob( os.path.join(traj_dir, output_prefix))
        self.thermo_list.sort()
        
        self.mode = mode
        
        if printf:
            print('%.d thermo file found :'%(len(self.thermo_list)))
            self._print_list(self.thermo_list)
            
        self.dump_list = glob.glob( os.path.join(traj_dir, dump_prefix))
        self.dump_list.sort()
        
        if printf:
            print('%.d dump file found :'%(len(self.dump_list)))
            self._print_list(self.dump_list)                

    def _print_list(self, file_list):
        
        skip = len(self.traj_dir)
        for file in file_list:
            #print(file.split('/')[-1])    
            print(file[skip:])
    
    def _creat_label(self, prefix, mode='thermo'):
        self.label_list = []

        if mode == 'thermo':
            for output in self.thermo_list:
                label = prefix+ output.split('/')[-1].split('.')[0].split('_')[-1]
                #label = prefix + output.split('_')[-1][:-4]
                self.label_list.append(label)
        elif mode == 'dump':
            for output in self.dump_list:
                label = prefix + output.split('/')[-1].split('.')[1]
                #label = prefix + output.split('_')[-1][:-4]
                self.label_list.append(label)
        elif mode == 'dirs':
            for output in self.dump_list:
                label = prefix + output.split('/')[-2].split('.')[-1]
                #label = prefix + output.split('_')[-1][:-4]
                self.label_list.append(label)

    def _evolution_plot(self, INIT, END, delta_t=1):
        
        fig, ax = plt.subplots(2,1,figsize=(3,4),dpi=200)

        for idx in range(len(self.label_list)):

            times, temps, press = _get_thermo(self.thermo_list[idx], self.mode, delta_t)

            i = 0
            ax[i].plot(times[INIT:END], temps[INIT:END], label=self.label_list[idx])
            ax[i].set_ylabel('Temperature (k)')

            i = 1
            ax[i].plot(times[INIT:END], press[INIT:END])
            ax[i].set_ylabel('Pressure (GPa)')

            ax[i].set_xlabel('Time (ps)')

            print(np.average(temps[INIT:END]), np.average(press[INIT:END]))

        ax[0].legend(fontsize=8)
        
    def _sampling_plot(self, ax, INIT, END, INTERVAL, prefix, delta_t=1):
        
        for idx in range(len(self.label_list)):
            times, temps, press = _get_thermo(self.thermo_list[idx], self.mode, delta_t)
            
            print(times.shape)
            
            ax.plot(press[INIT:END:INTERVAL], temps[INIT:END:INTERVAL], 'o', 
            ms=3, alpha=0.5, mew=0.5,label=prefix+self.label_list[idx])
            
            print(press[INIT:END:INTERVAL].shape)     
            
    def _generate_vasp(self, type_map):
        
        for i in range(len(self.label_list)):
            
            work_dir = os.path.join(self.traj_dir, self.label_list[i])
            print(work_dir, ' is working')
            
            os.system('mkdir '+work_dir)
            
            sys = dpdata.System(self.dump_list[i],fmt='lammps/dump',type_map=type_map)
            
            times, temps, press = _get_thermo(self.thermo_list[i], self.mode, delta_t=1)
            sigma = temps*kb*J2eV
            Nf = sys.get_nframes()
            Nt = temps.shape[0]
            
            print('%.d Configurations , %.d Thermodynamic results'%(Nf,Nt))
    
            for kk in range(Nt):
                sys.to('vasp/poscar', os.path.join(work_dir, 'POSCAR.%.d'%kk), frame_idx=kk)
                
                incar_kk = os.path.join(work_dir,'INCAR.%.d'%kk)
                os.system('cp '+os.path.join(self.traj_dir,'INCAR')+' '+ incar_kk)
                os.system("sed -i 's/SIGMA = xx/SIGMA = %.10f/g' "%sigma[kk]+' '+incar_kk)
            
            print(self.label_list[i],' is done, %.d frames are added'%Nt)     
            
    def _collect_data_to(self, out_dir):
        
        temps = []
        ms = dpdata.MultiSystems()
        
        for i in range(len(self.label_list)):
            
            work_dir = os.path.join(self.traj_dir, self.label_list[i])
            print(work_dir, ' is working')

            times, tt, press = _get_thermo(self.thermo_list[i], self.mode, delta_t=1)
            temps = np.append(temps, tt)
            Nt = tt.shape[0]
            
            for kk in range(Nt):
                sys = dpdata.LabeledSystem( os.path.join(work_dir, '%.d'%kk, 'OUTCAR') )
                ms.append(sys)
                
            print(self.label_list[i],' is done, %.d frames are collected'%Nt)  
        
        os.system('mkdir '+out_dir)
        ms.to_deepmd_raw(out_dir)
        
        np.savetxt( os.path.join(out_dir,'fparam.raw'), temps)

        os.chdir(out_dir)
        os.system('mv ./C*/*.raw ./')
        os.system('~/raw_to_set.sh 5000')
            
            
def _get_thermo(DIR, mode='msst', delta_t=1):
    
    data = np.loadtxt(DIR,skiprows=1)
    
    if mode == 'msst':
        temps = data[:,1]
        press = data[:,4]*bar2Pa/1e9
    
        return data[:,0],temps, press

    elif mode == 'npt':

        temps = data[:,0]
        press = data[:,1]*bar2Pa/1e9
        times = np.arange(temps.shape[0]) * delta_t
        
        return times, temps, press      
