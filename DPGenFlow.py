from Constant import *
from Dataset import _get_EOS

from monty.serialization import loadfn,dumpfn

class DPGenSys():
    
    # ========================================= #
    # read basic information for DPGEN run
    # dpgen_dir : the working directory of DPGEN
    # ========================================= #
    def __init__(self, dpgen_dir, printf = True):
        
        self.dir = dpgen_dir
        self.param_file = os.path.join(dpgen_dir, 'param.json')
        self.jdata =loadfn(self.param_file)
        
        self.N_iter = len(glob.glob( os.path.join( dpgen_dir ,'iter.*' )))
        
        self.confs_list = self.jdata['sys_configs']
        self.label_list = []
        
        self.printf = printf
        
        print("DPGenerator System contains %.d Iterations"%(self.N_iter) )
        print("There are %.d initial configurations for exploration: \t "%(len(self.confs_list)) )
        for cc in self.confs_list:
            label = cc[0].split('.')[0]
            print(label,'\t\t')
            self.label_list.append(label)

    
    def _learning_curve(self, iter_idx, model_idx=0):
        
        iter_idx = 'iter.%.6d'%iter_idx
        model_idx = '%.3d'%model_idx
    
        path = os.path.join(self.dir,iter_idx, '00.train', model_idx, 'lcurve.out')

        return np.loadtxt(path)
        
    def _model_devi(self, iter_idx, sys_idx = 0, case_idx = 0):
        
        iter_idx = 'iter.%.6d'%iter_idx
        task_idx = 'task.%.3d'%sys_idx+'.%.6d'%case_idx  
        
        path = os.path.join(self.dir,iter_idx, '01.model_devi', task_idx, 'model_devi.out')

        return np.loadtxt(path)

    def _get_PT(self, iter_idx, sys_idx = 0, case_idx = 0):
        
        iter_idx = 'iter.%.6d'%iter_idx
        task_idx = 'task.%.3d'%sys_idx+'.%.6d'%case_idx  
        
        path = os.path.join(self.dir,iter_idx, '01.model_devi', task_idx, 'input.lammps')

        temp = float(os.popen("grep ' TEMP ' " + path).readlines()[0].split()[-1])
        press = float(os.popen("grep ' PRES ' " + path).readlines()[0].split()[-1]) * bar2Pa /1e9
        return temp, press
    
    def _sampling(self, iter_idx, sys_idx = 0):
        
        iter_idx = 'iter.%.6d'%iter_idx
        data_idx = 'data.%.3d'%sys_idx
        
        path = os.path.join(self.dir, iter_idx, '02.fp', data_idx)
        
        # return (temps, press, vol, energy)
        return _get_EOS(path)
    
    def _plot_lc(self, ax, iter_idx, model_idx=0):
        #fig, ax = plt.subplots(1,3, figsize=(8,2),dpi=200)

        data = self._learning_curve(iter_idx, model_idx)

        ll_list = ['Energy (eV)','Force (eV/$\\rm{\mathring A}$)','Virial (eV)']
        for idx in range(3):
            if idx == 0:
                ax[idx].plot(data[:,0],data[:,2+idx],'o',ms=1,mew=0.5, label='Iter.%.3d (DP%.3d)'%(iter_idx, model_idx))
            else:
                ax[idx].plot(data[:,0],data[:,2+idx],'o',ms=1,mew=0.5)

            ax[idx].set_xscale('log')
            ax[idx].set_yscale('log')
            ax[idx].set_title(ll_list[idx])
            ax[idx].set_xlabel('stopbatch')
            
    def _plot_md(self, ax, iter_idx, sys_idx = 0, case_idx = 0):
        
        lo= self.jdata['model_devi_f_trust_lo']
        hi= self.jdata['model_devi_f_trust_hi']

        data = self._model_devi(iter_idx, sys_idx , case_idx)
        
        t,p = self._get_PT(iter_idx, sys_idx , case_idx)
        
        ax.plot(data[:,0],data[:,4],'o', mfc='none',mew=0.5,ms=2)
        
        ax.axhline(lo,linestyle='--',lw=0.5,color='k')
        ax.axhline(hi,linestyle='--',lw=0.5,color='k')  
        
        ax.set_title('T = %.d K, p = %.2f GPa, '%(t,p)+self.label_list[sys_idx])

        ax.set_ylabel('model_devi ($\\rm{eV/\mathring A}$)')
        ax.set_xlabel('timestep')
        ax.set_xlim(data[0,0],data[-1,0])
        
    def _plot_sampling(self, ax, iter_idx,  sys_idx = 0, color='dodgerblue', is_label=False, label=''):
        
        self.temps, self.press, self.vol, self.energy = self._sampling(iter_idx, sys_idx)

        if is_label:
            ax.plot(self.press, self.temps, 
                'o', color=color, alpha=0.4, mew=0.5, mfc='none',ms=3, label=label)   
        else:
            ax.plot(self.press, self.temps, 
                'o', color=color, alpha=0.4, mew=0.5, mfc='none',ms=3)   
        
        if self.printf:
            output = 'Iter %.6d '%iter_idx
            output += ' add %.d frames '%self.temps.shape[0]
            output += '-- sys.%.3d '%sys_idx 
            output += '(' + self.label_list[sys_idx] +')'

            print(output)

        ax.set_xlabel('Pressure (GPa)')
        ax.set_ylabel('Temperature (K)')


    def _all_sampling_plot(self, ax, color):

        count = np.zeros(len(self.confs_list))
        
        for i in range(self.N_iter-1):

            for sys_idx in range(len(self.confs_list)):
                
                try:
                    # label only plots in the first time for single system
                    if count[sys_idx] == 0:
                        self._plot_sampling(ax, iter_idx = i, sys_idx = sys_idx, 
                                           color=color[sys_idx], is_label=True, 
                                            label='DPGEN ('+self.label_list[sys_idx]+')')
                        count[sys_idx] +=1 
                    else:
                        self._plot_sampling(ax, iter_idx = i, sys_idx = sys_idx, 
                                           color=color[sys_idx])
                except:
                    pass        

    def _collect_data_to(self, out_dir, set_numb = 10000):
        
        for sys_idx in range(len(self.label_list)):

            
            file_list = glob.glob(os.path.join(self.dir, 'iter.*', '02.fp','data.%.3d'%sys_idx))

            file_list.sort()
            
            print('Sys.%.3d is working, there are %.d sub-datasets '%(sys_idx, len(file_list)))
            
            fparam = []
            ms = dpdata.MultiSystems()
        
            for file in file_list:

                sys = dpdata.LabeledSystem( file, fmt='deepmd/raw' )
                ms.append(sys)
                
                temps = np.loadtxt( os.path.join(file,'fparam.raw') )
                fparam = np.append(fparam, temps)
            
            outdir_ss = os.path.join(out_dir, self.label_list[sys_idx])
            os.mkdir( outdir_ss)
            
            ms.to_deepmd_raw(outdir_ss)
            np.savetxt( os.path.join(outdir_ss, 'fparam.raw'), fparam)
            
            os.chdir(outdir_ss)
            os.system('mv ./*/*.raw ./')
            os.system('~/raw_to_set.sh %.d'%(set_numb))
            
            print('Sys.%.3d is done'%(sys_idx))
            

        