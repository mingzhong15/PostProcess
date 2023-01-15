import cmaps
from Constant import *

class MetadSys():

    def __init__(self, meta_dir):
        
        self.dir = meta_dir
        
        self.cv_file = os.path.join(self.dir ,'COLVAR')
        self.hill_file = os.path.join(self.dir, 'HILLS')

        self.label_list = os.popen('head -n 1 '+ self.cv_file).readlines()[0].split()[2:]

        
    def _plot_cv_evol(self):
        
        data = np.loadtxt(self.cv_file)
        
        num_cv = int(data.shape[1]-3)
        
        fig, ax1 = plt.subplots(num_cv,1,figsize=(3,1.5*num_cv),dpi=200)
        
        for i in range(num_cv):    
            ax = ax1[i]

            ax.plot(data[:,0]/1e3, data[:,1+i])
            ax.set_ylabel(self.label_list[1+i])

        ax.set_xlabel(self.label_list[0])
        
        return ax1
        
    def _cv_scatter(self):
        
        data = np.loadtxt(self.cv_file)
        
        fig, ax = plt.subplots(figsize=(3,2),dpi=200)

        cs=ax.scatter(data[:,1], data[:,2], c=data[:,0], s=2,  
                      cmap=plt.cm.Spectral_r)

        plt.colorbar(cs, label='Bias')

        ax.set_xlabel(self.label_list[1])
        ax.set_ylabel(self.label_list[2])
        
        return ax

    def _bias_decay(self):

        try:
            data = np.loadtxt(self.hill_file)
            
            fig, ax = plt.subplots(figsize=(4,1.5),dpi=200)

            ax.plot(data[:,0]/1e3, data[:,-2],'-o',ms=2,lw=0.5)

            ax.set_xlabel('timestep (ns)')
            ax.set_ylabel('Bias')
        except:
            pass
        
        
    def _get_fes(self, n):
        
        self.ny = int(os.popen('grep nbins_'+self.label_list[1]+' '+self.dir+'fes_0.dat').readlines()[0].split()[-1])
        self.nx = int(os.popen('grep nbins_'+self.label_list[2]+' '+self.dir+'fes_0.dat').readlines()[0].split()[-1])

        data = np.loadtxt(self.dir + 'fes_%.d.dat'%n)
    
        x = data[:,0].reshape(self.nx, self.ny)
        y = data[:,1].reshape(self.nx, self.ny)
        z = data[:,2].reshape(self.nx, self.ny)
        
        return x,y,z

    def _plot_fes(self, x,y,z):
                    
        fig, ax = plt.subplots(figsize=(3,2),dpi=200)
        
        cs_line = ax.contour(x, y, -z, 15, linestyles='-', vmin=0, 
                         linewidths=0.1, colors='white')
        ax.clabel(cs_line, inline=5, fontsize=5, colors='white')
    
        cs = ax.contourf(x, y, -z, 40, alpha=1, cmap=cmaps.cmp_b2r)

        ax.set_xlabel(self.label_list[1])
        ax.set_ylabel(self.label_list[2])
        plt.colorbar(cs, label="Free Energy (eV)") #/atom

        plt.subplots_adjust(hspace=0.2, wspace=0.3)
        
        return ax
        