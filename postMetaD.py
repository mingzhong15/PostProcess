from Constant import *

class metaSys():
    
    def __init__(self, DIR):
        
        try:
            self.plm_data = np.loadtxt(DIR + 'plm.out')
        except:
            self.plm_data = np.loadtxt(DIR + 'out.plm')
            
    def _plot_evo_cv(self, cv1, cv2, label):
        
        
        
        