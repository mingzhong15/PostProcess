from Constant import *

class AlaSys():
    
    def __init__(self, DIR, element, latt_type, a0, mass_mole, ALAMODE_ROOT, LAMMPS, VASP=''):
        
        self.DIR = DIR
        self.element = element
        self.latt_type = latt_type
        
        if latt_type == 'fcc':
            self.natoms_uc = 4
        elif latt_type == 'bcc':
            self.natoms_uc = 2
        
        self.a0 = a0
        self.mass_mole = mass_mole
        
        print('System : '+latt_type+'-'+element+' (a0 = %.3f ang)'%a0)

        self.LAMMPS = LAMMPS
        self.VASP = VASP
        self.ALAMODE_ROOT = ALAMODE_ROOT

        self.pattern_type = 'lammps'
        
    def _generate_lammps_input(self, nx, model_type, model_file, in_type, strain = 1.0, fparam=300):
        
        output_acc = '%20.12f'
        
        self.nx = nx
        self.prefix = self.element+"%.d%.d%.d"%(self.nx,self.nx,self.nx)
        
        self.model_file = model_file
        self.model_type = model_type
        self.natoms = self.nx**3*self.natoms_uc
        self.strain = strain
        self.fparam = fparam
        
        # ======= basic ======== #
        ret=''
        ret+='units           metal\n'
        ret+='atom_style      atomic\n'
        ret+='boundary        p p p\n'
        ret+='\n'

        if in_type == 'minimize':
            # ======= geometric ====== #
            ret+='lattice         '+self.latt_type+' %.3f\n'%self.a0
            ret+='region          box block 0 1 0 1 0 1 units lattice\n'
            ret+='create_box      1 box\n'
            ret+='create_atoms    1 box\n'
            ret+='replicate       %.d %.d %.d\n'%(nx,nx,nx)
            ret+='\n'
            
        elif in_type == 'force':
            ret+='read_data       tmp.lammps\n'
            ret+='\n'
            

        # ======= potential ====== #
        if model_type == 'dp':
            ret+='pair_style      deepmd '+model_file+' fparam %.d \n'%fparam
            ret+='pair_coeff      * *\n'
            ret+='mass            1 %.2f\n'%self.mass_mole

        elif model_type == 'eam':
            ret+='pair_style      eam/alloy\n'
            ret+='pair_coeff      * * '+model_file+' '+self.element+'\n'

        ret+='\n'

        if in_type == 'minimize':
            ret+='timestep        1e-3\n'
            ret+='reset_timestep  0\n'
            ret+='\n'

            ret+='run             1\n'
            ret+='write_data      '+self.prefix+'.data\n'

            ret+='dump            1 all custom 1 coord.dat id type xu yu zu\n'
            ret+='dump_modify     1 sort id format float "'+output_acc+'"\n'
            ret+='dump            2 all custom 1 coord_scale.dat type xs ys zs\n'
            ret+='dump_modify     2 sort id format float "'+output_acc+'"\n'
            ret+='\n'

            ret+='run             0\n'
            ret+='clear           \n'
            #print(ret)

        
            file=open(os.path.join(self.DIR,'in.minimize'),'w')
            
            print('Creating in.minimize ...')
            print('%.d*%.d*%.d supercell'%(nx,nx,nx)+' (%.d atoms)'%(self.natoms))


        elif in_type == 'force':
            ret+='dump            1 all custom 1 XFSET id xu yu zu fx fy fz\n'
            ret+='dump_modify     1 format float "'+output_acc+'"\n'
            ret+='run             0\n'
    
            file=open(os.path.join(self.DIR,'in.force'),'w')
        
            print('Creating in.force ...')
            
        file.writelines(ret)
        file.close()

    def _run_minimize(self, output_poscar=False):

        os.chdir(self.DIR)
        
        os.system('mpirun -np 32 '+self.LAMMPS+' -i in.minimize')

        # =================== #
        # generate coord_scale for alm_input
        # =================== #        
        cmd = ''

        cmd += "natom=`echo %.d"%(self.natoms)+" | awk '{printf(\"%4d\",$1)}'`\n"
        cmd += "l1=`wc -l "+self.prefix+".data | awk '{printf(\"%4d\",$1)}'`\n"

        cmd += "sed -i \"16,${l1}d\" "+ self.prefix+'.data\n'

        cmd += "l3=`wc -l coord.dat | awk '{printf(\"%4d\",$1)}'`\n"

        cmd += "l4=`echo $[l3-natom] | awk '{printf(\"%4d\",$1)}'`\n"

        cmd += "sed -i \"1,${l4}d\" coord.dat\n"

        cmd += "cat "+self.prefix+".data coord.dat > "+self.prefix+".lammps\n"

        cmd += "sed -i \"1,${l4}d\" coord_scale.dat\n"

        #print(cmd)
        os.system(cmd)
        
        # =================== #
        # generate poscar
        # =================== #
        if output_poscar:
            frame = dpdata.System(self.DIR+self.prefix+'.lammps', 
                                  fmt='lammps/lmp', type_map=[self.element])
            frame.to('vasp/poscar', self.DIR+self.prefix+'.POSCAR')

    def _generate_alm_input(self, calc_type, k_path='', extra=''):

        outfile = 'alm_'+calc_type+'.in'
        self.alm_type = calc_type
        
        prt = '=======================================\n'
        prt += 'generating alamode input file for :'+self.alm_type+' calculation \n'
        print(prt)                  
        
        ret = ''
        ret += "&general\n"
        if calc_type == 'pattern':
            ret += "  PREFIX = "+self.prefix + "\n"
            ret += '  MODE = suggest\n'
            ret += '  NAT = %.d\n'%self.natoms
            
        elif calc_type == 'IFC':
            ret += "  PREFIX = "+self.element+"_harmo\n"
            ret += '  MODE = optimize\n'
            ret += '  NAT = %.d\n'%self.natoms
            
        elif calc_type == 'phonon' or calc_type == 'phdos' :
            ret += "  PREFIX = "+self.element+"_harmo\n"
            ret += '  MODE = phonons\n'
            ret += '  MASS = %.2f\n'%self.mass_mole
            ret += '  FCSXML = '+self.element+'_harmo.xml\n'
            
            
        ret += extra
        
        ret += '  NKD = 1\n'
        ret += '  KD = '+self.element+'\n'
        ret += '/\n'
        ret += '\n'

        if calc_type == 'IFC':
            ret += "&optimize\n"
            ret += "   DFSET = displace/DFSET_harmonic\n"
            ret += '/\n'
            ret += '\n'
        else:
            pass

        if calc_type == 'phonon' or calc_type == 'phdos':
            pass
        else:
            ret += "&interaction\n"
            ret += "  NORDER = 1\n"
            ret += "  NBODY = 2\n"
            ret += "/\n"
            ret += "\n"

        
        if calc_type == 'phdos':
            ret += "&analysis\n"
            ret += "  PRINTMSD = 1\n"
            ret += "/\n"
            ret += "\n"
        
        ret += "&cell\n"
        ret += "  %.11f\n"%A2bohr   
        
        if calc_type == 'phonon' or calc_type == 'phdos' :
            # primitive cell
        
            if self.latt_type == 'fcc':
                lpc = self.a0 * 0.5
                ret +=  '  0.0 %.6f %.6f\n'%(lpc,lpc)
                ret += '  %.6f 0.0 %.6f\n'%(lpc,lpc)
                ret += '  %.6f %.6f 0.0\n'%(lpc,lpc)
                
            elif self.latt_type == 'bcc':
                lpc = self.a0 * 0.5
                ret += '  %.6f   %.6f  -%.6f \n'%(lpc,lpc,lpc)
                ret += ' -%.6f   %.6f   %.6f \n'%(lpc,lpc,lpc)
                ret += '  %.6f  -%.6f   %.6f \n'%(lpc,lpc,lpc)
                
            ret += "/\n"
            ret += "\n"
 
            ret += "&kpoint \n"
            ret += k_path
            ret += '/\n'
    
            file = open( os.path.join(self.DIR, outfile), 'w' )
            file.writelines(ret)
            file.close()
    
        else:
            ret += "  %.6f 0.0 0.0 \n"%(self.a0*self.nx*self.strain)
            ret += "  0.0 %.6f 0.0 \n"%(self.a0*self.nx)
            ret += "  0.0 0.0 %.6f \n"%(self.a0*self.nx)
            
            ret += "/\n"
            ret += "\n"


            ret += "&cutoff \n"
            ret += "  *-* None\n"
            ret += "/\n"
            ret += "\n"

            ret += "&position \n"
            #print(ret)

            #file = open( os.path.join(ss.dir, 'alm0temp.in'), 'w' )
            file = open( os.path.join(self.DIR, outfile), 'w' )
            file.writelines(ret)

            coord = open( os.path.join(self.DIR, 'coord_scale.dat'), 'r' )
            lines = coord.readlines()

            file.writelines(lines)
            file.writelines('/\n')

            file.close()
            coord.close()

    def _run_alm(self, suffix=''):
        
        os.chdir(self.DIR)

        alm_file = 'alm_'+self.alm_type+'.in'
        
        prt = '=======================================\n'
        prt += 'performing : '+self.alm_type+' calculation ('+alm_file +')\n'
        print(prt)             
        
        if self.alm_type == 'phonon' or self.alm_type == 'phdos':
            os.system(self.ALAMODE_ROOT +'/anphon/anphon '+alm_file+' >> alm.log')
            
            if suffix != '':
                if self.alm_type == 'phonon':
                    
                    output_pre = self.element+'_harmo.bands'
                    output = output_pre + '.'+suffix
                    os.system('mv '+output_pre+' '+output)
                    
                    #output_pre = self.element+'_harmo.thermo'
                    #output = output_pre + '.'+suffix
                    #os.system('mv '+output_pre+' '+output)
                    
                elif self.alm_type == 'phdos':
                    output_pre = self.element+'_harmo.dos'
                    output = output_pre + '.'+suffix
                    os.system('mv '+output_pre+' '+output)
                    
                    output_pre = self.element+'_harmo.msd'
                    output = output_pre + '.'+suffix
                    os.system('mv '+output_pre+' '+output)
                    
                    output_pre = self.element+'_harmo.thermo'
                    output = output_pre + '.'+suffix
                    os.system('mv '+output_pre+' '+output)
                    
        else:
            cmd = self.ALAMODE_ROOT +'/alm/alm '+alm_file+' >> alm.log'
            print(cmd)
            os.system(cmd)
            
        if self.alm_type == 'pattern':

            os.system('mkdir displace')
            os.chdir( os.path.join(self.DIR,'displace') )

            INPUT = self.ALAMODE_ROOT+'/tools/displace.py'
            cmd = 'python ' + INPUT 

            if self.pattern_type == 'lammps':      
                cmd += ' --LAMMPS ../'+self.prefix+'.lammps'
            elif self.pattern_type == 'vasp':
                cmd += ' --VASP ../'+self.prefix+'.POSCAR'

            cmd += ' --prefix harmo --mag 0.01 '
            cmd += ' ../'+self.prefix+'.pattern_HARMONIC >> run.log \n'
            cmd += '\n'

            if self.pattern_type == 'lammps':      
                cmd += 'cp ../'+self.model_file+' .\n'
                cmd += 'cp ../in.force .\n'
                cmd += '\n'

            print(cmd)

            os.system(cmd)   
        
            
    def _cal_pattern(self, np = 28, sbatch_title=''):

        prt = '=======================================\n'
        prt += 'Force calculation : '+self.pattern_type+'\n'
        print(prt)     
        
        if self.pattern_type == 'lammps':   

            os.chdir( os.path.join(self.DIR,'displace') )

            harmo_list = glob.glob(self.DIR+'/displace/harmo*.lammps')
            self.npattern = len(harmo_list)

            for idx in range(self.npattern):

                if self.npattern < 10:
                    cmd = 'cp harmo%01d.lammps tmp.lammps\n'%(idx+1)
                elif self.npattern >= 10:
                    cmd = 'cp harmo%02d.lammps tmp.lammps\n'%(idx+1)
                cmd += 'mpirun -np %.d '%np+self.LAMMPS+' -i in.force > run.log\n'
                cmd += 'mv XFSET XFSET.harmo%01d\n'%(idx+1)
                
                print(cmd)
                os.system(cmd)
                
        elif self.pattern_type == 'vasp':

            os.chdir( os.path.join(self.DIR,'displace') )    

            harmo_list = glob.glob(self.DIR+'/displace/harmo*.POSCAR')
            self.npattern = len(harmo_list)

            for idx in range(self.npattern):

                pp = 'harmo%01d'%(idx+1)

                os.system('mkdir '+pp)
                os.chdir( os.path.join(self.DIR,'displace',pp) )

                cmd = 'cp ../../INCAR ./\n'
                cmd +='cp ../../POTCAR ./\n'
                cmd +='cp ../'+pp+'.POSCAR POSCAR\n'

                cmd +='sed -i "s/xx/%.11f/g" INCAR\n'%(self.fparam*kb*J2eV)

                print(cmd)
                os.system(cmd)
                
                sbatch_title += 'mpirun --np %d '%np+self.VASP+'\n'
                sbatch_title +='mv vasprun.xml ../vasprun_'+pp+'.xml\n'

                print(sbatch_title)

                file=open('./job.sbatch','w')
                file.writelines(sbatch_title)
                file.close()
                os.system('sbatch < job.sbatch')

                #os.system(cmd)     
        
    def _collect_pattern(self):
        
        os.chdir( os.path.join(self.DIR,'displace'))
        
        prt = '=======================================\n'
        prt += 'Extracing force results : '+self.pattern_type+'\n'
        print(prt)
        
        INPUT = self.ALAMODE_ROOT+'/tools/extract.py'
        cmd = 'python ' + INPUT 
        
        if self.pattern_type == 'lammps':
            cmd += ' --LAMMPS ../'+self.prefix+'.lammps'
            cmd += ' XFSET.harmo* > DFSET_harmonic\n'
        elif self.pattern_type == 'vasp':
            cmd += ' --VASP ../'+self.prefix+'.POSCAR'
            cmd += ' vasprun_harmo*.xml > DFSET_harmonic\n'
            
        print(cmd)

        os.system(cmd)