# PostProcess

逐步整理各类后处理的程序

`Constant.py` : Constant & Unit Trans
- 基本库的import
- 一些常用的科学常数
- 用于单位换算的常数
- 作图的一些基本函数（生成colormap，画Phonon spectra）

`Dataset.py`：训练集的分析和处理（结合dpdata使用）
- 核心是`Class SingleSys()`，通过`_get_EOS`读入压强、体积、温度和密度信息
- 数据的筛选

`DPGenFlow.py`：DPGEN过程的全分析
- 核心是`Calss DPGenSys()`，返回DPGEN的基本迭代、系统信息
- 返回训练曲线、模型偏差、相空间采样
- 返回所有采样结果在PT上的分布
- collect所有数据（包含`fparam.raw`）

`Traj.py`：处理msst的lammps轨迹，转化为vasp输入文件
- 核心是`msst_traj()`，读入两个文件，一个是fix print输出的热力学量，一个是轨迹
- 展示轨迹在热力学相空间的采样情况
- 把lammps轨迹转化为vasp的INCAR,POSCAR,以用于打标签
- collect所有数据（包含`fparam.raw`）
