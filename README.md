# PostProcess

逐步整理各类后处理的程序

`Constant.py` : Constant & Unit Trans
- 基本库的import
- 一些常用的科学常数
- 用于单位换算的常数
- 作图的一些基本函数

`Dataset.py`：训练集的分析和处理（结合dpdata使用）
- 基于`SingleSys()`进行分析，通过_get_EOS读入压强、体积、温度和密度信息
- 数据的提取和整合
- (待补充)原子环境特征空间分析

`DPGenFlow.py`：DPGEN过程的全分析
- 基于DPGenSys()，返回DPGEN的基本迭代、系统信息
- 返回训练曲线、模型偏差、相空间采样
- 返回所有采样结果在PT上的分布
