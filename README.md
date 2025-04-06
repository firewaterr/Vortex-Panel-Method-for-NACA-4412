# Vortex-Panel-Method-for-NACA-4412
coding in MATLAB
# 代码说明
本代码分为两部分。'AirfoilNACA4412.m'文件将弦长等分为100个点，构建出了NACA4412翼型。直接运行可得到翼型效果图。'VortexPanelMethod.m'文件运用'AirfoilNACA4412.m'处理好的翼型点，计算了攻角为0°，9°两种情形中翼型上各点压力系数分布，直接运行可得到结果。所有的结果图都已导出并存储在'figure'文件夹中。
# 代码精度
本代码为一阶精度格式。与'xfoil'图像对比发现，压力系数在翼型前缘及中部大部分区域内差别不大。最初计算得到的面涡强度不稳定，将面涡强度两两平均之后，得到的结果较好。而且可以发现，压力系数在翼型后缘几个点上仍旧非常不稳定，这与代码精度较低有关。