# My_SFM
 yang_xy's SFM test


Windons10 VS2019

requires:
Eigen (my version 3.2.10)
OpenCV (my version 3.4.5)
Ceres

## 2020.02.01
## Version1.0
### 之前功能：
1. 手动输入图片之间的关系（R,c,由E在OpenCV中分解得到）计算对极线
2. 利用对极线找图片之间的对应点;
3. 实现sfm：固定光心，或者给出光心的初值较好的情况下可以做出SFM的较好的结果
4. 仅能实现一张图作为参考，剩下两张图做他的neighbor的效果；

### 更新-------------------
### 添加功能：
1. 手动输入图片的对应点和匹配关系（都不在本程序的计算范围），得到图片中许多特征点之间的匹配关系
2. 利用这些匹配关系做3D点的初值以及各个相机的平移初值（从E的分解方向，X做BA得到）
3. 利用这些良好的初值整体BA得到场景中特征点的稀疏重建；
4. 利用重投影误差判断某个X重建的是否合理；
5. 移除不合理的特征点

### 待实现：
1. 多线程或者CUDA优化，在利用对极几何找特征点、遍历特征点找对应关系（structure的初始）、BA的时候都可以优化；
2. 一些滤波效果
3. 多张图分别作rref，不只是图1做ref而其他图仅做neighbor的重建，这样可以重建更多的点
3. 自己找neighbor而不仅仅是给定neighbor'

![](https://github.com/onepiece1997426/My_SFM/raw/master/图片1.png)

