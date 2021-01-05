# 文件目录说明
- frontend 是前端
- backend 是后端
- utils 是工具
- notes 是笔记
有啥需要的随便加

## 编译说明：

```c++
cd BA_schur  \\ 进入文件夹
mkdir build   
cd build
cmake ..
make -j4    
```

## 开发记录
```c++
1. ubuntu 16.04 列清楚数据和整个架构
2. ubuntu 18.04 继续理清思路

```

```
1.shared_ptr: hared_ptr是一种智能指针（smart pointer），作用有如同指针，但会记录有多少个shared_ptrs共同指向一个对象。这便是所谓的引用计数（reference counting）
2. 
```

| camera pose |相机观测值 Pc| points| 
|---|---|---|
|  相机位姿，3 个 pose|featurePerId |特征点，数目是20|
|元素包括Rwc qwc twc,featurePerId|Pw经过平移和旋转投影到归一化平面上Pc|Pw,都被3帧frame观测到| 
 
 


| problem  | problem|vertexCams_vec|allPoints|verticies_|
|---|---|---|---|---|
|AddVertex : pose 和逆深度|AddEdge: pose和逆深度|pose的7个量|逆深度|顶点数:size 23 前面3个6维的pose 后面20个是1维的逆深度|
|顶点：优化的变量 pose和观测值Pc|边：重投影误差项|四元素和平移量|Pc的倒数|所有的观测：逆深度和相机pose|
