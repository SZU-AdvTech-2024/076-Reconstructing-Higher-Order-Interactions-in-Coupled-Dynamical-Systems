本项目是论文Reconstructing higher-order interactions in coupled dynamical systems的复现与拓展。
论文提供了部分代码：https://github.com/LValentinaGambuzza/Code-for-Reconstructing-higher-order-interactions-in-coupled-dynamical-systems
本人进行一些修改与添加，代码功能如下：
|代码文件 |	功能说明 |
|---------|---------|
|lasso.m（提供）|套索Lasso方法的坐标下降发实现|
|LVeqComplete.m（提供）	|Lotka-Volterra模型方程的代码实现|
|LVmodel.m（提供）|绘制Fig.1-b，7个种群的丰度随时间的演化过程|
|LVexact.m（提供）|绘制Fig.1-c，在导数已知时OLS的重构误差随采样次数增加的变化|
|LVapproximations.m（提供）|绘制Fig.1-d，在导数未知时使用不同阶有限差分近似对OLS重构误差的影响|
|ZackaryNet.mat（提供）|数值模拟2使用的数据，包括邻接矩阵、一阶单纯形和二阶单纯形|
|roessler_hoi.m（提供）|Rossler模型方程的代码实现，包括一阶单纯形和二阶单纯形的作用|
|signal_lasso.m（提供）|signal lasso方法的代码实现，采用坐标下降法|
|RosslerExact.m（提供）|绘制Fig.2-b，在导数已知时不同方法随采样次数增加对重构误差的影响|
|RosslerApprox.m（提供）|绘制Fig.2-c，在导数未知、采用有限差分近似时不同方法随采样次数增加对重构误差的影响|
|signal_lasso_2.m（添加）|使用CVX包实现signal lasso方法|
|RosslerExact_noise.m（添加）|在Rossler模型中添加均匀分布的随机量，研究不同强度的噪声对重构误差的影响|
|RosslerExact_noise_2.m（添加）|在Rossler模型中添加均匀分布的随机量，研究在统一强度的噪声下重构误差随采样次数的变化|
|NNLS_cvx.m（添加）|使用CVX包实现NNLS方法|
|lasso\_2.m（添加） | 使用CVX包实现Lasso方法|
