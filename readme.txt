# Model----包括AR、MA、ARMA、ARIMA模型：
## AR model.m  
## MA model.m  
## ARMA model.m  
## ARIMA model.m

# Dependent function----包括AR、MA、ARMA、ARIMA模型运行所需要的函数：
## ar.m
## armax.m
## autocorr.m
## iddata.m
## parrcorr.m
## predict.m

# 运行环境：Windows 7及以上

# 运行平台：MATLAB

# 输入：AR model、MA model、ARMA model、ARIMA model中的参数‘y'均为时序数据。
# 输出：AR model、MA model、ARMA model、ARIMA model中的参数‘PF’均为时序预测值。

# 模型说明
- （1）AR：Auto Regressive Model，自回归模型
该模型是通过自身前面部分的数据与后面部分的数据之间的相关关系（自相关）来建立回归方程，该模型一般形式如下所示：

其中为白噪声，为模型系数，上述方程称为p阶的自回归模型，记为AR(p)。

- （2）MA：Moving Average Model，移动平均模型
MA模型的一般形式如下所示：

其中为白噪声，为模型系数，上述方程称为q阶移动平均模型，记为MA(q)。

- （3）ARMA：Auto Regressive and Moving Average Model，自回归移动平均模型
自回归移动平均模型是自回归和移动平均模型两部分组成的，ARMA的一般形式如下所示：
其中p是自回归阶数，q是移动平均阶数，上述模型可以记为ARMA(p,q)。

- （4）ARIMA：Auto Regressive Integrate Moving Average Model，差分自回归移动平均模型
AR/MA/ARMA模型适用于平稳时间序列的分析，而ARIMA模型能够用于齐次非平稳时间序列的分析，该模型可表示为ARIMA(p,d,q)，其中p为自回归阶数，q为移动平均阶数，d为时间成为平稳时所做的差分次数。

# 注：代码中有模型运行的具体注释。
