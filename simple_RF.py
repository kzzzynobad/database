import pandas as pd  # 加载pandas软件包
from sklearn.model_selection import train_test_split  # 导入相关库
from sklearn.ensemble import RandomForestRegressor  # 导入相关库
from sklearn.metrics import mean_squared_error, r2_score  # 导入相关库
import matplotlib.pyplot as plt  # 导入绘图相关库
data = pd.read_csv('YoungModulus_simple.csv')   # 读取名为YoungModulus_simple，格式为csv的文件，如果文件与代码不在同一个文件夹内，则需要修改为r'F:\pychram\YoungModulus_simple.csv',仅为举例，请类比。
# 将data切块，分成features和label
# []内按'，'分开，然后前面是行后面是列，0：表示从0至最后一行，0：6表示从第0列至第5列。单独一个数，则代表单独一行/列。注意，python的区间全部是左闭右开。
features = data.iloc[0:, 0:6]
label = data.iloc[0:, 6]

# test_size表示测试集的大小，分数(0-1),通常取0.2，random_state为随机数，确定随机数能让每次的split相同。
features_train, features_test, label_train, label_test = train_test_split(features, label, test_size=0.2, random_state=1)

model = RandomForestRegressor()  # 实例一个随机森立回归模型
model.fit(features_train, label_train)  # 将训练集的特征和标签fit进模型进行训练

label_train_pred = model.predict(features_train)  # 预测训练集
label_test_pred = model.predict(features_test)  # 预测测试集
print(r2_score(label_train, label_train_pred))  # 训练集的R2和MSE
print(mean_squared_error(label_train, label_train_pred))
print(r2_score(label_test, label_test_pred))  # 测试集的R2和MSE
print(mean_squared_error(label_test, label_test_pred))

fig, ax = plt.subplots()  # 绘图
# 训练集和测试集的散点
ax.scatter(label_train, label_train_pred, color="blue")
ax.scatter(label_test, label_test_pred, color="red")
a = 0
b = 200
# 绘制对角线
ax.plot([a, b], [a, b], "--k")
ax.set_ylabel("target predict")
ax.set_xlabel("true target")
# ax.set_title("")
# 图注
ax.text(
    25, 180,
    r"train $R^2$=%.2f,MSE=%.2f "
    % (r2_score(label_train, label_train_pred), (mean_squared_error(label_train, label_train_pred)))

)
ax.text(
    25, 160,
    r"test $R^2$=%.2f,MSE=%.2f "
    % (r2_score(label_test, label_test_pred), (mean_squared_error(label_test, label_test_pred)))

)
# x,y轴上下限
ax.set_xlim([a, b])
ax.set_ylim([a, b])
plt.show()
