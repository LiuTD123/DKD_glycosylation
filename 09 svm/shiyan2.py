import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn.preprocessing import StandardScaler
from pickle import load
import pickle


# 配置Matplotlib库的默认参数
plt.rcParams['font.family'] = 'SimHei'  # 中文
plt.rcParams['figure.dpi'] = 150  # 绘图精细程度
plt.rcParams['savefig.dpi'] = 600  # 保存的精细程度
plt.rcParams['axes.unicode_minus'] = False  # 确定负号的正确显示

# 读取数据
df = pd.read_excel(r"D:\\xqm2\\xqm\\xqm\\2024\\9月\\0906\\2 Z008-L002模型评估\\ROC+MSE综合下来觉得 SVM 0.2 42 linear效果最好\\SVM 0.2 42 linear可视化\\test_size 0.2 random_state 42 linear\\train\\5_gene_logistic_train.xlsx", index_col=0)

# 选择特征和目标变量
X = df.iloc[:, 1:6]  # 从第7列到第8列的数据作为特征
y = df.iloc[:, -1]   # 第一列的数据作为目标变量

# 转换为NumPy数组
X = X.values
y = y.values

# 数据分割
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, shuffle=True ##这里的分割在这个代码中没有起作用
)

# 载入之前训练好的SVM模型
#with open('path_to_your_model.pkl', 'rb') as file:  # 替换为您的模型文件路径
#    clf_poly_loaded = load(file)

# 定义绘制SVM决策边界的函数
def plot_svc_decision_function(clf, ax=None, plot_support_vectors=True):
    if ax is None:
        ax = plt.gca()
    
    if plot_support_vectors:
        ax.scatter(
            clf.support_vectors_[:, 0],
            clf.support_vectors_[:, 1],
            facecolors='none',
            edgecolors='k',
            s=100,
            label='Support Vectors'
        )
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    xx = np.linspace(xlim[0], xlim[1], 30)
    yy = np.linspace(ylim[0], ylim[1], 30)
    YY, XX = np.meshgrid(yy, xx)
    Z = clf.decision_function(np.c_[XX.ravel(), YY.ravel()])
    Z = Z.reshape(XX.shape)
    ax.contour(XX, YY, Z, colors='k', levels=[-1, 0, 1], alpha=0.5, linestyles=['--', '-', '--'])
    ax.legend()

#with open('path_to_your_model.pkl', 'wb') as file:
#    pickle.dump(clf_poly, file)
#with open('path_to_your_model.pkl', 'rb') as file:
#    model = pickle.load(file)
#    print(model)
    
with open('svm_model_linear.pkl', 'rb') as file:
    clf_poly_loaded = pickle.load(file)
#print(type(clf_poly_loaded))



#with open('path_to_your_model.pkl', 'rb') as file:  # 替换为您的模型文件路径
#    clf_poly_loaded = load(file)
# 创建一个三维图形，用于可视化SVM分类器的决策边界
# 使用加载的模型进行三维散点图的绘制
f2 = plt.figure()#创建一个三维的scatter plot图，用于可视化SVM分类器的决策边界
ax2 = plt.subplot(projection='3d')
ax2.view_init(elev=0,azim=30)
ax2.scatter3D(X_train[y_train == 0, 0], #使用scatter3D函数在三维空间中绘制了三个散点图
              X_train[y_train == 0, 1],
              clf_poly_loaded.decision_function(X_train[y_train == 0, :]),
              c='blue', label="Control", edgecolors='k')#第一个散点图（蓝色）表示y_train标签为0的样本点，其Z轴坐标由SVM的decision_function计算得出
ax2.scatter3D(X_train[y_train == 1, 0],
              X_train[y_train == 1, 1],
              clf_poly_loaded.decision_function(X_train[y_train == 1, :]),
              c='orange', label="ND", edgecolors='k')#第二个散点图（橙色）表示y_train标签为1的样本点，同样使用decision_function计算Z轴坐标
ax2.scatter3D(
    clf_poly_loaded.support_vectors_[:, 0], #X
    clf_poly_loaded.support_vectors_[:, 1], #Y
    clf_poly_loaded.decision_function(clf_poly_loaded.support_vectors_), #Z
    facecolor='None',
    edgecolor='k',
    s=100,
    label='Support vector'#第三个散点图表示支持向量，使用clf_poly.support_vectors_获取支持向量的坐标，并用decision_function计算其在决策边界上的位置
)
ax2.legend()
plt.show()