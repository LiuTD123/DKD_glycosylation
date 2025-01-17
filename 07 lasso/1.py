import pandas as pd
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LassoCV
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import StandardScaler


# 假设 CSV 文件的路径已经给出，并且 CSV 文件包含两列：一列是患病情况，其余列是基因表达数据
data = pd.read_csv('26_gene_exp_lasso_shuru.csv')

# 分离特征（基因表达数据）和目标变量（患病情况）
X = data.iloc[:, 1:]  # 从第二列开始，所有列都是基因表达数据
y = data.iloc[:, 0]   # 第一列是患病情况

# 划分训练集和测试集（此处只使用训练集进行 Lasso 回归）
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=21)



# 使用 LassoCV 进行交叉验证以选择最佳的正则化参数 alpha
lasso_cv = LassoCV(cv=5, random_state=21)  # cv 是交叉验证的折数
# 创建 StandardScaler 对象
scaler = StandardScaler()

# 拟合 scaler 并转换 X_train
X_train_scaled = scaler.fit_transform(X_train)

# 由于LassoCV内部会处理训练集和验证集，我们不需要对X_test进行变换
# 如果你需要使用相同的scaler变换X_test，可以调用scaler.transform(X_test)

# 现在使用标准化后的数据来训练 LassoCV
lasso_cv.fit(X_train_scaled, y_train)
# 获取最优的 alpha 值
print("Optimal alpha:", lasso_cv.alpha_)

# 获取非零系数对应的基因
selected_genes_cv = X_train.columns[lasso_cv.coef_ != 0]
print("Selected genes by LassoCV:", selected_genes_cv)

# 查看每个基因的系数
print("Coefficients of the selected genes by LassoCV:", lasso_cv.coef_)

# 绘制 LassoCV 的 MSE 图，并添加每个点对应的上下限的竖线标，横坐标使用 alpha 的对数
plt.figure(figsize=(10, 6))
mse_paths = lasso_cv.mse_path_
mean_mse = mse_paths.mean(axis=1)
mse_std = mse_paths.std(axis=1)  # 计算每个 alpha 的 MSE 的标准差

# 由于 semilogx 已经在 x 轴上应用了对数刻度，我们不需要再次取对数
plt.semilogx(lasso_cv.alphas_, mean_mse, label='Mean MSE (cross-validation)', marker='o')  # 使用 semilogx

# 添加误差棒
for alpha, mean, std in zip(lasso_cv.alphas_, mean_mse, mse_std):
    plt.semilogx([alpha, alpha], [mean - std, mean + std], 'k-', alpha=0.5)  # 上下限竖线

# 标记最优的 alpha 值
plt.axvline(lasso_cv.alpha_, linestyle='--', color='r', label='Optimal alpha')

plt.xlabel('Alpha (log scale)')
plt.ylabel('Mean Squared Error')
plt.title('Mean Squared Error vs Alpha (log scale)')
plt.legend()
plt.show()

# 绘制 Lasso 系数路径图
coef_paths = lasso_cv.path(X_train_scaled, y_train)[1].T  # 获取系数路径
n_alphas = len(lasso_cv.alphas_)
n_features = X_train_scaled.shape[1]

# 创建颜色列表，用于绘制不同的系数路径
colors = ListedColormap(['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'black'])

plt.figure(figsize=(10, 6))
for i in range(n_features):
    plt.semilogx(lasso_cv.alphas_, coef_paths[:, i], label=f'Gene {i+1}', color=colors(i))

# 由于您不希望有竖线，这里删除了 plt.axvline() 函数的调用
# plt.axvline(lasso_cv.alpha_, linestyle='--', color='black', label='alpha CV')

plt.xlabel('Alpha')
plt.ylabel('Coefficients')
plt.title('Lasso Coefficient Path')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()