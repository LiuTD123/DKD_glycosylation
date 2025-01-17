from sklearn.svm import SVC  # 导入SVC类用于SVM  
import pandas as pd  
import numpy as np  
from sklearn.model_selection import train_test_split  
import pickle  
  
# 数据加载部分保持不变  
dataFile = 'D:\\xqm2\\xqm\\xqm\\2024\\9月\\0905\\5 对残差筛选出来得gene进行诊断模型构建\\SVM建模\\输入数据\\5_gene_logistic.csv'  
data = pd.read_csv(dataFile)  
X = data.iloc[:, 0:-1]  
y = data['label']  
  
# 索引和数据集分割部分保持不变  
indices = list(range(len(data)))  
X_train, X_test, y_train, y_test, train_indices, test_indices = train_test_split(X, y, indices, test_size=0.5, random_state=21, stratify=y)  
  
# 打印索引部分保持不变  
print("训练集样本索引:", train_indices)  
print("测试集样本索引:", test_indices)  
  
# 使用SVM模型  
clf_svm_linear = SVC(kernel='linear',probability=True, max_iter=1000)  # 使用线性核的SVM  
clf_svm_linear.fit(X_train, y_train)  
acc_train_svm_linear = clf_svm_linear.score(X_train, y_train)  
acc_test_svm_linear = clf_svm_linear.score(X_test, y_test)  
print(acc_train_svm_linear, acc_test_svm_linear)  
  
# 使用非线性核的SVM（可选）  
clf_svm_rbf = SVC(kernel='rbf',probability=True, max_iter=1000)  # 使用RBF核的SVM  
clf_svm_rbf.fit(X_train, y_train)  
acc_train_svm_rbf = clf_svm_rbf.score(X_train, y_train)  
acc_test_svm_rbf = clf_svm_rbf.score(X_test, y_test)  
print(acc_train_svm_rbf, acc_test_svm_rbf)  
  
# SVM模型没有coef_属性，而是有coef_和intercept_，但这里只打印coef_  
print("Linear SVM Coefficients:", clf_svm_linear.coef_)  
# 如果需要查看RBF SVM的权重，注意RBF SVM没有coef_，因为核函数将特征映射到更高维空间  
  
# 保存SVM模型到文件  
with open('svm_model_linear.pkl', 'wb') as file:  
    pickle.dump(clf_svm_linear, file)  
  
# 重新加载模型  
with open('svm_model_linear.pkl', 'rb') as file:  
    clf_svm_linear_loaded = pickle.load(file)  
  
# 打印类型以确认是模型对象  
print(type(clf_svm_linear_loaded))  
  
# 如果你还想保存和加载RBF SVM模型，可以重复上述步骤，但文件名要不同  
with open('svm_model_rbf.pkl', 'wb') as file:  
    pickle.dump(clf_svm_rbf, file)  