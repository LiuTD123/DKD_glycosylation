import numpy as np    
import pandas as pd    
from sklearn.decomposition import PCA    
from mpl_toolkits.mplot3d import Axes3D  # 导入Axes3D来创建三维坐标轴  
import matplotlib.pyplot as plt  
  
# 假设您已经有一个名为df的pandas DataFrame，它包含已经标准化处理过的数据集  
df = pd.read_csv('hebing_tissue_shunxu2_zhuanzhi.csv')  # 使用这行代码从CSV文件读取数据  
  
# 根据列的范围设置component列  
df['component'] = 'Unknown'  # 初始化一个默认值的列  
df.loc[0:68, 'component'] = 'A'  # 1-10列设置为A  
df.loc[69:129, 'component'] = 'B'  # 11-22列设置为B  
#df.loc[91:151, 'component'] = 'C'  # 23-43列设置为C  
#df.loc[294:353, 'component'] = 'D'  # 44-54列设置为D  
  
# 选择要进行PCA的列（假设所有数值列都用于PCA）  
# 注意：这里我们假设除了'component'列之外的所有列都是数值列  
pca_data = df.drop('component', axis=1)  
  
# 定义PCA对象，指定要保留的主成分数量  
pca = PCA(n_components=3)  
  
# 对标准化后的数据进行PCA  
pca_results = pca.fit_transform(pca_data)  
  
# 将PCA结果转换为pandas DataFrame  
pca_df = pd.DataFrame(pca_results, columns=['Principal Component 1', 'Principal Component 2', 'Principal Component 3'])  
  
# 将'component'列添加回DataFrame，以便在散点图中使用  
pca_df['component'] = df['component']  
  
# 创建一个颜色映射，为每个组分分配一个不同的颜色  
component_colors = {'A': 'red', 'B': 'blue'}  
# 绘制散点图，使用不同的颜色和标记形状  
#plt.figure(figsize=(10, 8))  # 增大图形大小  
#for component, color in component_colors.items():  
#    subset = pca_df[pca_df['component'] == component]  
#    plt.scatter(subset['Principal Component 1'], subset['Principal Component 2'], subset['Principal Component 3'],  
#                c=color, label=component, alpha=0.2)  # 设置透明度  )[component_colors.index(color)])  # 使用不同的标记形状  
  
# 添加图例、坐标轴标签和标题  
#plt.xlabel('Principal Component 1')  
#plt.ylabel('Principal Component 2')
#plt.zlabel('Principal Component 3')  
#plt.title('PCA Results with Component Groups')  
#plt.legend()  # 显示图例  
#plt.grid(True)  
#plt.show()

# 创建一个新的figure和Axes3D对象  
fig = plt.figure(figsize=(10, 8))  
ax = fig.add_subplot(111, projection='3d')  # 添加一个三维坐标轴  
  
# 绘制三维散点图  
for component, color in component_colors.items():    
    subset = pca_df[pca_df['component'] == component]    
    ax.scatter(subset['Principal Component 1'], subset['Principal Component 2'], subset['Principal Component 3'],    
               c=color, label=component, alpha=0.2)  # 注意这里不再需要plt.scatter  
  
# 添加图例、坐标轴标签和标题  
ax.set_xlabel('Principal Component 1')    
ax.set_ylabel('Principal Component 2')  
ax.set_zlabel('Principal Component 3')  # 使用Axes3D对象的set_zlabel方法  
ax.set_title('PCA Results with Component Groups')    
ax.legend()  # 显示图例  
ax.grid(True)  # 注意这可能会使三维图看起来有些混乱，通常不建议在三维图中添加网格  
plt.show()