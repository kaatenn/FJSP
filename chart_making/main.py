import json
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# 读取 JSON 文件
with open('../cmake-build-debug/result.json', 'r') as file:
    data = json.load(file)

# 提取数组数据
myArray = data

# 创建横坐标（迭代次数）和纵坐标（工作流程时长）的列表
x_values = [i + 1 for i in range(len(myArray))]
y_values = myArray

# 设置中文支持的字体（例如SimHei）
font_properties = FontProperties(fname='Deng.ttf', size=12)

# 绘制折线图
plt.plot(x_values, y_values, marker='o', linestyle='-')

# 添加标题和标签
plt.title('工作流程时长随迭代次数的变化', fontproperties=font_properties)
plt.xlabel('迭代次数', fontproperties=font_properties)
plt.ylabel('工作流程时长', fontproperties=font_properties)

plt.savefig('../cmake-build-debug/result.png')
plt.close()
