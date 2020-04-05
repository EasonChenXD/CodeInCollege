from sympy import *
import numpy as np
from scipy.integrate import quad, dblquad, nquad
import matplotlib.pyplot as plt
from pylab import mpl
from matplotlib.font_manager import FontProperties

mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.sans-serif'] = ['SimHei']
mpl.rcParams['axes.unicode_minus'] = False
import warnings

warnings.filterwarnings('ignore')


def gauss_two(x):
    return 1 / np.sqrt(2 * np.pi) * np.exp(-x ** 2 / 2)


theta = np.linspace(-2, 2, 100)
a = [0.7, 0.8, 1.2, 1.6, 2.0, 2.4]
a = list(reversed(a))
b = [-1, -0.89, -0.78, -0.68, -0.59, -0.47]


def gauss_list(a, b, theta):
    y = []
    for the in theta:
        res, _ = quad(gauss_two, -np.inf, a * (the - b))
        # if a * (the - b) == 0:
        print(res)
        y.append(res)

    return y


# Figure 并指定大小
# plt.figure(figsize=(8,8))
for i in range(len(a)):
    plt.plot(theta, gauss_list(a[i], b[i], theta), linewidth=1.2, linestyle='-')

for i in range(len(a)):
    plt.plot([b[i], b[i]], [0, 0.5], linewidth=1.4, linestyle="--")

for i in range(len(a)):
    plt.plot([-2, b[i]], [0.5, 0.5], linewidth=1.4, linestyle="--")

x_ticks = np.linspace(-4, 4, 9)
y_ticks = np.linspace(0, 1, 5)

# x_ticks = list(x_ticks
# x_ticks = sorted(x_ticks)
plt.xticks(x_ticks)
plt.yticks(y_ticks)

# 设置 x，y 轴的范围以及 label 标注
plt.xlim(-2, 2)
# plt.ylim(-2,2)
plt.xlabel('theta')
plt.ylabel('p(theta)')
# plt.title('正态双卵模型',fontproperties=FontProperties(fname='/System/Library/Fonts/PingFang.ttc'))
ax = plt.gca()
ax.grid(True)
# plt.set(style="darkgrid")
mpl.rcParams['font.sans-serif'] = ['SimHei']
mpl.rcParams['axes.unicode_minus'] = False
# ax.set_title('正态双卵模型',fontsize=14,color='r')
ax.set_title('正态双卵模型', fontsize=50, fontproperties=FontProperties(fname='/System/Library/Fonts/PingFang.ttc'))
plt.legend(['a=2.4 b=0.05', 'a=2.0 b=0.17', 'a=1.6 b=0.35', 'a=1.2 b=0.56', 'a=0.8 b=0.74', 'a=0.7 b=0.88'])
plt.show()
# plt.savefig("zhengtai.png", dpi=300)
