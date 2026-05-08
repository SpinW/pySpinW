import matplotlib.pyplot as plt
import numpy as np

from pyspinw.util import colormap_alpha_plot

x = np.linspace(0, 2*np.pi, 101)
y = np.sin(x)
intensity = y**2

colormap_alpha_plot(x,y,intensity, color_or_colormap='jet')
plt.xlim([0, 2*np.pi])
plt.ylim([-1, 1])

plt.show()