import numpy as np
import matplotlib.pyplot as plt

x= np.linspace(0,100,100)

y= .5*np.random.randn(len(x))

plt.plot(x,y,color='k')
plt.ylim(-5,5)
plt.axis('off')
plt.show()
