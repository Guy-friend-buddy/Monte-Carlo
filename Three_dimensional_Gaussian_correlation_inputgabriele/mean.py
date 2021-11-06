import matplotlib.pyplot as plt
import numpy as np
mean = np.loadtxt("fort.10")
plt.title('Mean flow profile')
plt.ylabel(r"$U(y)$, $U''(y)$")
plt.xlabel(r"$y$")
plt.tick_params(axis='x',direction='in')
plt.tick_params(axis='y',direction='in')
plt.xticks(np.arange(0,18,2))
plt.yticks(np.arange(-0.2,1.2,.2))
plt.plot(mean[:,0],mean[:,1],'b-',label=r"$U(y)$")
plt.plot(mean[:,0],mean[:,2],'r-',label=r"$U''(y)$")
plt.legend()
plt.savefig("mean.png")
plt.show(block=False)