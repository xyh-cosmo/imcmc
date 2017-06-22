from numpy import *
from pylab import *

print 'loading chain ...'
chain = loadtxt('rosenbrock_1.txt')

print 'chain size:', chain.shape

figure()
hist(chain[:,1], 50);
xlabel('$\\chi^2$')

figure()
hist2d(chain[:,2], chain[:,3], 200, weights=chain[:,0],normed=True,cmin=0.04)
xlim(-4,6)
ylim(-1,30)
xlabel(r'$x_1$',fontsize=16)
ylabel(r'$x_2$',fontsize=16)
colorbar()

show()
