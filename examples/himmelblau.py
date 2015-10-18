from numpy import *
from pylab import *

print 'loading chain ...'
chain = loadtxt('himmelblau.txt')

print 'chain size:', chain.shape


figure()
hist(chain[:,1], 50);
xlabel('$\\chi^2$')

figure()
hist2d(chain[:,2], chain[:,3], 100, weights=chain[:,0], cmin=0)
colorbar()

show()
