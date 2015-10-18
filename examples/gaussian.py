from numpy import *
from pylab import *

d = loadtxt('gaussian.txt')

####    un-comment the following if you re-run the test using mpirun -np 4 ./gaussian_test
#d = loadtxt('gaussian_0.txt')

#for i in range(1, 10):
#    fname = 'gaussian_' + str(i) + '.txt'
#    print '--> loading ' + fname
#    dd = loadtxt(fname)
#    d = vstack((d, dd))

subplot(2,2,1)
hist2d(d[:,2], d[:,3], 50, weights=d[:,0], );

subplot(2,2,2)
hist2d(d[:,2], d[:,4], 50, weights=d[:,0], );

subplot(2,2,3)
hist2d(d[:,3], d[:,4], 50, weights=d[:,0], );

subplot(2,2,4)
hist(d[:,2], 50, weights=d[:,0], histtype='step', normed=True);
hist(d[:,3], 50, weights=d[:,0], histtype='step', normed=True);
hist(d[:,4], 50, weights=d[:,0], histtype='step', normed=True);

figure()
subplot(1,2,1)
hist(d[:,0], 30);

subplot(1,2,2)
hist(d[:,1], 50, weights=d[:,0]);

show()
