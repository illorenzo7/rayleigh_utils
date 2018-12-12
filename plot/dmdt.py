import matplotlib.pyplot as plt
import numpy as np
import sys, os

dirname = sys.argv[1]
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'

if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

dmdt_up, dmdt_down = np.load(datadir + 'dmdt.npy')
rr,tt,cost,sint,rr_depth,ri,ro,d=np.load(datadir + 'grid_info.npy')
mshell,mshell_msun = np.load(datadir + 'Mshell.npy')

f = 86400.*365/mshell # converts g/s --> Mshell/yr

plt.plot(rr/1.e8, f*dmdt_up, 'r', label='upflow')
plt.plot(rr/1.e8, f*dmdt_down, 'b', label='downflow')
plt.plot(rr/1.e8, f*(dmdt_up-dmdt_down),'k', label='total')

# Mark 5% line
rr5 = (ro - .05*d)/1.e8
ax = plt.gca()
ymin,ymax = ax.get_ylim()
plt.plot(np.ones(100)*rr5, np.linspace(ymin,ymax,100),'k--')
plt.ylim((ymin,ymax))

plt.legend()

plt.xlabel('r (Mm)')
plt.ylabel('dM/dt (shell masses/yr)')
plt.title ('Spherically-integrated mass flux')
plt.xlim((ri/1.e8, ro/1.e8))
plt.savefig(plotdir + 'dmdt.png', dpi=300)
plt.close()
