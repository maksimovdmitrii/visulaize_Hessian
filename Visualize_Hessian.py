#!/usr/bin/python3
import os
import sys
import shutil
import numpy as np
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Slider

from ase.visualize import view
from ase import Atom, Atoms
from ase.io import read, write, Trajectory
from ase.visualize.plot import plot_atoms

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-t", "--traj", dest="trajectory", help="Trajectory of relaxation")
(options, args) = parser.parse_args()



def separate_hessians(hes_file):
    
    folder = 'Hessian_progression'
    if not os.path.exists(folder):
        os.mkdir(folder)
    else:
        shutil.rmtree(folder)
        os.mkdir(folder)
    
    with open(hes_file) as hes:
        lines = hes.read()
        hessians = lines.split('Hessian\n')
        for i in range(len(hessians)):
            if len(hessians[i]) > 0:
                with open(os.path.join(folder, '{:05d}_hess.csv'.format(i)), 'w') as hh:
                    hh.write(hessians[i])


def read_multiple_csv(folder):
    
    dataframes = []
    for i in sorted(os.listdir(folder)):
        df = pd.read_csv(os.path.join(folder, i), header=None, sep='\s+')
        dataframes.append(df)
    return dataframes
        

folder = 'Hessian_progression'
hes_file = 'optrun_hessian.dat'
separate_hessians(hes_file)
dataframes = read_multiple_csv(folder)
hes_series = pd.concat(dataframes, keys=[str(i+1) for i in range(len(dataframes))])
MIN=hes_series.values.min()
MAX=hes_series.values.max()

traj = Trajectory(options.trajectory)
          

fig, ax = plt.subplots(2,2)
plt.subplots_adjust(left=0.25, bottom=0.25)

hes = 1
step = 1
l = ax[0,0].imshow(hes_series.loc[str(hes)].values, vmin = MIN, vmax = MAX)
#ax[0].margins(x=0)

axcolor = 'lightgoldenrodyellow'
steps = plt.axes([0.25, -0.0, 0.65, 0.03], facecolor=axcolor)
relaxatiostep = Slider(steps, 'Step', 1, len(dataframes), valinit=int(hes), valstep=step)

#ax[0,0].set_aspect('equal')
#ax[0,1].set_aspect('equal')
#ax[1,0].set_aspect('equal')
#ax[1,1].set_aspect('equal')
#elf.ax.set_xlim(0, self.w)
#self.ax.set_ylim(0, self.h)

def fmax(atoms):
    forces=atoms.get_forces()
    return np.sqrt((forces ** 2).sum(axis=1).max())

energies = [traj[i].get_potential_energy() for i in range(len(traj))]
forces = [fmax(traj[i]) for i in range(len(traj))]
#print(max(traj[0].get_forces()))
#sys.exit(0)

#print(forces)
ax[1,0].scatter(list(range(len(traj))), energies, c ='k')
ax[1,0].set_ylim(min(energies), max(energies))
ax[1,1].scatter(range(len(traj)), forces, c ='k')

#ax[0,0].set_aspect('equal')
#ax[0,1].set_aspect('equal')
#ax[1,0].set_aspect('equal')
#ax[1,1].set_aspect('equal')
           
def update(val):
    RelaxStep = relaxatiostep.val
    
    Step = int(RelaxStep)
    
    ax[0,0].imshow(hes_series.loc[str(Step)].values, vmin = MIN, vmax = MAX)
    ax[0,1].cla()
    plot_atoms(atoms=traj[Step],  ax=ax[0,1], rotation=('0x,0y,0z'), offset=(0, 0))
    ax[1,0].scatter(list(range(len(traj))), energies, c ='k')
    ax[1,1].scatter(range(len(traj)), forces, c ='k')
    ax[1,0].scatter(Step, traj[Step].get_potential_energy(), c ='r')
    ax[1,1].scatter(Step, fmax(traj[Step]), c ='r')
    fig.canvas.draw_idle()

relaxatiostep.on_changed(update)


divider = make_axes_locatable(ax[0,0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(l, cax=cax)
plt.tight_layout()
#ax.set_title("Test")


plt.show()











# Iteration over series:
#for i in range(len(dataframes)):
    #print(hes_series.loc[str(i+1)])


#df = hes_series.loc['1']
#fig = px.density_heatmap(df, x=hes_series.keys(), y=hes_series.keys())
#fig.show()



#fig, ax = plt.subplots()
#im = ax.imshow(df.values)
## create an axes on the right side of ax. The width of cax will be 5%
## of ax and the padding between cax and ax will be fixed at 0.05 inch.
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#plt.colorbar(im, cax=cax)
##ax.set_title("Test")
#fig.tight_layout()
#plt.show()


























