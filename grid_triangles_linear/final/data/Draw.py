from mpl_toolkits.axisartist.axislines import SubplotZero
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import matplotlib.tri as mtri
import matplotlib.patches
import csv


X, Y = [], []
fig, ax = plt.subplots()
plt.grid(True)
arrowprops={ 'arrowstyle': '-','ls':'--'}
try:
    file =  open('elem.txt', 'r')
    next(file)
    
    triangles = [line.split() for line in file]
    
    file =  open('node.txt', 'r')
    next(file) 
    
    coords = [line.split() for line in file]
    for  (c1, c2) in coords:
        X.append(float(c1))
        Y.append(float(c2))
        
    triang = mtri.Triangulation(X, Y, triangles)
    ax.set_aspect('equal')
    ax.use_sticky_edges = False
    ax.margins(0.07)
    ax.triplot(triang, color='0.8')
    ax.grid(False)
 
    
except IOError:
    print("Cant open")

try:
    file =  open('node.txt', 'r')
    next(file)
    i=0
    
    for line in file:
        values = [float(s) for s in line.split()]
        
        x=values[0]
        y=values[1]
        
        plt.annotate('%.2f'%x, xy=(x,y), xytext=(x, -0.01), 
             textcoords=plt.gca().get_xaxis_transform(),
             arrowprops=arrowprops,
             va='top', ha='center', rotation = 70, size = 10)
        plt.annotate(i, xy=(x,y), xytext=(x+0.01, y+0.05),size=3)
        plt.annotate('%.2f' % y, xy=(x,y), xytext=(-0.01, y), 
             textcoords=plt.gca().get_yaxis_transform(),
             arrowprops=arrowprops,
             va='center', ha='right', size = 10)
        plt.scatter(x,y, s=10)
        i=i+1
    
except IOError:
    print("Cant open")

ax.set_xlabel("r",fontsize = 16, loc = "right")
ax.set_ylabel("z", fontsize = 16, loc = "top", rotation='horizontal')



ax.tick_params(
    axis='both',          # changes apply to the both axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)

plt.yticks(color='w')

plt.show()

