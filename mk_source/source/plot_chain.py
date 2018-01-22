import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from optparse import OptionParser

def parse_to_list(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

if __name__=="__main__":
    parser=OptionParser()
    parser.add_option('-i','--input', default='chain_1000_1234.txt',type='string', help='input chain file (chain_1000_1234.txt)')
    parser.add_option('-p',default=None,type='string',dest='p',help='x axis parameter', action='callback', callback=parse_to_list)
    parser.add_option('--3D', default=False, action='store_true', dest='threeD',help='plot a 3D scatter plot')
    parser.add_option('-N',default=0,type='int',dest='N',help='number of chain samples to plot')
    (opts,args)=parser.parse_args()

    x = np.genfromtxt(opts.input, names=True, skip_footer=1)

    N = opts.N
    p = opts.p
    ThreeD = False
    try:
        ThreeD = opts.threeD
    except:
        pass
    Np = len(p)
    if Np > 4:
        print "Five dimensional plotting unsupported"
        exit()
    if Np == 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        S = ax.scatter(range(len(x[p[0]][-N:])),x[p[0]][-N:],c=x['logL'][-N:],s=8)
        C = plt.colorbar(S)
        C.set_label('logL')
        ax.set_xlabel(p[0])
        plt.show()
    if Np == 2:
        if ThreeD:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            S = ax.scatter(x[p[0]][-N:],x[p[1]][-N:],x['logL'][-N:],c=x['logL'][-N:],s=8)
            ax.set_xlabel(p[0])
            ax.set_ylabel(p[1])
            ax.set_zlabel('logL')
            C = plt.colorbar(S)
            C.set_label('logL')
            plt.show()
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            S = ax.scatter(x[p[0]][-N:],x[p[1]][-N:],c=x['logL'][-N:],s=8)
            ax.set_xlabel(p[0])
            ax.set_ylabel(p[1])
            C = plt.colorbar(S)
            C.set_label('logL')
            plt.show()

    if Np == 3:
        if ThreeD:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            S = ax.scatter(x[p[0]][-N:],x[p[1]][-N:],x[p[2]][-N:],c=x['logL'][-N:],s=8)
            ax.set_xlabel(p[0])
            ax.set_ylabel(p[1])
            ax.set_zlabel(p[2])
            C = plt.colorbar(S)
            C.set_label('logL')
            plt.show()
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            S = ax.scatter(x[p[0]][-N:],x[p[1]][-N:],c=x[p[2]][-N:],s=8)
            ax.set_xlabel(p[0])
            ax.set_ylabel(p[1])
            C = plt.colorbar(S)
            C.set_label(p[2])
            plt.show()
    if Np == 4:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        S = ax.scatter(x[p[0]][-N:],x[p[1]][-N:],x[p[2]][-N:],c=x[p[3]][-N:],s=8)
        ax.set_xlabel(p[0])
        ax.set_ylabel(p[1])
        ax.set_zlabel(p[2])
        C = plt.colorbar(S)
        C.set_label(p[3])
        plt.show()
