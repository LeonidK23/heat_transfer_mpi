import os
import argparse

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

parser = argparse.ArgumentParser(description='Process parameters of the grid.')
parser.add_argument('--n', metavar='N', default=216,
                    help='dimension of the grid')

args = parser.parse_args()
N = int(args.n)

for file in os.listdir("data/"):
    if file.endswith(".txt"):
        with open("data/"+file) as f:
            data = f.readlines()

        data = [line.strip().split(" ") for line in data]
        data = np.array(data).astype(float)

        X = data[:, 0]
        Y = data[:, 1]
        Z = data[:, 2]

        X = X.reshape(X.shape[0]//N, N)
        Y = Y.reshape(Y.shape[0]//N, N)
        Z = Z.reshape(Z.shape[0]//N, N)

        # fig = plt.figure()
        # ax = Axes3D(fig)
        # ax.plot_surface(X, Y, Z, cmap=cm.coolwarm)
        # ax.set_title("Heat transfer")
        # ax.set_zlabel("Temperature", rotation=2)

        # plt.imshow(Z)

        # plt.show()
        plt.imsave('images/'+file.split('.')[0]+'.png', Z)
