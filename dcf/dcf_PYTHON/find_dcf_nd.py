from scipy.spatial import Voronoi, ConvexHull
import numpy as np


def voronoi_volumes(points):
    v = Voronoi(points)
    vol = np.zeros(v.npoints)
    for i, reg_num in enumerate(v.point_region):
        indices = v.regions[reg_num]
        if -1 in indices: # some regions can be opened
            vol[i] = np.inf
        else:
            vol[i] = ConvexHull(v.vertices[indices]).volume
    return vol


if __name__ == '__main__':
    # Test 2D
    points = [[1,2],[2,3],[3,4],[6,7],[0,2]]
    v = Voronoi(points)
    print("End")


    # Test 3D