from scipy.spatial import Voronoi
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
# Generate simple density correction function with voronoi areas from k-space trajectory
# Gehua Tong, Oct 2019


# TODO test this code against MATLAB

def find_dcf_2d(ktrajs, opt_save=True, savename='dcf.npy'):
    """
    Parameters
    ----------
    ktrajs : np.ndarray
        (n x 2) or (m x n x 2) array of 2D k-space trajectory. Last dimension is x vs. y
    opt_save : bool, optional
        Whether to save dcf areas into file. Default is True.
    savename : str, optional
        Name for saved dcf file. Default is 'dcf.npy'.

    Returns
    -------
    dcf : np.ndarray
        Length-n or (m x n) array of density correction values
        Equal to voronoi areas in voronoi diagram calculated from k-space points

    """

    # Load k-space trajectory
    sizek = np.shape(ktrajs)
    if sizek[-1] != 2:
        raise ValueError("Last dimension of k-trajs should always be 2")

    # Convert into complex array
    kk = ktrajs
    if len(sizek) == 3:
        kk = np.reshape(ktrajs, (sizek[0] * sizek[1], sizek[2]))

    kxy = kk[:,0] + 1j * kk[:,1]

    # Scale k-space points to prevent "non-unique points excluded" operation
    dcf_out = voronoi_area(kxy)
    print(dcf_out)
    #Normalize and reshape
    dcf_out = dcf_out / np.max(dcf_out)
    if len(sizek) == 3:
        dcf_out = np.reshape(dcf_out, (sizek[0],sizek[1]))


    if opt_save:
        np.save(savename, 'dcf_out')

    return dcf_out


def voronoi_area(kxy):
    """ Calculate the area(density) of a k-space trajectory via Voronoi triangulation

    Parameters
    ----------
    k : np.ndarray
        k-space trajectory ([-0.5 0.5]*res)

    Returns
    -------
    v_areas : np.ndarray
        Voronoi areas associated with the k-space points

    Adopted from Florian Wiesinger's 07/2006 MATLAB code
    Gehua Tong 10/2019
    """

    fufa = 0.7 # Fudge factor to reduce outside areas (why?)

    kxy = np.array(kxy)
    ku, ki  = np.unique(kxy, return_index = True)
    nupts = len(kxy) - len(ku) # Number of non-unique points (?)

    nsimi = np.sum(np.absolute(np.diff(kxy[np.argsort(np.absolute(kxy))])) < 0.02)


    if nsimi > nupts:
        print('voronoi_area:nsimi', 'Detected {} similar points -> check dcf'.format(nsimi))


    if nupts > 0:
        print('{0:d} non-unique points excluded'.format(nupts))


    kx = np.real(ku)
    ky = np.imag(ku)
    kxy = np.transpose([kx, ky])

    # Create Voronoi object

    vor = Voronoi(kxy)

    dcf = np.zeros(len(ku))

    makx = np.max(np.absolute(kx))
    maky = np.max(np.absolute(ky))
    makxy = np.min([makx, maky])

    #outside_cnt = 0


    all_vertices = vor.vertices
    regions = vor.regions

    # Loop through all polygons
    for reg_ind in range(np.shape(regions)[0]):
        x = all_vertices[regions[reg_ind],0]
        y = all_vertices[regions[reg_ind],1]

        # Check for points outside the sampled area
        if any(np.isinf(x)) or any(np.isinf(y)) or any(np.absolute(x) > makx) \
            or any(np.absolute(y) > maky) or any(np.sqrt(x**2 + y**2) > makxy):
            # outside: radius -> approximate area
            #outside_cnt += 1 # count for number of "outside" points, it seems
            radii = np.sort(np.sqrt((x - kx[reg_ind])**2 + (y - ky[reg_ind])**2))
            dcf[reg_ind] = radii[0] * radii[1] * fufa
        # dcf(i1) = mean(radii(1:2)) ^ 2
        # dcf(i1) = mean(radii(1:2)) ^ 2 * pi
        else:
            # Inside: correct area
            dcf[reg_ind] = polygon_area(np.transpose([x,y]))


    # Set infinite areas to zero and warn
    if any(np.isinf(dcf)):
        dcf[np.isinf(dcf)] = 0
        print('Warning!: ', 'voronoi_area:inf;', 'dcf contains inf -> setting to zero')


    # Only set points at indices returned by np.unique() to the correct areas; others remain zero, it seems
    dcf_out = np.zeros(len(kxy))
    dcf_out[ki] = dcf
    mean_dcf_unsamp = np.mean(dcf_out[dcf_out > 1])

    if mean_dcf_unsamp > 1:
        print('Warning!: ', 'voronoi_area:unsamp; ', 'Numerically stable?')

        print("mean(dcf(dcf>1))={0:.2f}".format(mean_dcf_unsamp))

        plt.figure()
        plt.plot(dcf_out)

    return dcf_out

def polygon_area(vertices):
    # Calculates area of a 2D polygon from its vertices
    """
    Parameters
    ----------
    vertices : array_like
         (n x 2) array of (x,y) coordinates of the 2D polygon's vertices

    Returns
    -------
    area : float
        Numerical area of polygon

    """



    if np.shape(vertices)[1] != 2:
        raise ValueError("Second dimension must be 2!")

    # TODO Perform checks



    # Use set formula to calculate areas
    x = vertices[:,0]
    y = vertices[:,1]
    area = np.absolute(np.sum((np.roll(x,1) - x) * (np.roll(y,1) + y))/2)

    return area


if __name__ == '__main__':
    ktrajs = np.transpose(np.load('./data/Sep27_2019/rad_2d_fid_interleave2_N256_Ns5_ktrajs.npy'))
    dcf = find_dcf_2d(ktrajs, opt_save=True, savename='./data/Sep27_2019/rad_2d_fid_interleave2_N256_Ns5_dcf.npy')
    dcf0 = sio.loadmat('./data/Sep27_2019/rad_2d_fid_interleave2_N256_Ns5_dcf.mat')['dcf_out']
    # TODO dcf and dcf0 are too different! Fix this
    plt.show()