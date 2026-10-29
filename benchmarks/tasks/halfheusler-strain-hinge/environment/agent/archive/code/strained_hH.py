"""
----------------------------------------------------
Hinge electronic structure of strained half-Heuslers
----------------------------------------------------

In this module we reproduce some of the results presented in the paper:

S. K. Das, I. C. Fulga, R. Dhawan, H. C. Kandpal, J. van den Brink,
and  J. I. Facio
"Hinge electronic structure of strained half-Heuslers"
arXiv:2504.20190.

For examples of usage, see the main() function, which reproduces
some of our numerical results. This script can be imported in a
python interface or simply run as:

python3 strained_hH.py
"""
import numpy as np
import pickle
import kwant
import tinyarray as ta
import sys

import pylab as py
py.ion()

p = dict(k_x=0, k_y=0, k_z=0)


def update_progress(progress, decimalpoints=0):
    """ Make an interactive progress bar as described on:
    https://stackoverflow.com/questions/3160699/python-progress-bar
    """
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + 
            "-"*(barLength-block), round(progress*100, decimalpoints), status)
    sys.stdout.write(text)
    sys.stdout.flush()


def process_hamdata(infile, outfile, keep_only_atom_name=True):
    """ Read a "+hamdata" file containing information on a Wannier Hamiltonian,
    as obtained using density-functional theory methods. Using this data, save
    the relevant information about the Wannier Hamiltonian to a binary "ham"
    file. This file can then be further used in tight-binding model
    construction.

    Parameters
    ----------
    infile : string
        Name of the "+hamdata" file.
    outfile : string
        Name of the output "ham" file.
    keep_only_atom_name : boolean
        If True, remove the spin component and orbital information from the
        'wannames' stored in the +hamdata file.

    Notes
    -----
    This function assumes that the +hamdata file contains the following
    information: 'lattice_vectors', 'nwan', 'wannames', 'wancenters', as well
    as all of the Hamiltonian matrix elements, marked using 'Tij, Hij:'
    and 'end Tij, Hij:'.
    """
    with open(infile, 'r+') as f:
        lines = f.readlines()
        i = 0
        while not lines[i].startswith('lattice_vectors'):
            i += 1

        # extract 3 lattice vectors
        a = np.zeros((3, 3), dtype=float)
        for j in range(3):
            a[j,:] = [float(s) for s in lines[i+j+1].split()]

        inva = np.linalg.inv(a.T) # used to determine 'integer_distance_ij'

        while not lines[i].startswith('nwan'):
            i += 1

        nwan = int(lines[i+1]) # total number of Wannier orbitals

        while not lines[i].startswith('wannames'):
            i += 1

        wannames = [] # names of Wannier orbitals
        for ind in range(nwan):
            if lines[i+1+ind][-1:] == '\n':
                wannames.append(lines[i+1+ind][:-1])
            else:
                wannames.append(lines[i+1+ind])

        while not lines[i].startswith('wancenters'):
            i += 1

        wcen = np.loadtxt(lines[i+1:i+1+nwan]) # positions of all orbitals

        # remove the spin component and orbital name from Wannier names
        if keep_only_atom_name:
            for ind in range(nwan):
                wannames[ind] = wannames[ind][:wannames[ind].find(' ')]

        atom_pos = [] # atom positions
        atom_names = [] # atom names
        atom_orbs = [] # nr. orbitals for each atom
        orb_ind = [] # to which atom does each orbital belong
        orb_ind_start = [] # starting orbital index of each atom
        orb_ind_stop = [] # final orbital index of each atom

        orb_ind_start.append(0) # first atom starts with orbital nr. 0

        atom_pos.append(wcen[0, :])
        atom_names.append(wannames[0])
        orb_current = wcen[0, :]
        orb_counter = 0
        nratoms = 1
        for ind in range(nwan):
            orb_counter += 1
            if np.linalg.norm(np.abs(wcen[ind, :] - orb_current)) >= 1e-5:
                atom_pos.append(wcen[ind, :])
                atom_names.append(wannames[ind])
                orb_current = wcen[ind, :]
                nratoms +=1
                atom_orbs.append(orb_counter-1)
                orb_counter = 1
                orb_ind_start.append(ind)
                orb_ind_stop.append(ind)

            orb_ind.append(nratoms - 1)

        # last atom
        atom_orbs.append(orb_counter)
        orb_ind_stop.append(nwan)

        # determine which parts of the file contain the Wannier data
        data_start = []
        data_stop = []
        for ind in range(len(lines)):
            if lines[ind].startswith('Tij, Hij:'):
                data_start.append(ind)

            if lines[ind].startswith('end Tij, Hij:'):
                data_stop.append(ind)

        if len(data_start) != nwan**2:
            raise ValueError('The number of entries in the data file' +
                             ' does not mach the number of' +
                             ' Wannier orbitals.')

        # begin extracting Hamiltonain elements (hamiltonian_submatrix)
        total_hamdata = {}

        datalen = len(data_start)

        for ind in range(datalen):
            update_progress((ind+1)/datalen)
            row, col = [int(s)-1 for s in lines[data_start[ind]+1].split()]
            tempdata = np.loadtxt(
                        lines[data_start[ind]+2:data_stop[ind]])

            for ind2 in range(tempdata.shape[0]):
                nj = np.dot(inva, tempdata[ind2, :3].reshape((3, 1)) -
                                  wcen[col, :].reshape((3, 1)) + 
                                  wcen[row, :].reshape((3, 1))).T

                # sanity check: integer_distance_ij is integer
                for jj in nj[0]:
                    if np.abs(jj - np.round(jj)) > 1e-4:
                        print(jj, np.round(jj))
                        raise ValueError('Non-integer neighbor hopping...')

                nj = [int(j) for j in np.round(nj[0])]

                which_hop = str(orb_ind[row]) + ' ' + \
                            str(orb_ind[col]) + ' ' + \
                            str(nj)

                bigtmat = np.zeros((nwan, nwan), dtype=complex)
                bigtmat[row, col] += tempdata[ind2, 3] + 1j * tempdata[ind2, 4]

                tmat = bigtmat[orb_ind_start[orb_ind[row]]:
                                    orb_ind_stop[orb_ind[row]],
                               orb_ind_start[orb_ind[col]]:
                                    orb_ind_stop[orb_ind[col]]]

                if which_hop in total_hamdata:
                    total_hamdata[which_hop] += tmat
                else:
                    total_hamdata[which_hop] = tmat

    f = open(outfile, 'wb')
    pickle.dump([a, atom_pos, atom_names, atom_orbs, total_hamdata], f)
    f.close()


def build_template_from_ham_file(myfile):
    """ Load a Hamiltonian file as generated by process_hamdata() and return
    a kwant template.

    Parameters
    ----------
    myfile : string
        Name of the "ham" file.

    Returns
    -------
    template : kwant.builder.Builder
        Kwant template.
    """
    f = open(myfile, 'rb')
    a, atom_pos, atom_names, atom_orbs, di = pickle.load(f)
    f.close()

    lat = kwant.lattice.general(a, np.array(atom_pos), name=atom_names,
                                norbs=atom_orbs)

    subl = lat.sublattices

    ditems = []

    for key, val in di.items():
        ditems.append([key, val])

    template = kwant.Builder(kwant.lattice.TranslationalSymmetry(*(a)))

    onsite_terms = []
    hopping_terms = []

    for ind in range(len(ditems)):
        mystring = ditems[ind][0]

        # index of first sublattice
        s1 = int(mystring[:mystring.find(' ')])
        mystring = mystring[mystring.find(' ')+1:]

        # index of second sublattice
        s2 = int(mystring[:mystring.find(' ')])
        mystring = mystring[mystring.find(' ')+1:]

        # vector of integers specifying which unit cell to hop to
        vec = eval(mystring)

        if np.linalg.norm(vec) == 0 and s1 == s2: # onsite
            onsite_terms.append([s1, vec, ditems[ind][1]])

        else: # hopping
            hopping_terms.append([s1, s2, vec, ditems[ind][1]])

    for ind in range(len(onsite_terms)):
        s1, vec, ons = onsite_terms[ind]
        template[subl[s1](*vec)] = ons

    for ind in range(len(hopping_terms)):
        s1, s2, vec, hop = hopping_terms[ind]
        template[kwant.builder.HoppingKind(vec, subl[s1], subl[s2])] = hop

    return template


def build_bulk_H(template):
    """ Build the bulk Hamiltonian.

    Parameters
    ----------
    template : kwant.builder.Builder
        Hamiltonian template, as returned by build_template_from_ham_file().

    Returns
    -------
    fsys : kwant.builder.FiniteSystem
        Kwant bulk system.
    """
    fsys = kwant.wraparound.wraparound(template).finalized()
    return fsys


def shape_main_text(site):
    """ System shape used in Fig. 5 of the main text.

    Parameters
    ----------
    site : kwant.builder.Site
        A site of the system.

    Returns
    -------
    in_shape : boolean
        Returns True if the site is inside the shape, and False otherwise.
    """
    shapesize = 20
   
    return -0 <= (site.tag[0]) <= shapesize and \
           -0 <= (site.tag[1]) <= shapesize and \
           -0 <= (site.tag[2]) <= shapesize


def shape_appendix(site):
    """ System shape used in Appendix B of the paper.

    Parameters
    ----------
    site : kwant.builder.Site
        A site of the system.

    Returns
    -------
    in_shape : boolean
        Returns True if the site is inside the shape, and False otherwise.
    """
    mysizeL = -14
    mysizeU = -mysizeL + 1

    if (site.tag[1] + site.tag[2] == mysizeU) and \
       (site.family.name=='Zn2'):
        return False

    if (site.tag[1] + site.tag[2] == mysizeL) and \
       (site.family.name in ['Sb3', 'Bi3']):
        return False

    if (site.tag[1] == mysizeU) or (site.tag[1] == mysizeL):
        return False

    return -4 <= (site.tag[0]) <= 4 and \
           mysizeL <= (site.tag[1] - site.tag[2]) <= mysizeU and \
           mysizeL <= (site.tag[1] + site.tag[2]) <= mysizeU


def build_prism(template, shape=shape_main_text,
                inf_dir=(1, 0, 0), coord_name='x'):
    """ Construct the system in a prism geometry.

    Parameters
    ----------
    template : kwant.builder.Builder
        Kwant template.
    shape : function
        Function that specifies the shape of the unit cell. Two examples
        included in this file are shape_main_text() and shape_appendix().
    inf_dir : tuple of integers
        Translationally-invariant direction of the prism, expressed as an
        an integer multiple of the three Bravais vectors.
    coord_name : string
        Label given to the conserved momentum along the prism. For example,
        coord_name = 'x' will lead to a momentum variable labeled 'k_x'.

    Returns
    -------
    sys : kwant.builder.FiniteSystem
        Kwant system in a prism geometry.
    """

    periods = np.array(template.symmetry.periods)
    trsym = np.zeros(periods.shape[0])

    for ind, mydir in enumerate(inf_dir):
        trsym += mydir * periods[ind]

    sym = kwant.lattice.TranslationalSymmetry(trsym)
    sys = kwant.Builder(sym)

    sys.fill(template, shape, (0.,0.,0.))
    sys = kwant.wraparound.wraparound(sys, 
                    coordinate_names=coord_name).finalized()

    return sys


def get_projected_positions(fsys, normalvec=np.array([0., 1., 1.]), 
                 planarvecs=(np.array([1., 0., 0.]), np.array([0., 1., -1.]))):
    """ Return the coordinates of the lattice sites of fsys, projected onto a
    two-dimensional plane. This is useful for plotting wavefunction probability
    densities.

    Parameters
    ----------
    fsys : kwant.builder.FiniteSystem
        Kwant system, as returned by build_prism().
    normalvec : 1d array
        Coordinates of the 3d vector perpendicular to the projection plane.
    planarvecs : tuple of two 1d arrays
        Coordinates of the two 3d vector parameterizing the projection plane.

    Returns
    -------
    coord1, coord2 : lists of floats
        Coordinates of the sites' projections on the 2d plane.
    """
    coord1 = []
    coord2 = []

    normalvec /= np.linalg.norm(normalvec)

    pv1, pv2 = planarvecs

    pv1 /= np.linalg.norm(pv1)
    pv2 /= np.linalg.norm(pv2)

    for ind, mysite in enumerate(fsys.sites):
        vec = mysite.pos
        dist = np.sum(normalvec * vec)
        projpoint = vec - dist * normalvec
        coord1.append(np.sum(pv1 * projpoint))
        coord2.append(np.sum(pv2 * projpoint))

    return coord1, coord2


class Amorphous(kwant.builder.SiteFamily):
    """ This class enables to create kwant systems with arbitrary arrangements
    of lattice sites. This is useful for plotting purposes.
    """
    def __init__(self, coords, norbs=None):
        self.coords = coords
        self.norbs = norbs
        super(Amorphous, self).__init__("amorphous", "", norbs)

    def normalize_tag(self, tag):
        try:
            tag = int(tag[0])
        except:
            raise KeyError

        if 0 <= tag < len(self.coords):
            return tag
        else:
            raise KeyError

    def pos(self, tag):
        return ta.array(self.coords[tag])


def plotter_map_proj(fsys, psi, 
                 normalvec=np.array([0., 1., 1.]), 
                 planarvecs=(np.array([1., 0., 0.]), np.array([0., 1., -1.])),
                 bkg_color='#FFFFFF', oversampling=12):
    """ Plot the real-space probability density associated to a state.

    Parameters
    ----------
    fsys : kwant.builder.FiniteSystem
        Kwant system, as returned by build_prism().
    psi : 1d array
        Probability density of the eigenvector.
    normalvec : 1d array
        Coordinates of the 3d vector perpendicular to the plane of the system
        unit cell.
    planarvecs : tuple of two 1d arrays
        Coordinates of the two 3d vector parameterizing the unit cell of the
        system.
    bkg_color : matplotlib color spec
        Areas without sites are filled with this color.
    oversampling : integer
        Number of pixels per reference length. This makes the plot look
        smoother.
    """
    # make dummy kwant system for map plot
    c1, c2 = get_projected_positions(fsys, normalvec, planarvecs)

    coords = []
    for ind in range(len(c1)):
        coords.append([c1[ind], c2[ind]])

    temp_lat = Amorphous(coords, norbs=1)
    temp_sys = kwant.Builder()

    for i in range(len(coords)):
        temp_sys[temp_lat(i)] = 0.0

    # reduce vector to number of sites
    newpsi = np.zeros(len(c1))

    ind = 0
    for ind2, mysite in enumerate(fsys.sites):
        if mysite.family.name == 'Zn2':
            newpsi[ind2] = np.sum(psi[ind:ind+2])
            ind += 2
        if mysite.family.name in ['Sb3', 'Bi3']:
            newpsi[ind2] = np.sum(psi[ind:ind+6])
            ind += 6

    temp_sys = temp_sys.finalized()
    kwant.plotter.map(temp_sys, newpsi,
                      background=bkg_color, oversampling=oversampling)


def plot_project(fsys, normalvec=np.array([0., 1., 1.]), 
                 planarvecs=(np.array([1., 0., 0.]), np.array([0., 1., -1.])), 
                 colors=dict(Sb3='tab:olive', Bi3='tab:olive',
                                              Zn2='royalblue'),
                 sizes=20, draw_lines=True, line_width=1.0, line_range=6.0):
    """ Plot the system, projected onto the plane of the unit cell of the
    prism.

    Parameters
    ----------
    fsys : kwant.builder.FiniteSystem
        Kwant system, as returned by build_prism().
    normalvec : 1d array
        Coordinates of the 3d vector perpendicular to the plane of the system
        unit cell.
    planarvecs : tuple of two 1d arrays
        Coordinates of the two 3d vector parameterizing the unit cell of the
        system.
    colors : dictionary
        Specifies which colors to use for the different atoms. The keys of the
        dictionary should be the Site.family.name, and the value should be a
        matplotlib color descriptor.
    sizes : float or 1d array of floats
        Specifies the size of each of the plotted sites. If a single float is
        used, the same size will be used for all sites.
    draw_lines : boolean
        Connect nearby sites with lines to better visualize the system.
    line_width : float
        Width of the lines connecting nearby sites.
    line_range : float
        If the distance between the projected positions of two sites is smaller
        than line_range, the two sites will be connected with a hopping.
    """
    py.figure()

    if not hasattr(sizes, 'len'):
        sizes = [sizes,] * len(fsys.sites)

    c1, c2 = get_projected_positions(fsys, normalvec, planarvecs)

    if draw_lines:
        for ind1 in range(len(c1)):
            pos1 = np.array([c1[ind1], c2[ind1]])
            for ind2 in range(ind1+1, len(c1)):
                pos2 = np.array([c1[ind2], c2[ind2]])
                if np.linalg.norm(pos1-pos2) <= line_range:
                    py.plot([c1[ind1], c1[ind2]], 
                            [c2[ind1], c2[ind2]], c='k', lw=line_width,
                            zorder=1)

    for ind, mysite in enumerate(fsys.sites):
        py.scatter([c1[ind],], [c2[ind],], 
                    c=colors[mysite.family.name], s=sizes[ind],
                    zorder=2)


def pos_H(fsys, p, coord=0):
    """ Calculate the position operator in the 'coord' direction of the
    Hamiltonian of fsys.
    """

    H, ton, fon = fsys.hamiltonian_submatrix(return_norb=True,
                                            params=p, sparse=True)
    x = np.zeros(H.shape[0])
    ind = 0
    for i in range(len(fsys.sites)):
        for j in range(ind, ind + ton[i]):
            x[j] = fsys.sites[i].pos[coord]

        ind += ton[i]

    return x


def tag_H(fsys, p, coord=0):
    """ Return a diagonal matrix containing the integer index characterizing
    the unit cells in the 'coord' direction of the Hamiltonian of fsys.
    """

    H, ton, fon = fsys.hamiltonian_submatrix(return_norb=True,
                                            params=p, sparse=True)
    x = np.zeros(H.shape[0])
    ind = 0
    for i in range(len(fsys.sites)):
        for j in range(ind, ind + ton[i]):
            x[j] = fsys.sites[i].tag[coord]

        ind += ton[i]

    return x


def translate_BZ_points(kstart):
    """ Translate commonly used names for high-symmetry momenta, such as
    'G', 'X', or 'M', into a list of floats containing the values of those
    momental in the 3D Brillouin zone.
    """
    if kstart == 'G':
        return [0, 0, 0]
    elif kstart == 'X':
        return [np.pi, 0, 0]
    elif kstart == '-X':
        return [-np.pi, 0, 0]
    elif kstart == 'Y':
        return [0, np.pi, 0]
    elif kstart == '-Y':
        return [0, -np.pi, 0]
    elif kstart == 'Z':
        return [0, 0, np.pi]
    elif kstart == '-Z':
        return [0, 0, -np.pi]
    elif kstart == 'Z':
        return [0, 0, np.pi]
    elif kstart == 'P':
        return [np.pi, np.pi, np.pi]
    elif kstart == '-P':
        return [-np.pi, -np.pi, -np.pi]
    elif kstart == 'M':
        return [np.pi, np.pi, 0]
    elif kstart == '-M':
        return [-np.pi, -np.pi, 0]
    elif kstart == 'U':
        return [np.pi, 0, np.pi]
    elif kstart == '-U':
        return [-np.pi, 0, -np.pi]

    else:
        return kstart


def plot_bands(fsys, kstart='G', kstop='X', npoints=101, Ef=0,
               ylims=[-0.5, 1.25], color_choice='shape_main_text'):
    """ Plot the bandstructure of the system in the prism geometry.

    Parameters
    ----------
    fsys : kwant.builder.FiniteSystem
        Kwant system, as returned by build_prism().
    kstart, kstop : string or 1d array
        Coordinates of the ends of the momentum interval for which to plot the
        bandstructure. These can be specified either as 1d arrays of three
        floats, e.g. [0.0, 0.0, 0.0], or as strings denoting high-symmetry
        momenta, e.g. 'G', 'X', or 'M'.
    npoints : integer
        Number of momenta to include in the plotted interval.
    Ef : float
        Shift the Fermi level by this ammount.
    ylims : list of floats
        Range of energies included in the plot.
    color_choice : 'shape_main_text' or 'shape_appendix'
        Specify which geometry was used to construct the prism. This is useful
        for determining which sites are labeled as hinge, surface, and bulk.
    """

    k_start = np.array(translate_BZ_points(kstart))
    k_stop = np.array(translate_BZ_points(kstop))
    evals = []
    kvals = []
    colors = []
    sizes = []

    if color_choice == 'shape_main_text':
        tag1 = tag_H(fsys, p, 1)
        tag2 = tag_H(fsys, p, 2)
        maxtag1 = np.max(tag1)
        maxtag2 = np.max(tag2)
        hinge_inds1 = np.intersect1d(np.where(tag1 == maxtag1)[0], 
                                     np.where(tag2 == maxtag2)[0])

        hinge_inds2 = np.intersect1d(np.where(tag1 == 0)[0], 
                                     np.where(tag2 == maxtag2)[0])

        hinge_inds3 = np.intersect1d(np.where(tag1 == maxtag1)[0], 
                                     np.where(tag2 == 0)[0])

        hinge_inds4 = np.intersect1d(np.where(tag1 == 0)[0], 
                                     np.where(tag2 == 0)[0])

        hinge_inds = np.concatenate([hinge_inds1, hinge_inds2, 
                                     hinge_inds3, hinge_inds4])

        surf_inds1 = np.where(tag1 == maxtag1)[0]
        surf_inds2 = np.where(tag1 == 0)[0]
        surf_inds3 = np.where(tag2 == maxtag2)[0]
        surf_inds4 = np.where(tag2 == 0)[0]

        surf_inds = np.union1d(np.concatenate([surf_inds1, surf_inds2]), 
                               np.concatenate([surf_inds3, surf_inds4]))

        surf_inds = np.setdiff1d(surf_inds, hinge_inds)

        bulk_inds = np.setdiff1d(range(len(tag1)),
                                 np.union1d(surf_inds, hinge_inds))

    if color_choice == 'shape_appendix':
        c1 = pos_H(fsys, p, 0)
        posy = pos_H(fsys, p, 1)
        posz = pos_H(fsys, p, 2)

        c2 = (posy - posz) / np.sqrt(2)

        # only works for the function shape_appendix() as currently written
        L_edge = -76
        R_edge = 85
        T_edge = 52
        B_edge = -57

        L_surf_inds = np.where(c1 < L_edge)[0]
        R_surf_inds = np.where(c1 > R_edge)[0]

        T_surf_inds = np.where(c2 > T_edge)[0]
        B_surf_inds = np.where(c2 < B_edge)[0]

        bulk_inds = np.setdiff1d(range(len(c1)), L_surf_inds)
        bulk_inds = np.setdiff1d(bulk_inds, R_surf_inds)
        bulk_inds = np.setdiff1d(bulk_inds, T_surf_inds)
        bulk_inds = np.setdiff1d(bulk_inds, B_surf_inds)

        surf_inds = np.union1d(np.concatenate([L_surf_inds, R_surf_inds]), 
                               np.concatenate([T_surf_inds, B_surf_inds]))

        hinge_inds = np.concatenate([
                        np.intersect1d(L_surf_inds, T_surf_inds),
                        np.intersect1d(L_surf_inds, B_surf_inds),
                        np.intersect1d(R_surf_inds, T_surf_inds),
                        np.intersect1d(R_surf_inds, B_surf_inds),
                        ])

        surf_inds = np.setdiff1d(surf_inds, hinge_inds)

    for ind in np.linspace(0, 1, npoints):
        update_progress(ind)
        kval = k_start + ind * (k_stop - k_start)
        p['k_x'] = kval[0]; p['k_y'] = kval[1]; p['k_z'] = kval[2]
        Ham = fsys.hamiltonian_submatrix(params=p)
        Ham += Ef * np.eye(Ham.shape[0])
        bands, evecs = np.linalg.eigh(Ham)

        for ind2 in range(len(bands)):
            vec = np.abs(evecs[:, ind2])**2

            colors.append(np.sum(vec[hinge_inds]) +
                          np.sum(vec[surf_inds]))

            sizes.append(np.sum(vec[hinge_inds]))
            kvals.append(ind)
            evals.append(bands[ind2])

    py.figure()
    py.title('kstart='+repr(kstart)+'. kstop='+repr(kstop))
    
    sizes = np.array(sizes)
    colors = np.array(colors)
    kvals = np.array(kvals)
    evals = np.array(evals)

    inds = (np.abs(colors).argsort())

    sizes = sizes[inds]
    colors = colors[inds]
    kvals = kvals[inds]
    evals = evals[inds]

    sizes = 1 + sizes / np.max(sizes) * (100 - 1) # rescale sizes in [1, 100]

    py.scatter(kvals, evals, c=colors, linewidth=0, vmin=0, vmax=1,
               cmap='viridis', s=sizes)
    py.colorbar()
    py.ylim(ylims)
    py.xlim([0, 1])


def main():
    """ Here we reproduce some of our numerical results. Given the system
    sizes, some of the calculations may take a long time to complete. For
    reference, on an Intel(R) Xeon(R) W-2145 machine (8 cores, 16 threads,
    produced in 2017), running the main() function took 25 minutes.

    The code can be stopped at any time by pressing Ctrl+C.
    """

    print('# Extracting data from Wannier Hamiltonian files.')
    print()
    print('# LiBiZn, no strain')
    process_hamdata('../wannier/+hamdata_libizn_eq_orig_setup.txt',
                    '../kwant_ham/LiBiZn_nostrain.ham')
    print()
    print()

    print('# LiBiZn, strained')
    process_hamdata('../wannier/+hamdata_libizn_eq_strained_setup.txt',
                    '../kwant_ham/LiBiZn_strain.ham')
    print(); print()

    print('# LiSbZn, no strain')
    process_hamdata('../wannier/+hamdata_lisbzn_eq_orig_setup.txt',
                    '../kwant_ham/LiSbZn_nostrain.ham')
    print(); print()

    print('# LiSbZn, strained')
    process_hamdata('../wannier/+hamdata_lisbzn_eq_strained_setup.txt',
                    '../kwant_ham/LiSbZn_strain.ham')
    print(); print()

    print('# Constructing systems in the prism geometry of the main text')
    template_LiBiZn_nostrain = \
            build_template_from_ham_file('../kwant_ham/LiBiZn_nostrain.ham')

    template_LiSbZn_nostrain = \
            build_template_from_ham_file('../kwant_ham/LiSbZn_nostrain.ham')

    sys_Bi = build_prism(template_LiBiZn_nostrain)
    sys_Sb = build_prism(template_LiSbZn_nostrain)

    print()

    print('# Plotting prism unit cell')

    plot_project(sys_Bi)

    print()

    print('# Plotting hinge mode of unstrained LiSbZn')

    H = sys_Sb.hamiltonian_submatrix(params=p)
    E, V = np.linalg.eigh(H)
    plotter_map_proj(sys_Sb, np.abs(V[:, 2646])**2, bkg_color='#f3f3f3ff')

    print()

    print('# Now the code will begin plotting a low-resolution version of')
    print('# the bandstructures shown in panels (b) and (c) of Fig. 5. This')
    print('# may take a long time, depending on your machine. It is always')
    print('# possible to interrupt the code by pressing Ctrl+C.')

    print()

    asdasd = input('# Press Enter to begin calculating the bandstructure'+
                   ' of Fig. 5(b).')
    print()
    plot_bands(sys_Bi)
    print(); print()


    asdasd = input('# Press Enter to begin calculating the bandstructure'+
                   ' of Fig. 5(c).')
    print()
    plot_bands(sys_Sb)
    print(); print()

    asdasd = input('# Press Enter to exit...')


if __name__ == '__main__':
    main()
