#!/usr/bin/env python
"""Searches crystal structures (cif files) for particular d-ligand-d angles.

    For each structure, we employ the following strategy:

        1. Find all sites fully occupied with a given d-element.
        2. For every d-site, find its ligands and neighboring d-atoms.
        3. For every neighbor, we check if any of the ligands surrounding
           the central atom is also a ligand of this neighboring atom.
           If it is the case, we calculate the central-ligand-neighbor angle. 
        4. If the number of angles is equal to `nligands` or more, we add
            the respective angles to the preliminary list.
        5. If any of the angles in the preliminary list falls into the given
           range, we put the respective material into the final output list.

    The output is given as a markdown and a pdf file (the latter is generated
    using pandoc with the `eisvogel <https://github.com/Wandmalfarbe/pandoc-latex-template/blob/master/eisvogel.tex>`_
    template. Structures are imported via pymatgen.

    Example:

        ./anglefinder.py -d Cu -l O --ddmax 3.1 --dlmax 2.1 -a 93.6 94.0

    In this example, the scripts finds all cuprates with Cu--O--Cu angles
    falling in the range between 93.6 and 94 degrees. Oxygen atoms are
    treated as ligands only if the Cu--O distance is smaller or equal to
    2.1 Angstroem. Cu atoms are treated as neighbors if the respective
    distance does not exceed 3.1 Angstroem.

"""

from __future__ import absolute_import, division, with_statement

__author__  = 'Oleg Janson'
__date__    = 'September 30, 2019'
__email__   = 'olegjanson@gmail.com'
__version__ = '0.11'

import os
import argparse
import numpy as np
from datetime import datetime

import pypandoc as pd

from numpy.linalg import norm

import pymatgen as mg
from pymatgen.util.plotting import format_formula

description_text = "Searches for specific d-ligand-d angles in cif ciles."
epilog_text = "Please report any bugs to Oleg Janson <olegjanson@gmail.com>."

parser = argparse.ArgumentParser(description = description_text,
                                 conflict_handler = "resolve",
                                 formatter_class =\
                                     argparse.RawTextHelpFormatter,
                                 epilog = epilog_text)

parser.add_argument("cifdir", type=str, default=".",
                    help='the directory with the cif files.')

# Required arguments
parser.add_argument('-d', type=str, required=True, help='d atom type.')
parser.add_argument('-l', type=str, required=True, help='ligand atom type.')
parser.add_argument('--ddmax', type=np.float64, required=True,
                    help='maximal d-d distance.')
parser.add_argument('--dlmax', type=np.float64, required=True,
                    help='maximal d-ligand distance.')
parser.add_argument('-a', '--angles', nargs=2, type=np.float64, required=True,
                    help=('minimal and maximal d-ligand-d angles.'))

# Optional arguments 
parser.add_argument('--dvalence', type=int, default=None,
                    help='valence of the d atom. Default: None.')
parser.add_argument('--nligands', type=int, default=2,
                    help='maximal number of ligands. Default: 2.')
parser.add_argument('--nangles', type=int, default=6,
                    help='number of angles per line (output). Default: 6.')
parser.add_argument('--prefix', type=str,
                    default=datetime.now().strftime("%Y%m%d_%H%M%S%f"),
                    help='file names for MD and PDF output files.')

# Set the parameters
args = parser.parse_args()
d_atom_type = args.d
ligand_type = args.l
d_valence = args.dvalence
max_d_to_ligand_distance = args.dlmax
max_d_to_d_distance = args.ddmax
min_number_of_common_ligands = args.nligands
min_angle, max_angle = args.angles[0], args.angles[1]
max_angles_on_single_line = args.nangles
cifdir = args.cifdir
prefix = args.prefix

def get_angle(v1, v2):
    """Calculates the angle (in deg.) between two vectors.

    Args:
        v1 (np.array): The first vector.
        v2 (np.array): The second vector.

    Returns:
        (np.float): Angle in degrees.

    """
    return (180 / np.pi) * \
           np.arccos(np.dot(v1, v2) / norm(v1) / norm(v2))

def check_if_site_is_ordered(structure, site_index):
    """Checks if a given site in a given structure is fully ordered.

    Args:
        structure (pymatgen.Structure): A pymatgen structure.
        site_index (int): The site index.

    Returns:
        (bool)

    """
    return (len(structure[site_index].species.elements) == 1)

def get_environment(structure, site_index, max_distance,
                    symbol=None, oxi_state=None):
    """Returns the enviromenent of a given site.

    All incompletely ordered sites are discarded.

    Args:
        structure (pymatgen.Structure): A pymatgen structure.
        site_index (int): Index of the central atom.
        max_distance (float): Maximal distance to the cenral atom
            in Angstroem.
        symbol (str): Returns only atoms of this type.
        oxi_state (int): Returns only atoms with this oxidation state.

    Returns:
        (list):

    """
    environment = structure.get_neighbors(structure[site_index], max_distance)
    if symbol is not None:
        environment = [site for site in environment
                       if (site[0].species.elements[0].symbol == symbol)]
    if oxi_state is not None:
        environment = [site for site in environment
                       if (site[0].species.elements[0].oxi_state == oxi_state)]
    return [site.site.coords for site in environment
            if check_if_site_is_ordered(structure, site.index)]

# Output: title and header
out = "# anglefinder"
out += "\nversion {:s}".format(__version__)
out += "\nbased on pymatgen version {:s}\n".format(mg.__version__)

# Output: configuration
out += "\n~~~\nSearch criteria:\n"
out += "\nd atom: " + d_atom_type
out += "\nligand: " + ligand_type
if d_valence is not None:
    out += "\nvalence of the d atom: " + d_valence
out += ("\nmaximal distance between a d atom and a ligand: {:4.3f}"
        .format(max_d_to_ligand_distance))
out += ("\nmaximal distance between d atoms: {:4.3f}"
        .format(max_d_to_d_distance))
out += ("\nminimal number of common ligands: {:d}"
        .format(min_number_of_common_ligands))
out += ("\nminimal d-lingand-d angle: {:5.2f}".format(min_angle))
out += ("\nmaximal d-lingand-d angle: {:5.2f}\n~~~".format(max_angle))

# Output: header of the final table
out += "\n\n| compound | cif file | d-ligand-d angles |"
out += "\n|-----------------------------------------|---------|"
out += max_angles_on_single_line*"-------"
out += "|"

ciffiles = [filename for filename in os.listdir(cifdir)
            if filename.find('.cif') != -1]

if len(ciffiles) == 0:
    raise Exception("ERROR: no cif files found.")

err = "\n~~~\nSearching in a pool of {:d} cif files.".format(len(ciffiles))
err += "\n\nWarnings:\n"

for ciffile in ciffiles:
    try: structure = mg.Structure.from_file(os.path.join(cifdir, ciffile))
    except:
        err += "\nStructure {:s} is corrupt.".format(ciffile)
    else:
        site_indices \
            = [site_index for site_index in range(len(structure))
               if (structure[site_index]\
                       .species.elements[0].symbol == d_atom_type)
               and (check_if_site_is_ordered(structure, site_index))]
        if d_valence is not None:
            site_indices = [site_index for site_index in site_indices
                           if (structure[site_index]\
                                  .species.elements[0].oxi_state == d_valence)]
        bridging_angles = []
        for site_index in site_indices:
            central_atom = np.array(structure[site_index].coords)
            neighbors = get_environment(structure,
                                        site_index,
                                        max_d_to_d_distance,
                                        d_atom_type,
                                        d_valence)
            ligands = get_environment(structure,
                                      site_index,
                                      max_d_to_ligand_distance,
                                      ligand_type)
            for neighbor in neighbors:
                candidate_angles = []
                for ligand in ligands:
                    if norm(neighbor-ligand) < max_d_to_ligand_distance:
                        candidate_angles.append(get_angle(ligand-neighbor,
                                                          ligand-central_atom))
                if len(candidate_angles) >= min_number_of_common_ligands:
                    bridging_angles += candidate_angles

        if len(bridging_angles) != 0:
            angles = np.unique(np.around(np.array(bridging_angles),
                                         decimals=2))
            indices_of_relevant_angles \
                = np.array([angle_index for (angle_index, angle)
                            in enumerate(angles)
                            if (angle >= min_angle) and (angle <= max_angle)])
            if len(indices_of_relevant_angles) != 0:
                formula = format_formula(structure.composition.reduced_formula)
                out += "\n| {:s} | {:s} |".format(formula, ciffile)
                angles_on_current_line = 0
                for (angle_index, angle) in enumerate(angles):
                    if angles_on_current_line == max_angles_on_single_line:
                        out += "|\n|      |      |"
                        angles_on_current_line = 0
                    out += (" {:5.2f}".format(angle),
                            " **{:5.2f}**".format(angle))\
                           [angle_index in indices_of_relevant_angles]
                    angles_on_current_line += 1
                out += " |"

out += "\n" + err + "\n~~~"

with open(prefix + '.md', 'w') as mdfile:
    mdfile.write(out)

output = pd.convert_file(prefix + '.md', to='pdf', format='md',
                         outputfile=(prefix + '.pdf'),
                         extra_args=['--template=eisvogel.tex'])
