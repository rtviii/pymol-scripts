import pymol
from pymol import cmd, selector
import os
import shutil
import numpy as np


# Create object for every chain present:
# Some [namespace] issues with this one.
# def objectify_chains(obj_name):
#     chains_array = cmd.get_chains(obj_name)
#
#     print ("Chains present in {on}: ".format(on=obj_name))
#     print(" ".join(map(str, ["[" + str(chain) + "]" for chain in chains_array])))
#
#     for chain in chains_array:
#         cmd.create("chain" + "_" + chain, "chain {}".format(str(chain)))
#         cmd.zoom()
#     cmd.delete(obj_name)
#

def see(chain_name):
    cmd.select("chain " + chain_name)
    cmd.color("gray20", "all and (not chain " + chain_name + " )")
    cmd.color("cyan", "chain " + chain_name)
    cmd.deselect()


def nbr_prompt(fulcrum, nbr):
    return "☆" if str(nbr) == str(fulcrum) else "[" + nbr + "]"


def get_nbrs(chain_name, radius):
    cmd.set("specular_intensity", 0.1)
    # see(chain_name)

    cmd.select("diapason", "all within {r} of chain {c}".format(r=radius, c=chain_name))
    #   extrapolate chains present in the selection
    cmd.select("withinChainsOf{n}".format(n=chain_name), "bc. diapason")

    neighbours = cmd.get_chains("withinChainsOf{n}".format(n=chain_name))

    cmd.color("gray10", "all");
    cmd.color("green", "chain " + chain_name);
    for nbr in neighbours:
        if nbr != chain_name:
            cmd.color("white", "chain " + nbr)

    cmd.delete("diapason")
    cmd.delete("withinChainsOf{n}".format(n=chain_name))

    print (
        chain_name + " <--> ☆ | Present in", radius,
        "Å vicinity: "
        + " ".join(map(str, [nbr_prompt(chain_name, x) for x in neighbours])))


def save_xyz(chain_id):
    cmd.create("chain_" + chain_id, "chain {}".format(str(chain_id)))
    coordinate_set = cmd.get_coordset("chain_{}".format(str(chain_id)), 1)
    print("Got ", chain_id, "'s coordinates: ", coordinate_set)
    filename_string = 'xyz_chain_{name}'.format(name=chain_id)
    np.save(filename_string, coordinate_set, allow_pickle=True)
    move_to_coords(filename_string)
    print ("Coordinates saved successfully under ./Coordinates/", filename_string)


# To customize on further specification
# extension should be omitted
def move_to_coords(filename):
    filepath = filename + ".npy"
    cwd = os.getcwd()
    coord_path = "Coordinates"
    # Create target directory & all intermediate directories if don't exists
    try:
        os.makedirs(coord_path)
    except OSError:
        pass

    shutil.move(cwd + "\\" + filepath, cwd + "\\" + coord_path + "\\" + filepath)


def nbrmenu():
    print("___________________________________________________________________________________")
    print("|              SIGNATURE              |              EFFECT          \n")
    print("|get_nbrs [chain_identifier],[radius] | returns array of neighbor-chain identifiers  |\n")
    print("|see               [chain_identifier] | color-highlights a chain by identifier       |\n")
    # print("|objectify chains       [object_name] | partitions molecule into chain-objects       |\n")
    print("|save_xyz          [chain_identifier] | saves coordinates as .npy into /Coordinates  |\n")
    print("|nbrmenu                              | see this tableau again                       |")
    print("___________________________________________________________________________________\n")

cmd.extend("get_nbrs", get_nbrs)
# cmd.extend("objectify_chains", objectify_chains)
cmd.extend("see", see)
cmd.extend("save_xyz", save_xyz)
cmd.extend("nbrmenu", nbrmenu)

print ("Commands have been added:")
nbrmenu()
