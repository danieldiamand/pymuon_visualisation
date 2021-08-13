
import ase
from ase import io
import pathlib
import molgif
import yaml

yaml_settings_file_name = "MnSi-muairss-castep.yaml"
clustered_muons_file_name = "MnSi_MnSi_uep_clusters.dat"
relaxed_file_extension = ".xyz"
cluster_report_file_name = "MnSi_clusters.txt"
is_each_muon_different_atom = False #will set each muon to a different atom from atom_list depending on which cluster they join in relaxed and clustered stage.
is_create_gif = True #will look bad if is_each_muon_different_atom is true (this could be fixed but is a small issue)
atom_list = ["H","He","Ne","Ar", "Kr", "Xe"]

#loads yaml file as dict
def load_yaml(yaml_file_name):
    return(yaml.load(open(yaml_file_name)))

#fetches name from yaml file and opens file_name-out.cell with ase
def load_out_cell(yaml_settings_file_name):
    file_name = load_yaml(yaml_settings_file_name)["name"]
    return (ase.io.read(file_name+"-out.cell"))

#returns array of molecule [x,y,z] bounds
def get_lattice_params(molecule):
    return(molecule.get_cell().lengths())


#reads clustering report and generates 2d array where each array contains a list of all muon names of a cluster
def cluster_report_to_structures_array(cluster_report_file_name):
    with open(cluster_report_file_name, "r") as cluster_report:
        cluster_report_array= cluster_report.read().split("\tStructure list:\n")
    del cluster_report_array[0]
    structures = []
    for clusters in cluster_report_array:
        clusters = clusters.split("\n\n")
        structure = clusters[0].replace("\n","\t").split("\t")
        if structure != "": 
            structures.append(clusters[0].replace("\n","\t").split("\t"))
    return(structures)

#sets muon to differnt atom, depending on which array id it matches in clustering_array
def muon_unique_atom_from_clustering_array(muon, path, atom_list):
    clustering_array = cluster_report_to_structures_array(cluster_report_file_name)
    for i in range(1,len(clustering_array)):
        if (path.parts[-1] in clustering_array[i]):
            muon.symbol = atom_list[i]
    return(muon)


#turns molecule files into rotating gifs
def mol_to_gif(molecule, out_file_name):
    mol_gif = molgif.Molecule(molecule)
    mol_gif.remove_bonds()
    mol_gif.name = out_file_name
    mol_gif.save_rot_gif()
    
# fetches the random muon position from each folder and puts it into a single file 
def combine_random_muons(yaml_settings_file_name):
    out_file_name = "combined_random_muons.cif"
    yaml_settings_dict = load_yaml(yaml_settings_file_name)
    new_molecule = load_out_cell(yaml_settings_file_name)
    for path in pathlib.Path(yaml_settings_dict["out_folder"]+"/"+yaml_settings_dict["calculator"]).iterdir():
        if path.is_dir():
            muon = ase.Atom("H", load_yaml(str(path)+"/"+path.parts[-1]+".yaml")["mu_pos"])
            new_molecule.append(muon)
    ase.io.write(out_file_name, new_molecule)
    if (is_create_gif):
        mol_to_gif(new_molecule, out_file_name)


#fetches the relaxed muon positions, ensures they are within the bounds 
def combine_relaxed_muons(yaml_settings_file, relaxed_file_extension, cluster_report_file_name, is_each_muon_different_atom, atom_list):
    out_file_name = "combined_relaxed_muons.cif"
    yaml_settings_dict = load_yaml(yaml_settings_file_name)
    new_molecule = load_out_cell(yaml_settings_file_name)
    lattice_params = get_lattice_params(new_molecule)
    for path in pathlib.Path(yaml_settings_dict["out_folder"]+"/"+yaml_settings_dict["calculator"]).iterdir():
        if path.is_dir():
            muon = ase.io.read(str(path)+"/"+path.parts[-1]+relaxed_file_extension)[-1]
            for i in range(3):
                if muon.position[i] < 0:
                    muon.position[i]=muon.position[i]+lattice_params[i]
                elif muon.position[i] > lattice_params[i]:
                    muon.position[i]=muon.position[i]-lattice_params[i]
            if (is_each_muon_different_atom):
                muon = muon_unique_atom_from_clustering_array(muon, path, atom_list)
            new_molecule.append(muon)
    ase.io.write(out_file_name, new_molecule)
    if (is_create_gif):
        mol_to_gif(new_molecule, out_file_name)

#reads and combines the clustered muons into one file
def combine_clustered_muons(yaml_settings_file_name, clustered_muons_file_name, cluster_report_file_name, is_each_muon_different_atom, atom_list):
    out_file_name = "combined_clustered_muons.cif"
    new_molecule = load_out_cell(yaml_settings_file_name)
    with open(clustered_muons_file_name, "r") as clustered_muon_file:
        clustered_muon_array = clustered_muon_file.read().split("\n")
    for i in range(len(clustered_muon_array)-1):
        muon_array = clustered_muon_array[i].split("\t")
        muon = (ase.Atom("H",[muon_array[5], muon_array[6], muon_array[7]]))    
        if (is_each_muon_different_atom):
            muon.symbol = atom_list[i]
        new_molecule.append(muon)
    ase.io.write(out_file_name, new_molecule)
    if (is_create_gif):
        mol_to_gif(new_molecule, out_file_name)

if __name__ == "__main__":        
    combine_random_muons(yaml_settings_file_name)
    combine_relaxed_muons(yaml_settings_file_name, relaxed_file_extension, cluster_report_file_name, is_each_muon_different_atom, atom_list)
    combine_clustered_muons(yaml_settings_file_name, clustered_muons_file_name, cluster_report_file_name, is_each_muon_different_atom, atom_list)
