import math

import numpy
import Bio.PDB
import matplotlib.pyplot as plt

"""It might be a little confusing that ContactMap class is currently used to calculate contact maps AND distance maps. It is mostly done for efficiency, since 
   one needs to compute the distance map to derive the contact map anyway.... Maybe it should be split into two classes that inherent from a base"""

class ContactMap:
  def __init__(self, contact_threshold = 8):
    self.contact_map = None
    self.distance_map = None
    self.cutoff = contact_threshold
    
    self.resid_map = {}
    
    self.tp_predictions = None
    self.fp_predictions = None 

  def load_cm_from_pdb(self, id, chain_id, pdb_file ):
    structure = Bio.PDB.PDBParser().get_structure(id, pdb_file)
    model = structure[0]
    self.calc_contact_matrix(model[chain_id], model[chain_id])
 
  def calculate_cm_from_struct(self, bio_pdb_struct, input_seqsep ):
    self.calc_contact_matrix(bio_pdb_struct, bio_pdb_struct, input_seqsep )
  
  def load_tp_predictions(self, contact_list):
    if len(contact_list) >= 1:
        self.tp_predictions = contact_list

  def load_fp_predictions(self, contact_list):
    if len(contact_list) >= 1:
        self.fp_predictions = contact_list

  def calc_residue_dist(self,residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    atom_res1 = ''
    atom_res2 = ''
    if residue_one.get_resname() == "GLY":
        #print residue_one.get_resname()
        atom_res1 = "CA"
    else:
        atom_res1 = "CB" # CB here for CB definition!
    
    if residue_two.get_resname() == "GLY":
        #print residue_one.get_resname()
        atom_res2 = "CA"
    else:
        atom_res2 = "CB" # CB here for CB definition!
    
    try:
        diff_vector  = residue_one[atom_res1].coord - residue_two[atom_res2].coord
    except:
        diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))
        
    
  def calc_contact_matrix(self, chain_one, chain_two, seqsep = 12):
    """Returns a matrix of C-beta distances (C-alpha wheren residue is glycine) between two chains"""
    counter = 0
    row_count = 0
    col_count = 0
    self.contact_map = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    self.distance_map = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for row, residue_one in enumerate(chain_one) :

        for col, residue_two in enumerate(chain_two) :

            if (residue_two.get_resname() == 'HOH'):
                raise "PDB File contains HETATM entries. This is not currently supported. Delete HETATM lines in pdb."

            dist = self.calc_residue_dist(residue_one, residue_two)

            r_id_1 = int(residue_one. __repr__().split()[3].split('=')[1])
            r_id_2 = int(residue_two. __repr__().split()[3].split('=')[1])

            self.resid_map [( r_id_1,r_id_2)] = (row,col)
            self.distance_map [ row, col ] = dist

            #if dist <= float(self.cutoff) and abs(col_count - row_count) < 24 and abs(col_count - row_count) >= 12 :
            if dist <= float(self.cutoff) and abs(row - col) >= seqsep :

                self.contact_map [row, col] = 1.0
                counter += 1
            else:

                self.contact_map [row, col] = 0
            col_count += 1
        col_count = 0
        row_count += 1
  
  def number_of_contacts(self, seqsep_min = 6, seqsep_max = 999, resnum = None):
      num_contacts = 0
      for resnum_lower, residue_lower in enumerate(self.contact_map):
          for resnum_upper, residue_upper in enumerate(residue_lower):
              if resnum_lower < resnum_upper and self.contact_map[resnum_lower][resnum_upper] == 1 and abs(resnum_upper - resnum_lower) >= seqsep_min and abs(resnum_upper - resnum_lower) <= seqsep_max: 
                  if resnum == None:
                      num_contacts += 1
                  else:
                      if (resnum_lower + 1 == resnum):
                          num_contacts += 1
      return num_contacts


  def calc_centroid(self, chain_one):
      x = []
      y = []
      z = [] 
      len_chain = float(len(chain_one))
      for residue in chain_one:
          atom_res1 = ''
          if residue.get_resname() == "GLY":
            atom_res1 = "CA"
          else:
            atom_res1 = "CA" # CB here for CB definition!
          x.append( residue[atom_res1].coord[0])
          y.append( residue[atom_res1].coord[1])
          z.append( residue[atom_res1].coord[2])
      return [numpy.mean(x),numpy.mean(y),numpy.mean(z)]
  
  def print_cm(self):
    print self.contact_map
    
  def print_dm(self):
    print self.distance_map

  def contact_list(self):
      list = []
      for resnum_lower, residue_lower in enumerate(self.contact_map):
          for resnum_upper, residue_upper in enumerate(residue_lower):
              if resnum_lower <= resnum_upper and self.contact_map[resnum_lower][resnum_upper] == 1:
                  list.append((resnum_lower, resnum_upper, 1.0))
      return list

  def distances(self, seq_sep, include_residue_numbers):
      d_vector = []
      for resnum_lower, residue_lower in enumerate(self.distance_map):
          for resnum_upper, residue_upper in enumerate(residue_lower):
              if resnum_lower <= resnum_upper and abs(resnum_lower - resnum_upper) >= seq_sep:
                  if include_residue_numbers:
                      d_vector.append((resnum_lower+1, resnum_upper+1, self.distance_map[resnum_lower][resnum_upper]))
                  else:
                      d_vector.append(self.distance_map[resnum_lower][resnum_upper])
      return d_vector
      
  
  def draw_cm(self):
    draw_map = self.contact_map
    draw_map = draw_map *0.2
    draw_map[0][0] = 1
   
    m = plt.imshow(draw_map,cmap="binary", interpolation="nearest")
    m = plt.scatter(draw_map[0][:], draw_map[:][0],c="k", s="40")
    
    if self.tp_predictions != None:
        x_list = []
        y_list = []
        for x, y in self.tp_predictions:
            x_list.append(y-1)
            y_list.append(x-1)
        m = plt.scatter(x_list,y_list, c='#B40426', s=100,  linewidths=0.0, marker = (5,1))
        m = plt.scatter(x_list,y_list, c='#3B4CC0', s=60,  linewidths=0.0, marker = (5,1))
    
    if self.fp_predictions != None:
        x_list = []
        y_list = []
        for x, y in self.fp_predictions:
            x_list.append(y-1)
            y_list.append(x-1)
        m = plt.scatter(x_list,y_list, c='#B40426', s=100,  linewidths=0.0, marker = (5,1))
        m = plt.scatter(x_list,y_list, c='#B40426', s=60,  linewidths=0.0, marker = (5,1))
    
    plt.show()

  def get_mapped_distance(self, i, j):
        if self.resid_map.has_key((i,j)):
            mapped_tuple = self.resid_map[(i,j)]
            return self.distance_map[mapped_tuple[0]][mapped_tuple[1]]
        else:
            return 999.9

  def get_cm(self):
    return self.contact_map

  def get_dm(self):
    return self.distance_map

  def get_contacts(self, max_distance):
    distances = self.distances(12, True)
    for distance in distances:
        print distance[2]
    return 'a'