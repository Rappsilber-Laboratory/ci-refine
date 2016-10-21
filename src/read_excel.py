import xlrd
import numpy as np
import cPickle
from scipy.spatial.distance import euclidean
from sklearn import mixture

def quasi_distance(coords_1, coords_2):
    return sum([abs(i-j) for i, j in zip(coords_1, coords_2)])

def max_consecutive(chromatogram, val = 0.5):
    consecutive_count = 0
    max_consecutive_count = 0
    for value in chromatogram:
        if value > val:
            consecutive_count += 1
            if consecutive_count > max_consecutive_count:
                max_consecutive_count = consecutive_count
        else:
            consecutive_count = 0
    #print chromatogram
    return max_consecutive_count



def gmm_fit(chromatogram):
    best_com = 0
    best_bic = 999
    num_gauss_map = {3:1,
                     4: 1,
                     5: 1,
                     6: 1,
                     7: 1,
                     8: 1,
                     9: 2,
                     10: 2,
                     11: 2,
                     12: 3,
                     13: 3,
                     14: 3,
                     15: 4,
                     16: 5,
                     17: 5,
                     18: 5,
                     19: 5,
                     20: 5,
                     21: 5,
                     22: 5,
                     23: 5,
                     24: 5,
                     25: 5,
                     26: 5,
                     27: 5,
                     28: 5,
                     29: 5,
                     30: 5,
                     31: 5,
                     32: 5,
                     33: 5,
                     34: 5,
                     35: 5,
                     36: 5,
                     37: 5,
                     38: 5,
                     39: 5,
                     40: 5,
                     41: 5,
                     42: 5,
                     43: 5,
                     44: 5,
                     45: 5,
                     46: 5,
                     47: 5,
                     48: 5,
		     49: 5,
		     50: 5}

    for n_com in xrange(1, num_gauss_map[max_consecutive(chromatogram)]+1):
        gmm = mixture.GMM(n_components=n_com, covariance_type='full')
        gmm.fit(chromatogram[:, np.newaxis])
        bic = gmm.bic(chromatogram[:, np.newaxis])
        if bic < best_bic:
            best_bic = bic
            best_com = n_com
    gmm = mixture.GMM(n_components=best_com, covariance_type='full')
    gmm.fit(chromatogram[:, np.newaxis])
    means = np.round(gmm.means_, 3) # means of gaussians
    weights = np.round(gmm.weights_, 3) # weight of gaussians
    covars = np.round(gmm.covars_, 3) # covariances
    max_mean = 0
    max_w = 0
    max_covar = 0
    """Get the mean, covar of the gaussian with highest weight"""
    for m,w,c in zip(means,weights,covars):
        if w > max_w:
            max_w = w
            max_mean = m
            max_covar = c
        #print max_mean[0], max_covar[0],  gmm.bic(np.array(values))
    return max_mean[0]


def check_chromatogram(chromatogram, val = 0.0, con_count=5):
    consecutive_count = 0
    for value in chromatogram:
        if value > val:
            consecutive_count += 1
            if consecutive_count == con_count:
                # print chromatogram
                return True
        else:
            consecutive_count = 0

    #print chromatogram
    return False

workbook = xlrd.open_workbook("../data/protein_interaction_data/chromatogram_data.xls")
sh = workbook.sheet_by_name("Exp1")
# print sh.col_values(0)

for row in range(0, 1): # Exp 2, 3
#for row in range(29, 30): # Experiment 1
    for i in range(len(sh.row_values(row))):
        print i, sh.row_values(row)[i]

relevant_cols = [i for i in range(61, 109)] # Exp 1
#relevant_cols = [i for i in range(58, 103)] # Exp 2
#relevant_cols = [i for i in range(63, 113)] # Exp 3
print relevant_cols

data = []
for row in range(30, 1991): # Exp1
#for row in range(2, 2222): #Exp2
#for row in range(2, 2938): # Exp3
    uniprot = sh.row_values(row)[4]
    chromatogram_data = []

    for i in relevant_cols:
        if isinstance(sh.row_values(row)[i], float):
            chromatogram_data.append(sh.row_values(row)[i])
        else:
            chromatogram_data.append(0.0)

    if check_chromatogram(chromatogram_data, val=0.0, con_count=5) and check_chromatogram(chromatogram_data, val=0.5, con_count=3):

        if uniprot.count(";") == 0 and len(uniprot) > 1:
            #print uniprot
            data.append((str(uniprot), np.array(chromatogram_data)))#, gmm_fit(np.array(chromatogram_data))))
        else:
            #print uniprot
            protein_ids = [str(i.split("-")[0]) for i in uniprot.split(";") if len(str(i.split("-")[0])) > 1]
            #print protein_ids
            #print
            if len(list(set(protein_ids))) >= 1:
                data.append((list(set(protein_ids))[0], np.array(chromatogram_data)))#, gmm_fit(np.array(chromatogram_data))))
            #print (list(set(protein_ids)))                #print(list(set(protein_ids)))[0]



interaction_data = []
for i in range(0, len(data)):
    for j in range(i, len(data)):
        if i < j:
            euclidean_distance = euclidean(data[i][1], data[j][1])#euclidean(data[i][1], data[j][1])
            if euclidean_distance <= 4.0: #and abs(data[i][2]-data[j][2]) < 0.5:
                interaction_data.append((1.0-euclidean_distance/4.0, data[i][0], data[j][0]))
interaction_data.sort(reverse=True)

with open('interaction_data_exp_1_eu_4_final.pkl', 'wb') as outfile:
    cPickle.dump(interaction_data, outfile, cPickle.HIGHEST_PROTOCOL)

file = open("protein_list.txt", "w")
proteins = [protein for protein, _ in data]
proteins = list(set(proteins))
for protein in proteins:
    file.write(protein + '\n')
file.close()

# print sh.nrows
# print sh.ncols
