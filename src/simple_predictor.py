import numpy as np
from sklearn import cross_validation
from sklearn import datasets
from sklearn import svm
import cPickle
import sys
sys.path.append("experiments/pagerank_additional_test_sets/")
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import BernoulliNB
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_score
from sklearn.svm import LinearSVC
from scipy.sparse import coo_matrix, hstack, vstack
from sklearn.cross_validation import KFold
import argparse
from evaluation_helper import load_pdb_ids
import random
import itertools
from tune_pagerank_parameters import load_native_map, load_files, run_pagerank_classifier, evaluate_ranking
from compute_co_features import co_occurance_features_easy, co_occurence_features
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import Normalizer, StandardScaler, PolynomialFeatures
from sklearn import svm

def parse_arguments():
    """Specify and parse command line inputs
    """
    global options

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="contact_file_list", help="PDB ID list", required=True)
    parser.add_argument("-p", dest="prediction_folder", required=True)
    parser.add_argument("-s", dest="secondary_structure_folder", required=True)
    parser.add_argument("-m", dest="method_name", required=True)
    options = parser.parse_args()


def dataset_splits(protein_data, fold):
    # deepcopy protein data list
    tmp_protein_data = list(protein_data)
    # shuffle dataset
    random.shuffle(tmp_protein_data)
    # Get folds
    kf = KFold(len(tmp_protein_data), fold)
    return kf

    #index_train = int(float((fold-1))/float(fold) * len(protein_data))
    #x_train = protein_data[:index_train]
    #x_test = protein_data[index_train+1:]
    #return x_train, x_test


def load_feature_data(pdb_data_list, feature_folder):
    """
    Parameters
    ----------
    protein_data : List of proteins that you want to load features
    feature_folder : Folder where the features are stored

    Returns
    -------
    Dictionary : key->pdb_id, value->tuple(features, labels)
    """

    data_dict = {}

    for training_protein, training_protein_length in pdb_data_list:
        labels = np.genfromtxt("%s%s_co_occurance_labels.dat" % (feature_folder, training_protein))
        with open('%s%s_sparse_array_features.dat' % (feature_folder, training_protein), 'rb') as infile:
            features = cPickle.load(infile)
        data_dict[training_protein] = (features, labels)
        #training_labels = list(itertools.chain(list(training_labels), list(labels)))

        #training_features = vstack((training_features, features))
    #return training_features, training_labels
    return data_dict


def stack_features(protein_data, data, train_indices, test_indices):


    numpy_protein_data = np.array(protein_data)
    train_proteins = numpy_protein_data[train_indices]
    test_proteins =numpy_protein_data[test_indices]

    counter = 0
    for protein in train_proteins:
        features, labels = data[protein[0]]
        if counter == 0:
            training_features = features
            training_labels = labels
        else:
            training_features = vstack((training_features, features))
            training_labels = list(itertools.chain(list(training_labels), list(labels)))
        counter +=1

        """
        counter = 0
        for protein in test_proteins:
            features, labels = data[protein[0]]
            if counter == 0:
                test_features = features
                train_labels = labels
            else:
                test_features = vstack((training_features, features))
                train_labels = list(itertools.chain(list(training_labels), list(labels)))
            counter += 1
        """
    return training_features, training_labels, test_proteins


def find_pagerank_parameters(test_proteins, alpha_range, beta_range, clf, shift_dict, all_values, clf_parameters, norm):
    #print len(test_proteins)
    print "K-fold PageRank parameter tuning"
    for pdb_id, length in test_proteins:

        contacts, sec_struct = load_files(pdb_id, sec_struct_folder=options.secondary_structure_folder,
                                          pred_folder=options.prediction_folder)
        native_map = load_native_map(pdb_id)
        feature_matrix = []
        contact_pairs = []

        co_occurance_features_easy(contacts[0:int(int(length) * 2.5)], sec_struct, shift_dict, feature_matrix,
                                   contact_pairs)
        X = norm.transform(feature_matrix)
        probs = clf.predict_proba(X)
        #print probs
        #sum_probs = 0
        #for i in probs:
        #    sum_probs += i[1]
        contact_score_dict = {}

        count = 0
        class_indicator = 1
        for i in clf.classes_:
            if i == 1.:
                class_indicator = count
            count += 1
        print class_indicator
        for i in xrange(0, len(contact_pairs)):
            contact_score_dict[contact_pairs[i]] = probs[i][class_indicator] #/ sum_probs
        #for keys, values in contact_score_dict.iteritems():
        #    print values
        #contact_list = []
        #for keys, values in contact_score_dict.iteritems():
        #    contact_list.append((keys[0], keys[1], values))

        #contact_list.sort(reverse=True, key=lambda x: x[2])
        #print contact_list
        for alpha in alpha_range:
            for beta in beta_range:
                xl_ranked = run_pagerank_classifier(contacts, sec_struct, alpha, beta, shift_dict, length,
                                                    contact_score_dict)

                acc = evaluate_ranking(xl_ranked, native_map, int(float(length) * 0.2))
                all_values[alpha, beta, clf_parameters].append(acc)
                #print alpha, beta, acc

def best_parameters(all_values):
    parameter_ranking = []
    for keys, values in all_values.iteritems():
        parameter_ranking.append((np.mean(values), keys))
    parameter_ranking.sort()
    parameter_ranking.reverse()
    return parameter_ranking


#labels = np.genfromtxt("co_occurance_labels")
#with open('test_sparse_array.dat', 'rb') as infile:
#    features = cPickle.load(infile)
#features_dense = features.toarray()

#for i in xrange(0, 5):
#    X_train, X_test, y_train, y_test = cross_validation.train_test_split(features_dense, labels, test_size=0.4)
#    clf = BernoulliNB().fit(X_train, y_train)
#    x_predict = clf.predict(X_test)
#    print precision_score(y_test, x_predict)
    #print clf.score(X_test, y_test)

#with open('easy_model.clf', 'wb') as outfile:
#    cPickle.dump(clf, outfile, cPickle.HIGHEST_PROTOCOL)

def main():

    parse_arguments()
    protein_data = load_pdb_ids(options.contact_file_list)

    # Initialize and load data
    shift_dict = cPickle.load(open("probabilities/shifts_sigma_0.05.txt", "rb"))
    feature_folder = "/scratch/schneider/projects/pagerank_refinement/data/co-occurence_features_easy_both/"
    protein_data = protein_data
    data = load_feature_data(protein_data, feature_folder)
    kfold_split_iterator = dataset_splits(protein_data, 2)


    # PageRank Parameter tuning
    alpha_range = [0.50]#, 0.5, 0.6]#, 0.5,0.6]# , 0.35, 0.4, 0.45, 0.5]#, 0.5, 0.6]#, 0.6]#, 0.6]#, 0.6]
    beta_range = [2.5]#, 2.5]#, 2.5]#, 2.5]#, 2.0]#, 2.5]#, 2.5]#, 2.5, 3.0]#, 2.5]
    #cost_range = [0.0001]#, 0.01]
    cost_range = [0.01]#, 1e-5]#, 0.00001, 0.0001]#, 0.00001, 0.0001, 0.001]
    print alpha_range
    print beta_range
    all_values = {}

    for alpha in alpha_range:
        for beta in beta_range:
            for c in cost_range:
                all_values[(alpha, beta, c)] = []

    for train_indices, test_indices in kfold_split_iterator:
        print "Constructing data set"
        training_features, training_labels, protein_test = stack_features(protein_data, data, train_indices,
                                                                            test_indices)
        training_features = training_features.toarray()
        print "Training classifier with %s samples"%len(training_labels)
        for c in cost_range:
            norm = StandardScaler().fit(training_features)
            scaled_training = norm.transform(training_features)
            #alpha = 1.0
            clf = LogisticRegression(C=c).fit(scaled_training, training_labels)
            #clf = svm.SVC(C=c, kernel='linear', probability=True).fit(training_features, training_labels)
            #clf = RandomForestClassifier(n_estimators=c, min_samples_split=500, min_samples_leaf=500).fit(training_features, training_labels)
            #clf = KNeighborsClassifier(n_neighbors=c).fit(scaled_training, training_labels)
            print "Tuning PageRank parameters cost %s" % str(c)
            find_pagerank_parameters(protein_test, alpha_range, beta_range, clf, shift_dict, all_values, c, norm)



    best_c = best_parameters(all_values)[0][1][2]
    print "Best paramters are:", best_parameters(all_values)[0]

    # train classifier with full dataset parameters
    train_indices = [i for i in xrange(0, len(protein_data))]
    test_indices = []
    training_features, training_labels, protein_test = stack_features(protein_data, data, train_indices,
                                                                      test_indices)
    training_features = training_features.toarray()
    print "Train classifier with full set"
    #clf = KNeighborsClassifier(n_neighbors=best_c).fit(training_features, training_labels)
    norm = StandardScaler().fit(training_features)
    scaled_training = norm.transform(training_features)

    clf = LogisticRegression(C=best_c).fit(scaled_training, training_labels)
    #clf = RandomForestClassifier(n_estimators=best_c, min_samples_split=500, min_samples_leaf=500).fit(training_features, training_labels)
    #print clf.feature_importances_

    with open('easy_model_dummy.clf', 'wb') as outfile:
        cPickle.dump(clf, outfile, cPickle.HIGHEST_PROTOCOL)
    with open('normalizer', 'wb') as outfile:
        cPickle.dump(norm, outfile, cPickle.HIGHEST_PROTOCOL)

        # Tune alpha/beta parameters

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")
