import glob
import shutil, os, sys, re
import subprocess
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats
from scipy.spatial import distance

from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.cluster import HDBSCAN
from sklearn.preprocessing import minmax_scale

import MDAnalysis as mda
from MDAnalysis.analysis.dssp import DSSP

import pymol


class blind_screening():
    def cluster_structures(X):
        """
        loop through values of k and define best value of k with silhouette_score
    
        Input: 
            X : np.ndarray (n, m) | result of PCA
    
        Output: 
            cluster_labels : (n, 1) | list of optimal clusters for X
        """
    
        k_range = range(2,51)
        sil_score = []
        for k in k_range:
            clustering = HDBSCAN(min_cluster_size=k,min_samples=1)
            clustering.fit(X)
            if len(set(clustering.labels_)) > 1 and len(set(clustering.labels_)) < len(X):
                score = silhouette_score(X, clustering.labels_, metric='euclidean')
                sil_score.append(score)
            else:
                sil_score.append(-1)
    
        opt_k = k_range[np.argmax(sil_score)]
        clustering = HDBSCAN(min_cluster_size=opt_k)
        clustering.fit(X)
        return clustering.labels_
    
    def k_medoids(X, l, labels, k=3, max_iter=100):
        """
        K-Medoid algorithm to find suitable representative structures from each cluster defined by HDBSCAN.
    
        Input:
            X:        np.ndarray (n, m)  | all points from one HDBSCAN cluster
            k:        number of medoids  | 
            max_iter: maximum number of iterations allowed to minimize the distance
            l:        current HDBSAN label
            labels:   full list of HDBSCAN labels
    
        Output:
            medoids:     indices of the K medoids
            total_cost:  sum of distances of each point to its medoid
        """
        np.random.seed(42)
    
        #start with random k points
        temp = X.copy()
        mask = np.zeros(X.shape, dtype=bool)
        mask[np.argwhere(labels == l)] = True
        
        #check the number of points in a cluster
        #if less than 4 just return those indices
        _, cluster_count = np.unique(mask[:,0], return_counts=True) # count = False, True
        cluster_count = cluster_count[[idx for idx, val in enumerate(_) if val == True][0]]#<-account for the case of one cluster
    
        if cluster_count < 4:
            return np.ravel(np.argwhere(mask[:,0] == True)), np.nan
        # block out values that are not within the current HDBSCAN group
        temp[~mask] = 9999
    
        number_samples = temp.shape[0]
        medoids = np.random.choice(number_samples, k, replace=False)
    
        #distance matrix of randomly chosen points
        D = distance.cdist(temp, temp[medoids], metric='euclidean')
        tot_cost = np.sum(np.min(D, axis=1))
    
        itr = 0
        while itr < max_iter:
            reduced = False
    
            #loop through all possibilities
            for m_idx in range(k):
                for current_idx in range(number_samples):
                    if current_idx in medoids:
                        continue
    
                    new_medoids = medoids.copy()
                    new_medoids[m_idx] = current_idx
    
                    #new distance matrix
                    D_new = distance.cdist(temp, temp[new_medoids], metric='euclidean')
                    new_cost = np.sum(np.min(D_new, axis=1))
    
                    #if the cost has been reduced move onto the the next sample
                    if new_cost < tot_cost:
                        medoids = new_medoids
                        tot_cost = new_cost
                        reduced = True
                        break
                if reduced:
                    break
    
            if not reduced:
                #If there was no improvement we should be converged
                break
            itr+=1
        return medoids, tot_cost


    def __init__(self, pdb1_name, blind_path):
   # def main():
        """
        requires Foldseek and Pymol
    
        Find all pdb files from CF-Random generated directories.
        This script will automatically generate a Foldseek database of these structures 
        then calculate a similarity matrix of all structures based on bit-score.
        similarity matrix -> PCA -> HDBSCAN -> K-medoids -> structures of interest.
    
        The final output is then a png file showing the result of PCA and HDBSCAN
        a text file containing the coordinates of the structures of interest, file name, and group ID
        finally this script will automatically generate a pse file of the structures_of_interest
        """

    
        #_______________collect all pdb files that CF-Random generated_____________________________
        db_directory = blind_path +  "/pdbs_for_db/"
        #db_directory =  "/pdbs_for_db/"
        if not os.path.isdir(db_directory):
            os.mkdir(db_directory)
        #pdb_files = glob.glob("./**/*.pdb", recursive=True)
        pdb_files = glob.glob(blind_path + "/**/*.pdb", recursive=True)
        pdb_files = [file for file in pdb_files if db_directory not in file]
        print("Gathering pdb pdb files for self-search")
        for file in pdb_files:
            dest_name = file.replace('/','-')
            if not os.path.isfile(db_directory + dest_name[17:]):
                shutil.copyfile(file, db_directory + dest_name[17:])
        #__________________________________________________________________________________________
    
    
        print("Creating database...")
        create_db = ["foldseek", "createdb", db_directory, db_directory + "DB"]
        if not os.path.isfile(db_directory + "DB"):
            try:
                response = subprocess.run(create_db, capture_output=True, text=True, check=True )
            except subprocess.CalledProcessError as e:
                print("ERROR:\n", e.stderr)
            
            print('Succes database is up!')
        else:
            print("found an existing DB")
    
        #________________Calculate foldseek self comparison of all predicted structures____________
    
        for file in pdb_files:
            foldseek_run = ["foldseek", "easy-search", file, db_directory + "DB", file.replace(".pdb","-self.foldseek"), blind_path + "/tmp", "--format-mode", "0", "--format-output", "query,target,alntmscore,qaln,taln,alnlen,evalue,bits", "--exhaustive-search", "1", "-s", "9.5"]
            if not os.path.isfile(file.replace(".pdb","-self.foldseek")):
                response = subprocess.run(foldseek_run, capture_output=True, text=True, check=True)
                try:
                    response = subprocess.run(foldseek_run, capture_output=True, text=True, check=True)
                    print(response.check_returncode())
                except subprocess.CalledProcessError as e:
                    print("foldseek failed to run {:}".format(file))
                    print("Error:", e.stderr)
                print('{:} succeeded!!!'.format(file))
            else:
                print("{:} already exists".format(file.replace(".pdb","-self.foldseek")))
    
        #__________________________________________________________________________________________
    
    
        #__________Populate a correlation matrix with bit scores_______________________________________________
    
        #everything will be sorted by the text of the file name
        files = glob.glob(blind_path + "/**/*-self.foldseek")
    
        #first remove any outliers from the dssp loop distribution, they tend to be unfolded predictions
        files_dssp = [];files_count = [];
        for file in files:
            u = mda.Universe(file.replace("-self.foldseek",".pdb"))
            s = DSSP(u).run().results.dssp[0]
            dssp, count = np.unique(s, return_counts=True)
            # ['-' 'E' 'H']
            if len(dssp) < 3:
                if '-' not in dssp:
                    dssp  = np.insert(dssp, 0 ,'-')
                    count = np.insert(count,0, 0)
                if 'E' not in dssp:
                    dssp  = np.insert(dssp, 1 ,'E')
                    count = np.insert(count,1, 0)
                if 'H' not in dssp:
                    dssp  = np.insert(dssp, 2 ,'H')
                    count = np.insert(count,2, 0)
            files_dssp.append(dssp)
            files_count.append(count)
        files_dssp = np.array(files_dssp); files_count = np.array(files_count);
        z_scores = stats.zscore(files_count[:, 0])
        outlier_idx = np.argwhere(z_scores > 3)
    
        # remove unfolded proteins from file list 
        files = np.array(files)
        mask = np.zeros(files.shape, dtype=bool)
        mask[outlier_idx] = True
        for file in files[mask]:
            print("removed from analysis: ",file.replace("-self.foldseek",".pdb"))
        files = files[~mask]
        files = sorted(files)
        files_pdb = [file.replace('/','-')[17:].replace("-self.foldseek","") for file in files]
        #files_pdb = [file.replace("-self.foldseek",".pdb") for file in files]
        corr_mtx = []
    
        df = {}
        for file in files:
            # it is possible for predictions to be so different that it isn't returned with a bit_score
            # in that case we return a zero
            dict_with_all = {file:[0] for file in files_pdb}
            with open(file, 'r') as _:
                data = [l.rstrip().split('\t') for l in _]
            for d in data:
                dict_with_all[d[1]] = d
                print(dict_with_all[d[1]])
            #bug in foldseek occasionally returns -2,147,483,648
            _temp = []
            for pdb in files_pdb:
                print("testing", pdb)
                x = int(dict_with_all.get(pdb, 0)[-1])
                print(x)
                if x == -2147483648:
                    _temp.append(0)
                else:
                    _temp.append(x)
    
            corr_mtx.append(_temp)
        
        corr_mtx = np.array(corr_mtx)
    
        #normalize each row and subtract top model from full MSA depth to give more
        #specific meaning to variance
        norm_corr_mtx = minmax_scale(corr_mtx, axis=1)
        norm_corr_mtx = (norm_corr_mtx + norm_corr_mtx.T) /2
    
        sklearn_pca = PCA(n_components=4)
        pca = sklearn_pca.fit_transform(norm_corr_mtx)
        labels = blind_screening.cluster_structures(pca)
        
        plt.figure(figsize=(8,6))
        plt.scatter(pca[:,0], pca[:,1], c=labels, cmap='viridis', s=45)
        plt.savefig(blind_path + '/' + pdb1_name + '-cluster.png')
        plt.clf()
    
    
        #find the structures_of_interest
        files_of_interest = []
        pca_of_interest = []
        for l in np.unique(labels):
            kmed_idx, tot_cost  = blind_screening.k_medoids(pca, l, labels)
            for idx in kmed_idx:
                files_of_interest.append([files[idx], l])
                pca_of_interest.append(pca[idx])
    
        #create pse file with colors that match viridis colors in cluster.png
        viridis = plt.get_cmap('viridis',len(files_of_interest))
        largest_group_num = max(files_of_interest, key=lambda x: x[1])
        pymol.cmd.load(files[0].replace('-self.foldseek','.pdb'), 'Dominant')
        with open(pdb1_name + "-structures_of_interest.csv", "w") as file:
            file.write("group, file, pca_1, pca_2\n")
    
        with open(blind_path + '/' + pdb1_name + "-structures_of_interest.csv", "a") as file:
            for idx, foi in enumerate(files_of_interest):
                if largest_group_num[1] == -1:
                    color = 0
                else:
                    color = (foi[1] + 1) / (largest_group_num[1]+1)
                color = viridis(color)[:3]
                new_name = re.findall(r'(full)|(max\w+)|(rank_\d+)', foi[0])
                new_name = str(idx)+ '_' + '_'.join([i for n in new_name for i in n if i != ''])
                pymol.cmd.load(foi[0].replace('-self.foldseek','.pdb'), new_name)
                pymol.cmd.align(new_name,'Dominant')
                color_name = 'col_'+str(foi[1])
                pymol.cmd.set_color(color_name, color)
                pymol.cmd.color(color_name,new_name)
                file.write(f"{foi[1]}, {foi[0]}, {pca_of_interest[idx][0]}, {pca_of_interest[idx][1]}\n")
    
        pymol.cmd.save(blind_path + '/' + pdb1_name + '-structures_of_interest.pse', 'pse')
        pymol.cmd.delete('all')
        pymol.cmd.reinitialize()
    
        #save all data with clusters
        with open("structures_all.csv", 'w') as file:
            file.write("group, file, pca_1, pca_2\n")
            for idx, f in enumerate(files):
                file.write(f"{labels[idx]},{f},{pca[idx, 0]},{pca[idx, 1]}\n")
    
        sys.exit()
    
    
