import os
import numpy as np

def conf_prep(mu,beta,D,w):
    file = open("conf.txt","w")
    file.writelines("w = %d\nb = %d\nD = %d\nmu = %f" % (w,beta,D,mu))
    file.close()

# Update this with the path to msgsteiner
msgpath = "/home/a/Downloads/msgsteiner-1.3/msgsteiner"

# The interval stop value is not included in the range
mu_range = np.arange(2.1,2.4,0.3)
beta_range = np.arange(150,160,10)
# Increase the stop value to test two values of w
w_range = np.arange(2,3.1,1)
prize_file = "gbm_prize.txt"
edge_file = "/home/a/Downloads/OmicsIntegrator-master/data/iref_mitab_miscore_2013_08_12_interactome.txt"
conf_file = "conf.txt"
wt_path = "/home/a/Downloads/OmicsIntegrator-master/results/WT/"
ko_path = "/home/a/Downloads/OmicsIntegrator-master/results/KO/"

D =3

# Create output directories if needed
if not os.path.exists(wt_path):
    os.makedirs(wt_path)
if not os.path.exists(ko_path):
    os.makedirs(ko_path)

for mu in mu_range:
    for beta in beta_range:
        for w in w_range:
            conf_prep(mu,beta,D,w)
            out_label = "WT_w%f_beta%d_D%d_mu%f" %(w,beta,D,mu)
            os.system("python /home/a/Downloads/OmicsIntegrator-master/scripts/forest.py --prize %s --edge %s --conf conf.txt --msgpath %s --outpath %s --outlabel %s" %(prize_file,edge_file,msgpath,wt_path,out_label))
            out_label = "KO_w%f_beta%d_D%d_mu%f" %(w,beta,D,mu)
            os.system("python /home/a/Downloads/OmicsIntegrator-master/scripts/forest.py --prize %s --edge %s --conf conf.txt --msgpath %s --knockout EGFR --outpath %s --outlabel %s" %(prize_file,edge_file,msgpath,ko_path,out_label))
