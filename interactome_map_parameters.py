'''
8.3.2020
Ali Yurtseven

Interactome_map, Cytoscape output finder

Here, both SAINTExspress_v3_6.3, OmicsIntegrator1 was used.

'''
import pandas as pd
import numpy as np
import shutil
import os
from os.path import expanduser
from datetime import date
from math import log
from gprofiler import GProfiler
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as  mpl
import pandas
import pathlib
# Parameters needed to be entered carefully
path_parameter_file=os.path.join(pathlib.Path(__file__).parent.absolute(),'parameters.txt')

if os.path.exists(path_parameter_file)==False:
    print('parameters.txt file should be in the same directory with the interactome_parameter.py file! There will be an error!')

data=pd.read_csv(path_parameter_file,delimiter='\t',header=None)

# Let me take the parameters!

output_path= data[data.columns[1]][0]
output_name= data[data.columns[1]][1]
path= data[data.columns[1]][2]
path2= data[data.columns[1]][3]
bait_name= data[data.columns[1]][4]
ask_user= data[data.columns[1]][5]
number_of_conditions= int(data[data.columns[1]][6])
gene_names_path= data[data.columns[1]][7]
gene_names_path2= data[data.columns[1]][9]
label1= data[data.columns[1]][11]
label2= data[data.columns[1]][12]


# Let me create the output folder path using os system

all_output_path=os.path.join(output_path,output_name)

if os.path.exists(all_output_path):
    existence=input("You have the same folder in the given directory! Do you want to remove it? Type:yes if so, Type:no if you don't want to: ")
    if existence=="yes":
        shutil.rmtree(all_output_path)
        os.makedirs(all_output_path)
        os.makedirs(os.path.join(all_output_path,'SaintExpress'))
    else:
        output_name=input('Give me another name for output folder!')
        all_output_path = os.path.join(output_path, output_name)
        os.makedirs(all_output_path)
        os.makedirs(os.path.join(all_output_path,'SaintExpress'))
else:
    os.makedirs(all_output_path)
    os.makedirs(os.path.join(all_output_path, 'SaintExpress'))

'''
with a specific feature was converted into  a list!
'''

def reader(path,path2):
    df2 = pd.read_excel(path)
    df2_c = pd.read_excel(path2)
    df=df2.dropna()
    df_c=df2_c.dropna()
    Peptide_number=df['# AAs'].tolist()
    Peptide_number_c=df_c['# AAs'].tolist()
    PSM_val=df['Σ# PSMs'].tolist()
    PSM_val_c=df_c['Σ# PSMs'].tolist()
    Accession_name=df['Accession'].tolist()
    Accession_name_c=df_c['Accession'].tolist()
    Gene_all=df['Description'].tolist()
    Gene_all_c=df_c['Description'].tolist()
    return Gene_all,Gene_all_c, Accession_name,Accession_name_c,Peptide_number,Peptide_number_c,PSM_val,PSM_val_c
Gene_all,Gene_all_c,Accession_name,Accession_name_c,Peptide_number,Peptide_number_c,PSM_val,PSM_val_c= reader(path,path2)

'''
After having our features in a list, we need to get the correct name of the gene family/ Please get the one named
after GN=, don't get the the last names! As a result a function for the filtration and adjustment of the features were written 

Also, in the excell file, since there are 2 replicates, the name of  the features were written twice. Therefore, I removed this
feature names interfere with the name of my feature values. 

Here is the function that I wrote!!
'''


'''
Since there are some Uncharachterized proteins in the list, there is no GN name for them, so I wrote Uncharachterized
protein for that gene names. First, I found the indexes if these proteins. then, I put name of UNcharachterized protein to
these indexes and I extracted Description name from this list!!
'''



def filter(Gene_all,Gene_all_c,Peptide_number,Accession_name,Accession_name_c,Peptide_number_c,PSM_val,PSM_val_c):
    GeneID_list= []
    index=[]
    ind=0
    for i in Gene_all:
        counter=0
        t=[]
        z=[]
        for j in i:
            t.append(j)
            if t[-1]=='=' and t[-2]=='N' and t[-3]=='G':
                index.append(ind)
                gene_name=i[counter+1:-1]
                for k in gene_name:
                   if k!=' ':
                       z.append(k)
                   else:
                       break
                y=''.join(z)
                GeneID_list.append(y)
            counter=counter+1
        ind=ind+1
    # Let me find the indexes of the proteins that have the uncharacterized proteins!
    c2=[]
    for i in range(0,len(Gene_all)-1):
        if i not in index:
            c2.append(i)
# 'c2' is a list of the indexes, now, let's form our gene list again
    one_in_desc2=[]
    for i in c2:
        if Gene_all[i]=='Description':
            one_in_desc2.append(i)

    count=0

    for i in c2:
        GeneID_list.insert(i,'Uncharachterized_Protein')
    gene=GeneID_list



    Accession_ID = []
    for i in Accession_name:
        if i != 'Accession':
            Accession_ID.append(i)
    Peptide_num=[]
    for i in Peptide_number:
        if i!='# AAs':
            Peptide_num.append(i)

    PSM_filtered=[]
    for i in PSM_val:
        if i!='Σ# PSMs':
            PSM_filtered.append(i)

    '''
    Since there are some Uncharachterized proteins in the list, there is no GN name for them, so I wrote Uncharachterized
    protein for that gene names. First, I found the indexes if these proteins. then, I put name of UNcharachterized protein to
    these indexes and I extracted Description name from this list!!
    '''
    GeneID_list_c = []
    index=[]
    ind=0
    for i in Gene_all_c:
        counter=0
        t=[]
        z=[]
        for j in i:
            t.append(j)
            if t[-1]=='=' and t[-2]=='N' and t[-3]=='G':
                index.append(ind)
                gene_name=i[counter+1:-1]
                for k in gene_name:
                   if k!=' ':
                       z.append(k)
                   else:
                       break
                y=''.join(z)
                GeneID_list_c.append(y)
            counter=counter+1
        ind=ind+1
    # Let me find the indexes of the proteins that have the uncharacterized proteins!
    c=[]
    for i in range(0,len(Gene_all_c)-1):
        if i not in index:
            c.append(i)
    # 'c' is a list of the indexes, now, let's form our gene list again
    one_in_desc=[]
    for i in c:
        if Gene_all_c[i]=='Description':
            one_in_desc.append(i)

    count=0
    for i in c:
        GeneID_list_c.insert(i,'Uncharachterized_Protein')
    gene_c=GeneID_list_c
    '''
    for i in c:
    k=Gene_all_c[i]
    print(k)
    '''

    Accession_ID_c = []
    for i in Accession_name_c:
        if i != 'Accession':
            Accession_ID_c.append(i)

    Peptide_num_c = []
    for i in Peptide_number_c:
        if i != '# AAs':
            Peptide_num_c.append(i)

    PSM_filtered_c=[]
    for i in PSM_val_c:
        if i!='Σ# PSMs':
            PSM_filtered_c.append(i)

    return gene,gene_c,Peptide_num,Peptide_num_c,Accession_ID,Accession_ID_c,PSM_filtered,PSM_filtered_c

gene,gene_c,Peptide_num,Peptide_num_c,Accession_ID,Accession_ID_c,PSM_filtered,PSM_filtered_c=filter(Gene_all,Gene_all_c,Peptide_number,Accession_name,Accession_name_c,Peptide_number_c,PSM_val,PSM_val_c)

# I want to look at if there are any Accession IDs different in control!

def index_finder_for_AccessionID(Accession_ID,Accession_ID_c):
    Unique_Acccession_controls=[]
    ind_for_Accession_c=[]
    count=0
    for i in Accession_ID_c:
        if i not in Accession_ID:
            if i not in Unique_Acccession_controls:
                Unique_Acccession_controls.append(i)
                ind_for_Accession_c.append(count)
        count=count+1

    Unique_Acccession=[]
    ind_for_Accession=[]
    count2=0
    for i in Accession_ID:
        if i not in Accession_ID_c:
            if i not in Unique_Acccession:
                Unique_Acccession.append(i)
                ind_for_Accession.append(count2)
        count2=count2+1


    Mutual_Accession=[]
    ind_mutual=[]
    count3=0
    for i in Accession_ID:
        if i in Accession_ID_c:
            if i not in Mutual_Accession:
                Mutual_Accession.append(i)
                ind_mutual.append(count3)
        count3=count3+1

    return ind_for_Accession,ind_for_Accession_c,ind_mutual,Unique_Acccession_controls,Unique_Acccession,Mutual_Accession
ind_for_Accession,ind_for_Accession_c,ind_mutual,Unique_Acccession_controls,Unique_Acccession,Mutual_Accession=index_finder_for_AccessionID(Accession_ID,Accession_ID_c)
# Let me find the lists for each exclusive features wanted for prey file. Left one is Accession ID, middle is AA. number, and right one is Gene name.

# Let's first find the Accession IDs
Accession=[]
for i in Unique_Acccession:
    Accession.append(i)
for j in Unique_Acccession_controls:
    Accession.append(j)
for k in Mutual_Accession:
    Accession.append(k)

#Let's find petide numbers

AA_length=[]
for i in ind_for_Accession:
    each=Peptide_num[i]
    AA_length.append(each)
for j in  ind_for_Accession_c:
    each2=Peptide_num_c[j]
    AA_length.append(each2)
for k in ind_mutual:
    each3=Peptide_num[k]
    AA_length.append(each3)

Gene=[]
for i in ind_for_Accession:
    each=gene[i]
    Gene.append(each)
for j in  ind_for_Accession_c:
    each2=gene_c[j]
    Gene.append(each2)
for k in ind_mutual:
    each3=gene[k]
    Gene.append(each3)

# Let me write the bait file
def bait_writer():
    with open(os.path.join(os.path.join(all_output_path,'SaintExpress'),'bait.txt'),'w') as wr:
        wr.write('{}-1'.format(bait_name)+'\t'+'{}'.format(bait_name)+'\t'+'T'+'\n')
        wr.write('{}-2'.format(bait_name)+'\t'+'{}'.format(bait_name)+'\t'+'T'+'\n')
        wr.write('Mark_1'+'\t'+'Mark_1'+'\t'+'C'+'\n')
        wr.write('Mark_2'+'\t'+'Mark_2'+'\t'+'C')
    return
bait_writer()
# Now, let me write pray file

def pray_file_writer(Gene,AA_lenght,Accession):
    with open(os.path.join(os.path.join(all_output_path,'SaintExpress'),'pray.txt'),'w') as pdd:
        count=0
        for i in Gene:
            k='{}\t{}\t{}\n'.format(Accession[count],AA_lenght[count],i)
            pdd.write(k)
            count=count+1
    return
pray_file_writer(Gene,AA_length,Accession)

# Let me find interaction file, for this we need PSM_values, and Gene name, so called, pray name.

df=pd.read_excel(path)
Accession_name=df['Accession'].tolist()
Accession_name_c=pd.read_excel(path2)

count=0
for i in Accession_name:
    if i=='Accession':
        number=count
    count=count+1

count2=0
for j in Accession_name_c:
    if j=='Accession':
        number_c=count2
    count2=count2+1
print(number)

#Since I found each group in both condition case and control case, let me find the interaction file now.

def inter_file(number,number_c,Accession_ID,Accession_ID_c,PSM_filtered,PSM_filtered_c):
    Accession_sample1=[]
    Accession_sample2=[]
    PSM_value_sample1=[]
    PSM_value_sample2=[]

    for i in range(0,number):
        Accession_sample1.append(Accession_ID[i])
        PSM_value_sample1.append(PSM_filtered[i])
    for j in range(number,len(PSM_filtered)-1):
        Accession_sample2.append(Accession_ID[j])
        PSM_value_sample2.append(PSM_filtered[j])
    Accession_control1=[]
    Accession_control2=[]
    PSM_value_control1=[]
    PSM_value_control2=[]

    for i in range(0, number_c):
        Accession_control1.append(Accession_ID_c[i])
        PSM_value_control1.append(PSM_filtered_c[i])
    for j in range(number_c, len(Accession_ID_c)):
        Accession_control2.append(Accession_ID_c[j])
        PSM_value_control2.append(PSM_filtered_c[j])
    return Accession_control1,Accession_control2,PSM_value_control1,PSM_value_control2,Accession_sample1,Accession_sample2,PSM_value_sample1,PSM_value_sample2
Accession_control1,Accession_control2,PSM_value_control1,PSM_value_control2,Accession_sample1,Accession_sample2,PSM_value_sample1,PSM_value_sample2=inter_file(number,number_c,Accession_ID,Accession_ID_c,PSM_filtered,PSM_filtered_c)

def inter(Accession_control1,Accession_control2,PSM_value_control1,PSM_value_control2,Accession_sample1,Accession_sample2,PSM_value_sample1,PSM_value_sample2):
    with open(os.path.join(os.path.join(all_output_path,'SaintExpress'),'inter.txt'),'w') as pdd:
        count=0
        for i in (Accession_sample1):
            m='{}-1\t{}\t{}\t{}\n'.format(bait_name,bait_name,i,PSM_value_sample1[count])
            pdd.write(m)
            count=count+1
        count3=0
        for j in (Accession_sample2):
            m='{}-2\t{}\t{}\t{}\n'.format(bait_name,bait_name,j,PSM_value_sample2[count3])
            pdd.write(m)
            count3=count3+1
        count4=0
        for i in (Accession_control1):
            m='Mark_1\tMark_1\t{}\t{}\n'.format(i,PSM_value_control1[count4])
            pdd.write(m)
            count4=count4+1
        count5=0
        for i in (Accession_control2):
            m='Mark_2\tMark_2\t{}\t{}\n'.format(i,PSM_value_control2[count5])
            pdd.write(m)
            count5=count5+1
    return
inter(Accession_control1,Accession_control2,PSM_value_control1,PSM_value_control2,Accession_sample1,Accession_sample2,PSM_value_sample1,PSM_value_sample2)
# Now, you need to run Saint Express! for that, you need to move these 3 files into home directory

shutil.copy(os.path.join(os.path.join(all_output_path,'SaintExpress'),'inter.txt'),expanduser("~"))
shutil.copy(os.path.join(os.path.join(all_output_path,'SaintExpress'),'bait.txt'),expanduser("~"))
shutil.copy(os.path.join(os.path.join(all_output_path,'SaintExpress'),'pray.txt'),expanduser("~"))

# if you add your path into your system, you need to run saint express!
os.chdir(expanduser("~"))

path= data[data.columns[1]][13]
technical_replica_number=data[data.columns[1]][14]
os.system("{} {} inter.txt pray.txt bait.txt".format(path,technical_replica_number))

# Now, we created result of the SAINTexpress, in the home path. I will move this to my SAINTexpress file

shutil.move(os.path.join(expanduser("~"),"list.txt"),os.path.join(all_output_path,'SaintExpress'))

# I would like to change the name of my_file by using the bait_name as file name and saint express.
os.chdir(os.path.join(all_output_path,'SaintExpress'))

time= date.today()
os.rename('list.txt','{}_Saintresult_{}.txt'.format(bait_name,time))


# Now, let me show saint scores that are significantly changed in these 2 conditions and convert data to excell file!

df=pd.read_csv(os.path.join(os.path.join(all_output_path,'SaintExpress'),'{}_Saintresult_{}.txt'.format(bait_name,time)),delimiter='\t')

df.loc[(df['SaintScore']==1) &(df['BFDR']<0.05),'Significant Ones (+), NonSignificant Ones(-)']='+'
df.loc[df['Significant Ones (+), NonSignificant Ones(-)'].isnull(),'Significant Ones (+), NonSignificant Ones(-)']='-'

df= df.sort_values(by='Significant Ones (+), NonSignificant Ones(-)')

df=df.reset_index()
df=df.drop(['index'],axis=1)

df.to_excel(os.path.join(os.path.join(all_output_path,'SaintExpress'),'{}_Saintresult_{}.xlsx'.format(bait_name,time)), index=None)

# Now, I need to run omics integrator. For this, firstly, I need to create a directory for both input and output of omics integrator

saint_text_path=os.path.join(os.path.join(all_output_path,'SaintExpress'),'{}_Saintresult_{}.txt'.format(bait_name,time))

# For the further analysis, Saint score was taken as more than 0.5 and BFDR is less than or equal to  0.05
# Let me filter the data using the filter function, then, this data was demonstrated in an excell file and the further analysis was executed.

def filter(saint_text_path):
    df = pd.read_table(saint_text_path)
    df2=df[(df['SaintScore']>0.5)&(df['BFDR']<=0.05)]
    return df2

df2=filter(saint_text_path)
df2= df2.reset_index()
df2.to_excel(os.path.join(os.path.join(all_output_path,'SaintExpress'),'{}_Saintresult_only_saint_proteins{}.xlsx'.format(bait_name,time)))

# This part was written for both GO anotation and Omics Integrator analysis!

def file_formater(df2):
    gene_names=df2['PreyGene'].tolist()
    fold_change=df2['FoldChange'].tolist()
    saint_score=df2['SaintScore'].tolist()
    log2_FC=[]
    for i in fold_change:
        s= log(i,2)
        log2_FC.append(s)
    return gene_names,fold_change,saint_score,log2_FC
gene_names,fold_change,saint_score,log2_FC=file_formater(df2)

def gene_name_finder_for_GO_Annotation(gene_names): # This is written for having a text file to copy the genes and do the GO anotations in g-profiler!
    with open(os.path.join(os.path.join(all_output_path,'SaintExpress'),'{}_Saintresult_GO_gene_names{}.txt'.format(bait_name,time)),'w') as pddd:
        for i in gene_names:
            pddd.write('{}\n'.format(i))
    return
gene_name_finder_for_GO_Annotation(gene_names)


os.mkdir(os.path.join(all_output_path,'Omics Integrator'))


def prize_file_writer(gene_names,log2_FC):
    with open(os.path.join(os.path.join(all_output_path,'Omics Integrator'),'gbm_prize.txt'),'w') as pdd:
        count=0
        for i in gene_names:
            pdd.write('{}\t'.format(i))
            pdd.write('{}\n'.format(log2_FC[count]))
            count=count+1
    return

prize_file_writer(gene_names,log2_FC)

# Let me now move this file into .../OmicsIntegrator-master/example/GBM/ directory and run the python code in this diretory.

omics_file=data[data.columns[1]][15]
# Let me copy gbm_prize file into the directory of interest!
shutil.copy(os.path.join(os.path.join(all_output_path,'Omics Integrator'),'gbm_prize.txt'),os.path.join(omics_file,'example/GBM'))

# Let me run the python code given as GBM_case_study

'''
Output this run will be in KO and WT files. You should not  change the directory since this program takes 
the files into your result file in Omics Integrator
'''
removed_files=os.listdir(os.path.join(omics_file,'results/KO'))
for f in removed_files:
    os.remove(os.path.join(os.path.join(omics_file,'results/KO'),f))

os.chdir(os.path.join(omics_file,'example/GBM'))
os.system("python {}".format(os.path.join(omics_file,'example/GBM/GBM_case_study.py')))

# Let me copy the output into your result file, before open result file in Omics Integrator

os.makedirs(os.path.join(os.path.join(all_output_path,'Omics Integrator'),'output'))

source=os.path.join(omics_file,'results/KO')
destination=os.path.join(os.path.join(all_output_path,'Omics Integrator'),'output')

files=os.listdir(source)
for f in files:
    shutil.move(os.path.join(source,f), destination)

# Let me now find the GO anaotation graphs  for the proteins that have SaintExpress score>0.5 and BFDR<0.01

'''
This code asks you if you have one or more than one conditions. In case, you have one, it gives you only one conditional 
horizontal bar graphs. If not, it would compare the bar graphs.  

To have 2 conditions, you need to run this code 2 times with different outputs! 
'''

if ask_user=="YES":
    if number_of_conditions==1:
        # Only one condition! getting GO annotation profiles of proteins that have >0.5 saint score and <0.01 BFDR score

        gp = GProfiler(return_dataframe=True)

        profiler = gp.profile(organism='hsapiens', query=gene_names)

        BP_profiler = profiler[profiler["source"] == "GO:BP"]
        CC_profiler = profiler[profiler["source"] == "GO:CC"]
        MF_profiler = profiler[profiler["source"] == "GO:MF"]

        BP_profiled = BP_profiler.sort_values(by=["p_value"])
        CC_profiled = CC_profiler.sort_values(by=["p_value"])
        MF_profiled = MF_profiler.sort_values(by=["p_value"])

        location_BP = BP_profiled["name"].to_list()[0:10]
        p_BP = BP_profiled["p_value"].to_list()[0:10]
        logged_p_BP = []
        for i in p_BP:
            x = -log(i, 10)
            logged_p_BP.append(x)

        location_CC = CC_profiled["name"].to_list()[0:10]
        p_CC = CC_profiled["p_value"].to_list()[0:10]
        logged_p_CC = []
        for i in p_CC:
            y = -log(i, 10)
            logged_p_CC.append(y)

        location_MF = MF_profiled["name"].to_list()[0:10]
        p_MF = MF_profiled["p_value"].to_list()[0:10]
        logged_p_MF = []
        for i in p_MF:
            z = -log(i, 10)
            logged_p_MF.append(z)

        # Drawing the GO:MF annotation graph

        label_size = 8
        mpl.rcParams['ytick.labelsize'] = label_size

        df = pandas.DataFrame(dict(graph=location_MF,
                                   n=logged_p_MF))
        ind = np.arange(len(df))
        width = 0.35

        fig, ax = plt.subplots()

        ax.barh(ind, df.n, width, color='green')

        hfont = {'fontname': 'Arial'}

        plt.ylabel('ylabel', **hfont)
        ax.set(yticks=ind + 1 * width / 8, yticklabels=df.graph, ylim=[2 * width - 1, len(df)])
        plt.xlim(0, int(logged_p_MF[0]) + 4, 2)
        ax.set_title('GO: Molecular Function')
        plt.xlabel('-log10(Adjusted p Value) ')

        from math import log

        plt.axvline(-log(0.05, 10), color='y', linestyle='dashed', linewidth=1)
        plt.show()

        # Drawing the GO:CC annotation graph

        label_size = 8
        mpl.rcParams['ytick.labelsize'] = label_size

        df = pandas.DataFrame(dict(graph=location_CC,
                                   n=logged_p_CC))
        ind = np.arange(len(df))
        width = 0.35

        fig, ax = plt.subplots()

        ax.barh(ind, df.n, width, color='green')

        hfont = {'fontname': 'Arial'}

        plt.ylabel('ylabel', **hfont)
        ax.set(yticks=ind + 1 * width / 8, yticklabels=df.graph, ylim=[2 * width - 1, len(df)])
        plt.xlim(0, int(logged_p_CC[0]) + 4, 2)
        ax.set_title('GO: Cellular Component')
        plt.xlabel('-log10(Adjusted p Value) ')

        from math import log

        plt.axvline(-log(0.05, 10), color='y', linestyle='dashed', linewidth=1)
        plt.show()

        # Drawing the GO:BP annotation graph

        label_size = 8
        mpl.rcParams['ytick.labelsize'] = label_size

        df = pandas.DataFrame(dict(graph=location_BP,
                                   n=logged_p_BP))
        ind = np.arange(len(df))
        width = 0.35

        fig, ax = plt.subplots()

        ax.barh(ind, df.n, width, color='green')

        hfont = {'fontname': 'Arial'}

        plt.ylabel('ylabel', **hfont)
        ax.set(yticks=ind + 1 * width / 8, yticklabels=df.graph, ylim=[2 * width - 1, len(df)])
        plt.xlim(0, int(logged_p_BP[0]) + 4, 2)
        ax.set_title('GO: Biological Process')
        plt.xlabel('-log10(Adjusted p Value) ')

        from math import log

        plt.axvline(-log(0.05, 10), color='y', linestyle='dashed', linewidth=1)
        plt.show()
    elif number_of_conditions==2:

        first_sample_data = pd.read_excel(gene_names1_path)
        second_sample_data = pd.read_excel(gene_names2_path)

        # First, we need to take BP, MF, and CC of these. Then, we can merge these features!

        # g-profiler for the first sample

        gene_names1 = first_sample_data["PreyGene"].to_list()

        gp = Gprofiler(return_dataframe=True)

        profiler1 = gp.profile(organism='hsapiens', query=gene_names1)

        # Now, lets do this for the second sample!

        gene_names2 = second_sample_data["PreyGene"].to_list()

        gp = Gprofiler(return_dataframe=True)

        profiler2 = gp.profile(organism='hsapiens', query=gene_names2)

        # Let me get BP,CC,and MF for the first sample and sort them according to their p values

        BP_profiler1 = profiler1[profiler1["source"] == "GO:BP"]
        CC_profiler1 = profiler1[profiler1["source"] == "GO:CC"]
        MF_profiler1 = profiler1[profiler1["source"] == "GO:MF"]

        BP_profiled1 = BP_profiler1.sort_values(by=["p_value"])
        CC_profiled1 = CC_profiler1.sort_values(by=["p_value"])
        MF_profiled1 = MF_profiler1.sort_values(by=["p_value"])

        # Let me get BP,CC,MF for the second sample and sort them according to their p values

        BP_profiler2 = profiler2[profiler2["source"] == "GO:BP"]
        CC_profiler2 = profiler2[profiler2["source"] == "GO:CC"]
        MF_profiler2 = profiler2[profiler2["source"] == "GO:MF"]

        BP_profiled2 = BP_profiler2.sort_values(by=["p_value"])
        CC_profiled2 = CC_profiler2.sort_values(by=["p_value"])
        MF_profiled2 = MF_profiler2.sort_values(by=["p_value"])

        # Now, let me merge BPs, CCs, and MFs

        merged_BP = BP_profiled1[['name', 'p_value']].merge(BP_profiled2[['name', 'p_value']], how='inner',
                                                            left_on='name', right_on='name')
        merged_CC = CC_profiled1[['name', 'p_value']].merge(CC_profiled2[['name', 'p_value']], how='inner',
                                                            left_on='name', right_on='name')
        merged_MF = MF_profiled1[['name', 'p_value']].merge(MF_profiled2[['name', 'p_value']], how='inner',
                                                            left_on='name', right_on='name')


        # Let me take name of the intersected groups for each annotation and take -log10 of p values for these groups!
        location_merged_BP = merged_BP['name'].to_list()[0:10]
        location_merged_CC = merged_CC['name'].to_list()[0:10]
        location_merged_MF = merged_MF['name'].to_list()[0:10]

        p_BP_x = merged_BP['p_value_x'].tolist()[0:10]
        logged_bp_x = -np.log10(p_BP_x)

        p_BP_y = merged_BP['p_value_y'].tolist()[0:10]
        logged_bp_y = -np.log10(p_BP_y)

        p_CC_x = merged_CC['p_value_x'].tolist()[0:10]
        logged_cc_x = -np.log10(p_CC_x)

        p_CC_y = merged_CC['p_value_y'].tolist()[0:10]
        logged_cc_y = -np.log10(p_CC_y)

        p_MF_x = merged_MF['p_value_x'].tolist()[0:10]
        logged_mf_x = -np.log10(p_MF_x)

        p_MF_y = merged_MF['p_value_y'].tolist()[0:10]
        logged_mf_y = -np.log10(p_MF_y)

        # Let me draw cellular component graph
        label_size = 8
        mpl.rcParams['ytick.labelsize'] = label_size

        df = pandas.DataFrame(dict(graph=location_merged_CC,
                                   n=logged_cc_x, m=logged_cc_y))
        ind = np.arange(len(df))
        width = 0.45

        fig, ax = plt.subplots()

        ax.barh(ind, df.n, width, color='red', label='{}'.format(label1))
        ax.barh(ind + width, df.m, width, color='green', label='{}'.format(label2))
        hfont = {'fontname': 'Arial'}

        plt.ylabel('ylabel', **hfont)
        ax.set(yticks=ind + 1 * width / 2, yticklabels=df.graph, ylim=[2 * width - 1, len(df)])
        plt.xlim(0, max(logged_cc_x[0], logged_cc_y[0]) + 4, 2)
        ax.set_title('GO: Cellular Component')
        plt.xlabel('-log10(Adjusted p Value) ')

        ax.legend()
        from math import log

        plt.axvline(-log(0.05, 10), color='y', linestyle='dashed', linewidth=1)
        plt.show()

        # Let me draw biological process graph

        label_size = 8
        mpl.rcParams['ytick.labelsize'] = label_size

        df = pandas.DataFrame(dict(graph=location_merged_BP,
                                   n=logged_bp_x, m=logged_bp_y))
        ind = np.arange(len(df))
        width = 0.45

        fig, ax = plt.subplots()

        ax.barh(ind, df.n, width, color='red', label='{}'.format(label1))
        ax.barh(ind + width, df.m, width, color='green', label='{}'.format(label2))
        hfont = {'fontname': 'Arial'}

        plt.ylabel('ylabel', **hfont)
        ax.set(yticks=ind + 1 * width / 2, yticklabels=df.graph, ylim=[2 * width - 1, len(df)])
        plt.xlim(0, max(logged_bp_x[0], logged_bp_y[0]) + 4, 2)
        ax.set_title('GO: Biological Process')
        plt.xlabel('-log10(Adjusted p Value) ')

        ax.legend()
        from math import log

        plt.axvline(-log(0.05, 10), color='y', linestyle='dashed', linewidth=1)
        plt.show()

        # Let me draw molecular function graph
        label_size = 8
        mpl.rcParams['ytick.labelsize'] = label_size

        df = pandas.DataFrame(dict(graph=location_merged_MF,
                                   n=logged_mf_x, m=logged_mf_y))
        ind = np.arange(len(df))
        width = 0.45

        fig, ax = plt.subplots()

        ax.barh(ind, df.n, width, color='red', label='{}'.format(label1))
        ax.barh(ind + width, df.m, width, color='green', label='{}'.format(label2))
        hfont = {'fontname': 'Arial'}

        plt.ylabel('ylabel', **hfont)
        ax.set(yticks=ind + 1 * width / 2, yticklabels=df.graph, ylim=[2 * width - 1, len(df)])
        plt.xlim(0, max(logged_mf_x[0], logged_mf_y[0]) + 4, 2)
        ax.set_title('GO: Molecular Function')
        plt.xlabel('-log10(Adjusted p Value) ')

        ax.legend()
        from math import log

        plt.axvline(-log(0.05, 10), color='y', linestyle='dashed', linewidth=1)
        plt.show()
