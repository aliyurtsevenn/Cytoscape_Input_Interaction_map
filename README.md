# Cytoscape_input_generator


This program gets AP-MS interaction data as input and using both SAINTexpress, Omics Integrator and gprofiler, it both annotates significantly interacting proteins and generate cytoscape input files. Since SAINTExpress and Omics Integrator works only in Linux operating system, unfortunately, the code also works only in Linux operating system. You shouldn't try to run it in Windows or Mac operating system. 

3 programs are required for getting both interactome map and GO annotation graphs. These are SAINTexpress, OmicsIntegrator, and Cytoscape.

1. Downloading SAINTexpress

SAINTexpress helps you to find out which molecules are significantly interacting with each other in your interactome data. It gives you output in text format. Here, my code convert this text file into excell file and show the ones that have saint score values bigger than 0.5 and false discovery rate less than 0.01. 

For GO annotation analysis, obtained excell file in SAINTexpress folder is used as input for gprofiler package.  

This program is working in 3.6.3 version of SAINTExpress. You can download it here: https://sourceforge.net/projects/saint-apms/files/

After downloading it, you need to  install g++ or gcc to your computer, use the command

>sudo apt-get install build-essential

Then, you need to compile your SAINTexpress_v3.6.3__2018-03-09 file, typing:
>make

Now, the program can recognize SAINTexpress. You just need to give the path of the compiled files in bin folder: either SAINTexpress-spc or SAINT-express-int.

2. Downloading Omics Integrator

OmicsIntegrator predicts the interactome map of given proteins with their prize values. Here, prize values are taken as log2 fold change of the proteins. Omics Integrator generates the map by using some parameters, which you can change if you want more steiner proteins or not, or if you want to have more classes in your tree or not. 

You can download Omics Integrator from here: https://github.com/fraenkel-lab/OmicsIntegrator

After downloading it, you first need to download the required packages by typing:

> pip install -r requirements.txt

OmicsIntegrator-master file has GBM_case.py file, which  contains parameters, output paths and input paths. You need to change parameters and paths accordingly. Please, make the output path of knock-out as ...../result/KO. You can see an example of the changes in the attachment. 

3. Downloading Cytoscape, 

Cytoscape helps you to form the interaction map. It gets data from String, Omics Integrator, etc. 

You can download cytoscape from here: https://cytoscape.org/

After downloading it, you can upload the files in OmicsIntegrator/output folder. 



RUNNING THE PROGRAM 

As you can see, there are 2 files uploaded in the attachments. After downloading the programs, you need to download gprofiler package in python. The package is  found in https://pypi.org/project/gprofiler-official/

Then, you need to download upload 2 files. Before running the code, you need to change the parameters accordingly. 

You can give 2 different conditions at max. for GO annotation. More than that cannot be accepted, but if you have 3 or more conditions, please send us notifications. 

Then, you can use either a development environment like Pycharm or Jupyter to run interactome_map_parameters.py code, or you can directly run it by using command line. Please make sure all the required packages with their suitable versions are installed or imported into your environment. 

Then, you need to make sure that interactome_map_parameters.py and parameters.txt files are in the same directory. If they are not, the code is not going to work.  


You can reach input files in our example folder. 







