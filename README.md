# CF_CLUSTERSpy

 (c) Author  :  Woongkyu Jee, Scott M. Woodley;
 Contact     :  woong.jee.16@ucl.ac.uk / wldndrb1@gmail.com;
 Affiliation :  University College London, Department of Chemistry;
 Address     :  University College London, Kathleen Lonsdale Materials Chemsitry, 20 Gordon Street, London WC1H 0AJ.
 

#### Description
  
CF_CLUSTSERSpy is a code, designed for comparing structures of moleclues or clusters.

There are two different types of Python written sources, "CF_CLUSTERSpy_MAIN.py" and "CF_CLUSTERSpy_MOD.py".

"CF_CLUSTERSpy_MOD.py" is the source, which contains the methods used in "CF_CLUSTERSpy_MAIN.py".

"CF_CLUSTERSpy_MAIN.py" is the main executable that writes the output of the user requested work.

For instance, to execute the code user can type the following command line,

(***Before trying the following command lines, make sure that you have Python version above 3 and Numpy package is installed***)

$ python CF_CLUSTERSpy_MAIN.py   /path/to/1st/input/file   /path/to/2nd/input/file

In the above command line, you can find that the code requires two arguments, which are files (or their paths) including structural data of molecules or clusters that user wants to compare.

***The format of input files must be "~.xyz", which are generally used in many different structural visualisation packages, e.g., VESTA,
and the atoms/ions in those two different xyz files must be in the same order***

Once you successfully finshed cloning the repository, you can find two different sample input files, '1_AIMS_N6.xyz' and '1_SLAM_N6.xyz',
which are structural data of (SnO)N clusters with size of N=6 calculated by using DFT (density functional theory) and SLAM force field, respectively.

To compare those two different structures, one can try the following command line,

$ python CF_CLUSTERSpy_MAIN.py 1_AIMS_N6.xyz 1_SLAM_N6.xyz 

Afterwards there will be four new files generated in the same directory, which are 

1. CF.out      : main output file including general overview of structural anaylsis.

2. CF_1.xyz    : rotated / recentred structural data (in xyz format) of the first input file (for this case '1_AIMS_N6.xyz').

3. CF_2.xyz    : rotated / recentred structural data (in xyz format) of the second input file (for this case '1_SLAM_N6.xyz').

4. CF_3.xyz    : plot of both 'CF_1.xyz' and 'CF_2.xyz' in the same xyz file for visual comparison.

Note that the 'rotated / recentred' implies the original refrence frames (i.e., coordinate systems) of two input molecules / clusters are rotated and translated
according to the information of their inertia tensors, and those results are given in the files, CF_1.xyz and CF_2.xyz. 

There are three major figures in the general output file, 'CF.out', which users might want to check.

1. RMS CONFIG     : refers to the 'root-mean-square deviation' of two input atom/ion configurations, which are rotated / recentred.
		    If this figure is close to zero, then it implies two structures are in matching, i.e., if the value is < 0.0025.

2. RMS PMOI       : refers to the 'root-mean-square deviation' of principle moment of inertia of two input atom/ion configurations, which are rotated / recentred.
		    This figure can be used to find out how the atoms or ions are scatterd in the space or used as supportive data correcting 'RMS CONFIG'. 

3. Scaling Factor : refers to how the first structure is expanded with respect to the second structure.


######Thank you very much !
