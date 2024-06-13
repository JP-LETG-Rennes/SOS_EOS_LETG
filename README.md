Classification Start of Season, End of Season
======

Overview
-----


About
-------------------------


Installation
-------------

For ease of use, we recommend that you first install Git and Anaconda.

Pre-requis:

      Git
      Anaconda or miniconda

Git Windows:

      #For Windows distribution, you can install with this link :
      https://git-scm.com/download/win
We recommmend to choose the Standalone Installer

Git Linux (Debian/Ubuntu):

      # You can install git with this command:
      sudo apt install git-all
      
Git Linux (Fedora):
      
      # You can install git with this command:
      sudo dnf install git-all

Anaconda: 

      # You can download and install anaconda with this link :
      https://www.anaconda.com/download/success

Library Dependency
---------------
Once the prerequisites have been installed, you can launch the next section via Anaconda Prompt

```
# Clone the repo
git clone JP-LETG-Rennes/SOS_EOS_LETG.git
cd SOS_EOS_LETG

# Create a new python environnement with conda  
conda create --name SOS_EOS python=3.9
conda activate SOS_EOS

# Prepare pip
conda install pip
pip install --upgrade pip

# Install requirements
pip install -r requirements_envs.txt

```

Getting Started
---------------
If Using Terminal:

      1. Download build from source. 
      2. Activate Conda Environnement Hiclass
      3. Enter variables in configuration file : configuration_SOS_EOS_LETG.json
      4. Launching via terminal SOS_EOS.py

     
Simplified usage with notebook : 
---------------
For those unfamiliar with the use of git, it is possible to avoid it by typing the following in the conda prompt :  

```
conda create --name Hiclass python=3.8 --yes 

activate Hiclass

pip install hiclass numpy pandas geopandas matplotlib rasterio scipy scikit-learn pyproj scipy notebook seaborn xarray rioxarray shap 

jupyter notebook

```

Typing this in conda terminal will create a new environment named "Hiclass", activate it, install all necessary dependencies and open jupyter notebook. All that is left to do is download the .ipynb file, navigate to it through the jupyter interface and open it. 

Example dataset
---------------

The 5 X 5 km site study site is located along the Gironde estuary (France) in the Calupeyre catchment area (45.47°N, 1.08°W). The site comprises coastal dunes in the western part, a marsh in the central part and limestone hillsides in the eastern part.

The dataset consists of:
- a raster file containing the predictive variables; 
- a table giving the values of each predictive variable and the vegetation class for 1,629 field samples collected over the whole watershed.

The raster Calupeyre_pred_variables.tif was projected in the French projection system (Lambert 93, code EPSG 2154) at 10  10 m spatial resolution. It includes 6 bands characterizing the micro-topography (Panhelleux et al., 2023) - that were derived from an airborne digital terrain model (RGE ALTI ® IGN) -  as well as 27 bands characterizing the first three components of the functional principal component analysis for each of the 9 spectral bands (band 2, 3, 4, 5, 6, 7, 8, 10 & 11) derived from a pluriannual (2015-2021) multispectral satellite Sentinel-2 time series (® ESA). The method used to generate the functional principal component analysis components was based on Pesaresi et al (2022).

The table "Calupeyre_extracted_raster_values.csv" is in csv format (comma separator). It contains for each sample (row) the value of each predictor variable (column) as well as the corresponding natural habitat class in the European EUNIS typology (Davies et al., 2004) for each of the three hierarchical levels (level1, level2, level3). The vegetation samples were collected in the field and archived in various open access databases:
- The French national inventory of natural heritage (Poncet, 2013) accessed on the website of the national natural history museum (https://inpn.mnhn.fr/accueil/index?lg=en), 
- the French national forest inventory (Hervé, 2016), accessed from the IGN website (https://inventaire-forestier.ign.fr/dataifn/), 
- samples collected as part of wetland inventories (Gayet et al., 2022).

Each sample was automatically assigned to level 3 of the EUNIS typology based on their floristic composition and environmental characteristics according to the approach detailed in Chytrý et al. (2021). In addition, crop samples taken from the European land parcel information system (LPIS) - available on the IGN ® website (https://geoservices.ign.fr/rpg) - were assigned to the Intensive unmixed crops habitat (EUNIS code I1.1). It should be noted that, for the purposes of this tutorial and for genericity, the characterization of some samples was voluntary degraded to level 1 or 2 of the EUNIS typology.


References 
-------------

Citation
---------




