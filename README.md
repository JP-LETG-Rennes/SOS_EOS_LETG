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
      4. Launching via terminal SOS_EOS_Traitement_V2.py

     
-------------

Citation
---------




