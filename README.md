# CpHMD-container
Git repository for divulging the st-CpHMD code in container format. Having such tool will allow for CpHMD to be ran in any type of machine as long as <ins>**singularity**<\ins> is installed. 

The container has been built in Ubuntu 24.04 and containes the following:
- Ubuntu 24.04 - base of the container
- GROMACS 2024.3 (CPU only)
- CpHMD running scripts
- CpHMD Force Fields - GROMOS54a7, CHARMM36, AMBER14SB
- Delphi v
- Petit 1.6.1

Environments currently tested for CpHMD container calculations: 
- [X] Linux environment
- [X] Windows environment (using WSL)
- [ ] Mac OS

# How to run CpHMD 

Running CpHMD requires several steps to ensure your system is ready for production. 
Three steps are usually required for a safe CpHMD run: 
1. Preparing your starting pdb file [Wiki:System Preparation](https://github.com/mms-fcul/CpHMD-container/wiki/System-Preparation)
2. Running Minimization and Initialization
3. Preparing the inputs for the CpHMD run [Wiki:Editing-the-CpHMD.settings](https://github.com/mms-fcul/CpHMD-container/wiki/Editing-the-CpHMD.settings-to-your-taste)

With all files prepared, starting the simulation is as simple as running your production folder.
```
singularity exec --bind <your home directory> <CpHMD container> /CpHMD/scripts/CpHMD.sh  ./CpHMD-(basic|advanced).settings
```
How to manage the simulation submission, paralelization and background tasking is left on the hands of the user. This allows for a better adaptation on however each computing infrastructure is set, requiring only integration on the prefered workflow.

--------------------------------------------------------------------------------------------------------------------------------------
## Key available apps
### System preparation 
- **pdb2cphmd** [Wiki_page](https://github.com/mms-fcul/CpHMD-container/wiki/System-Preparation)

    - This app takes a pdb file and changes the residue names to be compliant with the CpHMD titratable residues. This change is needed since CpHMD uses tautomers. EX.: Asp has 4 tautomers on the carboxylic acid (1 front facing and 1 back facing hydrogen on each oxygen atom) hence ASP ranges from AS0 to AS4. Usually Asp is exchanged with AS4 (which is the charged tautomer).
  
    _How to use the app:_
    ```
    singularity run --app pdb2cphmd CpHMD.sif [-h for help] -p <pdb file> -f <force-field (G54a7pH;CHARMM36pH;Amber14SBpH)> -r  <residue names to change for CpHMD - ASP GLU HIS, etc>"
    ```
      
### File extraction 
- **extract-force-fields**

    - This app copies all the force fields available within the container to your current directory. This gives user access to the force fields used in the CpHMD container and allows for customization.
  
    _How to use the app:_
    ```
    singularity run --app extract-force-fields CpHMD.sif 
    ```
- **extract-tautomers**

    - This app copies the tautomer definition files (Sts) used in the CpHMD. These files are crucial to define the different charges and atoms of each residue tautomer.  
  
    _How to use the app:_
    ```
    singularity run --app extract-tautomers CpHMD.sif 
    ```
