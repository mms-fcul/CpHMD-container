# CpHMD-container
Git repository for divulging the st-CpHMD code in container fashion. 




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
