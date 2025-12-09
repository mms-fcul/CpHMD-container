In the wiki for this project:

- https://github.com/mms-fcul/CpHMD-container/wiki

There's a very nice and detailed tutorial that can be used to
  get started using the a Singularity container with the Sto-
  chastic titatrion CpHMD code on different system types.
So far, the tutorials overview:
- CpHMD: Protein in Solution
- CpHMD: Protein + Membrane system
- CpHMD: Protein + Solute as include
- CpHMD: Protein + Solute in the FF

Please refer to these tutorials if you need any guidance on the
  steps needed. To download them, ensure that you have the CpHMD.sif
  file in the current directory, and run:

singularity run --app get-tutorials CpHMD.sif
