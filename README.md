# The NEF (NMR Exchange Format)

    1.Gutmanas, A. et al. NMR Exchange Format: a unified and open standard for representation of NMR restraint data. Nat Struct Mol Biol 22, 433–434 (2015).
  
The purpose of NEF is to provide a unique implementation of NMR data (e.g. restraints). In order for this format to be used, computational software that uses NMR data should be able to read, process and write this type of format files.

# MELD 
Our purpose as developers of the MELD package has been to combine simulations and experimental data to solve problems that could not be tackled by either approach alone or by accelerate the time it gets to go from data to structural knowledge. We have participated in several CASP-NMR events with success and our purpose here is to make the intragration of standard NMR information into MELD easier.

    1.MacCallum, J. L., Perez, A. & Dill, K. Determining protein structures by combining semireliable data with atomistic physical models by Bayesian inference. Proc National Acad Sci 112, 6985–6990 (2015).
    2.Robertson, J. C. et al. NMR‐assisted protein structure prediction with MELDxMD. Proteins Struct Funct Bioinform 36, D402 (2019).
  
# MELDNMRParser
This is a set of scripts for reading, processing and writing NEF files. The scripts start by reading an initial NEF file provided by an NMR collaborator, processing it and producing the required MELD files to run structure determination simulations. The scripts also write out a new NEF file that contains the original information as well as how the restraints are being enforced inside the MELD simulation in a way that is compatible with the NEF format.

