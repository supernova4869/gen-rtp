# GEN-RTP

## Author

Jiaxing Zhang at Tianjin University (zhangjiaxing7137@tju.edu.cn)

## Discriptions

A `.rtp` and `.hdb` file generator to be used together with the [Sobtop](http://sobereva.com/soft/Sobtop/) program. It could handle the following problems:

- Fix the wrong H names in the `.mol2` file and generate `newXX.mol2`
- Generate the rtp file based on the itp file created by Sobtop program (following the sobtop instructions)
- Remove the atoms and corresponding bonds, angles, and dihedrals in the rtp file (TODO: by rules)
- Generate the hdb file based on the hydrogen rules
- TODO: Leave the atomtype and ffnonbonded instructions to separate files for user to add the items manually
