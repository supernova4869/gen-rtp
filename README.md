# GEN-RTP

## Author

Jiaxing Zhang at Tianjin University (zhangjiaxing7137@tju.edu.cn)

## Discriptions

An `rtp` and `hdb` file generator to be used together with the [Sobtop](http://sobereva.com/soft/Sobtop/) program. It could handle the following problems:

- Fix the wrong H names in the `mol2` file, change all heavy atom names to "element+id" and overwrite the origin file
- Generate the `rtp` file based on the `itp` file created by Sobtop program (follow the sobtop instructions)
- Remove the atoms and corresponding bonds, angles, and dihedrals in the rtp file by the rules in the next section
- Generate the hdb file based on the hydrogen rules
- TODO: Leave the atomtype and ffnonbonded instructions to separate files for user to add the items manually

### Bond, Angle, Dihedral, Improper Item Retain Rules

- Bond: for AMBER force field, retain "-C"; for GROMOS force field, retain "+N"
- Angle: Remove all exclude items
- Dihedral: Remove all exclude items
- Improper: Retain those only contain the connect atom
