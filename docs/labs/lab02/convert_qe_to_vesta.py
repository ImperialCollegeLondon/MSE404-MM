##############################################################################################################################################
##############################################################################################################################################
#                           This is a python file that will convert a QE input file to a VESTA readable POSCAR file.                         #
#                                                                                                                                            #
#                   How to use: You must set ibrav=0 and have the CELL_PARAMETERS defined in your input file. Additionally                   #
#                               your ATOMIC_POSITIONS card must be in cartesian coordinates in units of angstrom (Å).                        #
#                                                                                                                                            #
#                                                   How to run: python3 convert_qe_to_vesta filename.in                                      #
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
import argparse
import re

def read_qe_input(file_path):
    """Reads Quantum ESPRESSO input file and extracts lattice vectors and atomic positions."""
    with open(file_path, 'r') as f:
        lines = f.readlines()

    lattice_vectors = []
    atomic_positions = []
    atomic_species = []

    # Flags for reading
    reading_lattice = False
    reading_positions = False

    for line in lines:
        # Extract lattice vectors
        if 'CELL_PARAMETERS' in line:
            reading_lattice = True
            continue
        if reading_lattice:
            if len(line.strip()) == 0:  # Stop when empty line encountered
                reading_lattice = False
            else:
                lattice_vectors.append(line.strip())

        # Extract atomic species and positions
        if 'ATOMIC_POSITIONS' in line:
            reading_positions = True
            continue
        if reading_positions:
            if len(line.strip()) == 0 or 'End of' in line:  # Stop when empty line or section end encountered
                reading_positions = False
            else:
                tokens = line.strip().split()
                atomic_species.append(tokens[0])
                atomic_positions.append(tokens)

    return lattice_vectors, atomic_species, atomic_positions

def generate_poscar(lattice_vectors, atomic_species, atomic_positions, output_file='POSCAR'):
    """Generates POSCAR file from the extracted lattice vectors and atomic positions."""

    # Get unique atomic species and count how many atoms for each species
    species_count = {}
    for species in atomic_species:
        species_count[species] = species_count.get(species, 0) + 1

    # Write POSCAR file
    with open(output_file, 'w') as f:
        f.write("Generated from QE output\n")
        f.write("1.0\n")  # Scaling factor

        # Write lattice vectors
        for vec in lattice_vectors:
            f.write(f"{vec}\n")

        # Write atomic species and their counts
        f.write(" ".join(species_count.keys()) + "\n")
        f.write(" ".join(map(str, species_count.values())) + "\n")

        # Cartesian coordinates
        f.write("Cartesian\n")

        # Write atomic positions in Cartesian coordinates
        for position in atomic_positions:
            f.write(f"{position[1]} {position[2]} {position[3]}\n")

    print(f"POSCAR file generated: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='A script to calculate the charge center by integrating along lattice directions.')
    parser.add_argument('INFILE', nargs='?', default='', help="name of the quantum espresso input file.")

    args = parser.parse_args()

    # Path to the Quantum ESPRESSO output file
    qe_input_file = args.INFILE

    # Extract lattice vectors, atomic species, and atomic positions
    lattice_vectors, atomic_species, atomic_positions = read_qe_input(qe_input_file)

    # Generate POSCAR file
    generate_poscar(lattice_vectors, atomic_species, atomic_positions)

if __name__ == "__main__":
    main()
