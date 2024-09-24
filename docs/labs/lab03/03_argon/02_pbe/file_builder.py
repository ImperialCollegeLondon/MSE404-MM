##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#                           This is a python file that will generate multiple input files for a convergence test.                            #
#                                                                                                                                            #
#                   How to use: Copy and paste the text from your input file to common_content_template as shown below.                      #
#                               Make sure you are putting a space both sides of = sign when pasting.                                         #
#                                                                                                                                            #
#                                                   How to run: python3 analysis.py                                                          #
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################


# Directory where you want to create the files
output_directory = "./"

# Number of files to create
num_files = 20

# Common content template with a placeholder for the number
common_content_template = """ &CONTROL
    pseudo_dir = '.'
    disk_io = 'none'
 /

 &SYSTEM
    ibrav = 1
    A = 12.0
    nat = 2
    ntyp = 1
    ecutwfc = 30.0
 /

 &ELECTRONS
    ! This is the default
    conv_thr = 1.0E-6
 /

ATOMIC_SPECIES
 Ar  39.948  Ar.pbe-n-rrkjus_psl.1.0.0.UPF

ATOMIC_POSITIONS angstrom
 Ar  0.0       0.0       0.0
 Ar  0.0       0.0       {}

K_POINTS gamma
"""

# Loop to create the files
for i in range(1, num_files + 1):
    # Define the content for each file with the number replaced
    calc = round(3+0.1*i, 3)
    content = common_content_template.format(calc)

    # Generate the file name
    file_name = f"{output_directory}/scf.mol.{str(i).zfill(3)}.in"

    # Open and write to the file directly
    with open(file_name, 'w') as file:
        file.write(content)
