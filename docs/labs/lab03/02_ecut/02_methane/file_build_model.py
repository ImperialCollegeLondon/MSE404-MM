##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#                           This is a python file that will generate multiple input files for a convergence test.                            #
#                                                                                                                                            #
#                   How to use: Copy and paste the text from your input file to common_content_template as shown below.                      #
#                               Make sure you are putting a space both sides of = sign when pasting.                                         #
#                                                                                                                                            #
#                                                   How to run: python3 ecut_build.py                                                        #
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################


# Directory where you want to create the files
output_directory = "./"

# Number of files to create
num_files = 10

# Common content template with a placeholder for the number
common_content_template = """
 &CONTROL
    pseudo_dir = '.'
    disk_io = 'none'
 /

 &SYSTEM
    ibrav = 1
    A = 15.0
    nat = 5
    ntyp = 2
    ecutwfc = {}
 /

 &ELECTRONS
    ! This is the default
    conv_thr = 1.0E-6
 /

ATOMIC_SPECIES
 C  12.011  C.pz-vbc.UPF
 H   1.008  H.pz-vbc.UPF

ATOMIC_POSITIONS angstrom
 C  0.0       0.0       0.0
 H  0.0       0.0       1.089
 H  1.026719  0.000000 -0.363000
 H -0.513360 -0.889165 -0.363000
 H -0.513360  0.889165 -0.363000

K_POINTS gamma
"""

# Loop to create the files
for i in range(1, num_files + 1):
    # Define the content for each file with the number replaced
    calc = round(5+5*i)
    content = common_content_template.format(calc)

    # Generate the file name
    file_name = f"{output_directory}/scf.mol.{str(i).zfill(3)}.in"

    # Open and write to the file directly
    with open(file_name, 'w') as file:
        file.write(content)
