##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#                           This is a python file that will generate multiple input files for a convergence test.                            #
#                                                                                                                                            #
#                   How to use: Copy and paste the text from your input file to common_content_template.                                     #
#                                       Replace ecutwfc = (some number) with ecutwfc = {}                                                    #
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
