# open file as a file object
with open("C_diamond_base.in") as f:
    # read file to content [list]
    content = f.readlines()
    # the last line (-1) is where we specify k-point density
    kpt = content[-1].strip().split()

# we want to test the k-point density from 2 to 10:
# note that range function only outputs [2,11).
for num_kpt in range(2,11):
    # lets set the k-point density along ka; kb; kc.
    kpt[0] = str(num_kpt); kpt[1] = str(num_kpt); kpt[2] = str(num_kpt); 
    # construct the k-point line
    kpt_final = ' '.join(kpt)
    # set the k-point line to teh content list
    content[-1] = kpt_final

    # write the files 
    with open(f"C_diamond_{num_kpt}.in",'w') as f:
        f.writelines(content)

    print(f"{num_kpt}...done.")

