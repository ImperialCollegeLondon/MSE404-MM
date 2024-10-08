import glob

def process_file(filename, output_file):
    nkpt = None
    
    with open(filename, 'r') as file:
        for line in file:
            if 'number of k points' in line:
                nkpt = line.split()[4]
            if line.startswith('!') and 'total' in line:
                if nkpt is not None:
                    total = line.split()[4]
                    output_file.write(f"{nkpt} {total}\n")

def main():
    output_filename = 'etot_v_nkpt.dat'
    files = glob.glob('*out_conv')
    
    with open(output_filename, 'w') as output_file:
        for filename in files:
            process_file(filename, output_file)

if __name__ == '__main__':
    main()
