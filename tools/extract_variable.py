import argparse
import json
import numpy as np

variable_table = {
    "temperature" : 0,
    "vorticity" : 1,
    "streamfunction" : 2,
    "xi" : 7
}

def main():
    parser = argparse.ArgumentParser(description='Extract chosen variable')
    parser.add_argument('filename', help='file to open')
    parser.add_argument('--output',
                        help='Output file',
                        required=True)
    parser.add_argument('--constants', help='constants file',
                        required=True)
    parser.add_argument('--variable', help='variable to extract',
                        required=True)
    args = parser.parse_args()

    constants_file = open(args.constants, "r")
    constants = json.load(constants_file)

    n_modes = constants["nN"]
    n_z_gridpoints = constants["nZ"]

    data = np.fromfile(args.filename, dtype=np.dtype(np.cdouble))

    if args.variable in variable_table:
        varidx = variable_table[args.variable]
    else:
        print("Error: ", args.variable, "not one of", variable_table.keys())

    variable = data[
        (varidx)*n_modes*n_z_gridpoints:
        (varidx+1)*n_modes*n_z_gridpoints
    ]

    with open(args.output, "w") as fp:
        variable.tofile(fp)

main()

