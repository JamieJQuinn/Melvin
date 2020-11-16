import numpy as np
from PIL import Image
from matplotlib import cm
from scipy import interpolate
import argparse
import json
import colorsys

from helper_functions import extract_variable

# def color(x, base_color):
    

def main():
    """main function"""
    parser = argparse.ArgumentParser(description='Print variable contents')

    parser.add_argument('filename', help='file to open')
    parser.add_argument('--output',
                        help='Output file')
    parser.add_argument('--constants', help='constants file')
    parser.add_argument('--color', help='Base colour for shading')
    parser.add_argument('--x_resolution', help='Resolution in x', type=int)
    parser.add_argument('--z_resolution', help='Resolution in z', type=int)
    parser.add_argument('--levels', help='Number of levels to round data to (0 to disable)', type=int, default=0)
    parser.add_argument('--index', help='Index of variable to print', type=int, default=0)
    args = parser.parse_args()

    constants_file = open(args.constants, "r")
    constants = json.load(constants_file)

    is_ddc = constants["isDoubleDiffusion"]
    n_modes = constants["nN"]
    nZ_in = constants["nZ"]
    aspect_ratio = constants["aspectRatio"]

    if "horizontalBoundaryConditions" not in constants:
        constants["horizontalBoundaryConditions"] = "impermeable"

    data = np.fromfile(args.filename, dtype=np.dtype(np.cdouble))
    temp = extract_variable(data, args.index, n_modes, nZ_in)

    nX_out = args.x_resolution
    nZ_out = args.z_resolution

    if constants["horizontalBoundaryConditions"] == "periodic":
        temp_spatial = np.fft.irfft(temp, nX_out)*(2.0*nX_out-1)
    else:
        x_axis = np.linspace(0, aspect_ratio, num=nX_out)
        cosine = np.cos(np.pi/aspect_ratio*np.outer(np.arange(n_modes), x_axis))
        temp_spatial = np.real(np.dot(temp, cosine))

    z_in = np.linspace(0, 1, nZ_in)
    interpolation = interpolate.interp1d(z_in, temp_spatial, axis=0)
    z_out = np.linspace(0, 1, nZ_out)
    upscaled_temp = interpolation(z_out)

    levels = args.levels
    if args.levels:
        upscaled_temp = 1.0/levels*np.round(upscaled_temp*levels)

    cmap = cm.get_cmap('magma_r')

    im = Image.fromarray(cmap(upscaled_temp, bytes=True))
    im.save(args.output)

main()
