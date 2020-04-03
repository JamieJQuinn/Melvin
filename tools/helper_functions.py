import numpy as np

def extract_variable(data, varidx, n_modes, n_z_gridpoints):
    return np.transpose(data[(varidx)*n_modes*n_z_gridpoints:(varidx+1)*n_modes*n_z_gridpoints]\
                           .reshape(n_modes, n_z_gridpoints))

def get_spatial_data(data, varidx, constants, is_sine):
    n_modes = constants["nN"]
    n_z_gridpoints = constants["nZ"]
    aspect_ratio = constants["aspectRatio"]
    var = extract_variable(data, varidx, n_modes, n_z_gridpoints)
    if constants["horizontalBoundaryConditions"] == "periodic":
        var_spatial = np.fft.irfft(var)*(2.0*n_modes-1)
    else:
        n_x_gridpoints = int(aspect_ratio*n_modes)
        x_axis = np.linspace(0, aspect_ratio, num=n_x_gridpoints)
        if is_sine:
            sine = np.sin(np.pi/aspect_ratio*np.outer(np.arange(n_modes), x_axis))
            var_spatial = np.real(np.dot(var, sine))
        else:
            cosine = np.cos(np.pi/aspect_ratio*np.outer(np.arange(n_modes), x_axis))
            var_spatial = np.real(np.dot(var, cosine))

    return var_spatial

