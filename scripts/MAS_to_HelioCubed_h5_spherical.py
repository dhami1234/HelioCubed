import h5py
import numpy as np
from datetime import datetime as dt
import time
from glob import glob
import os
import matplotlib.pyplot as plt

def toYearFraction(date):
    def sinceEpoch(date): return time.mktime(date.timetuple())
    s = sinceEpoch
    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)
    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    return date.year + yearElapsed/yearDuration

def read_maps(path):
    """
    Reads vr, rho, p, and optional br HDF5 file and returns the data and coordinates.
    """
    with h5py.File(path, "r") as f:
        vr = f["vr"][:]
        rho = f["rho"][:]
        p = f["p"][:]
        br = f["br"][:] if "br" in f else None
        phi = f["phi"][:]
        theta = f["theta"][:]
        r = f["r"][:]

    return vr, rho, p, br, phi, theta, r

def code_units_to_cgs(vr, rho, p, br, r):
    physical_vr = vr * 481.3711 * 1e5    # cm / s
    physical_rho = rho * 1.6726 * 1e-16  # g / cm^3
    physical_rho = physical_rho/(1.676 * 1e-24) # number density in cm^-3
    physical_p = p * 0.3875717           # dyn / cm^2
    physical_p = physical_p * 1e12       # convert to picodyn / cm^2
    physical_br = br * 1e6 if br is not None else None  # G -> microG
    physical_r = r * 6.96 * 1e10         # cm
    return physical_vr, physical_rho, physical_p, physical_br, physical_r

def R_omega_cms(r0_au, rotation_period_days=25.38):
    """
    Returns R*Omega (cm/s) for a given inner-boundary radius in AU.
    Uses Carrington sidereal rotation period by default.
    """
    AU_cm = 1.495978707e13
    R_cm = r0_au * AU_cm
    omega = 2.0 * np.pi / (rotation_period_days * 86400.0)
    return R_cm * omega

def resize_2d_interp(data, target_rows, target_cols):
    """
    Resize a 2D array using separable linear interpolation.
    """
    src_rows, src_cols = data.shape
    row_old = np.linspace(0.0, 1.0, src_rows)
    col_old = np.linspace(0.0, 1.0, src_cols)
    row_new = np.linspace(0.0, 1.0, target_rows)
    col_new = np.linspace(0.0, 1.0, target_cols)

    temp = np.empty((target_rows, src_cols), dtype=np.float64)
    for j in range(src_cols):
        temp[:, j] = np.interp(row_new, row_old, data[:, j])

    out = np.empty((target_rows, target_cols), dtype=np.float64)
    for i in range(target_rows):
        out[i, :] = np.interp(col_new, col_old, temp[i, :])

    return out

def resize_1d_interp(data, target_size):
    """
    Resize a 1D coordinate array using linear interpolation.
    """
    if data.size == target_size:
        return data

    old = np.linspace(0.0, 1.0, data.size)
    new = np.linspace(0.0, 1.0, target_size)
    return np.interp(new, old, data)

def make_output_name(input_file):
    bc_dir = os.path.dirname(input_file)
    model = os.path.basename(bc_dir)
    cr = os.path.basename(os.path.dirname(bc_dir))
    return os.path.join(bc_dir, f"{cr}_{model}_b.h5")

def convert_file(file, out_file=None, save_png=True, png_output_dir=None):
    r0 = 0.1
    num_components = 8
    time = toYearFraction(dt(2024, 3, 16, 12, 0, 0))
    phys_domain = [0.1, 0, 0, 0.1, 6.28319, 3.14159]
    step_const = [0, 1, 0]
    Br_scale_factor = 1.0
    target_rows, target_cols = 128, 128

    if out_file is None:
        out_file = make_output_name(file)
    if png_output_dir is None:
        png_output_dir = os.path.join(os.path.dirname(file), "frame_images_b")
    if save_png:
        os.makedirs(png_output_dir, exist_ok=True)

    vr, rho, p, br, phi, theta, r = read_maps(file)
    if br is None:
        raise ValueError(f"Input file does not contain br: {file}")

    data_Vr, data_rho, data_P, data_Br, r_cgs = code_units_to_cgs(vr, rho, p, br, r)
    data_Vr = resize_2d_interp(data_Vr, target_rows, target_cols)
    data_rho = resize_2d_interp(data_rho, target_rows, target_cols)
    data_P = resize_2d_interp(data_P, target_rows, target_cols)
    data_Br = resize_2d_interp(data_Br, target_rows, target_cols) * Br_scale_factor
    phi = resize_1d_interp(phi, target_rows)
    theta = resize_1d_interp(theta, target_cols)

    dim_siz1, dim_siz2 = data_rho.shape
    rho_flat = data_rho.flatten('F')

    theta = theta.flatten('F')
    phi = phi.flatten('F')

    Vr_flat = data_Vr.flatten('F')
    Vp_flat = Vr_flat * 0.0
    Vt_flat = Vr_flat * 0.0

    P_flat = data_P.flatten('F')

    Br_flat = data_Br.flatten('F')
    R_OMEGA_cms = R_omega_cms(r0_au=r0)
    theta_2d = np.tile(theta, (dim_siz1, 1))
    sin_theta_flat = np.sin(theta_2d).flatten('F')

    Bp_flat = -(R_OMEGA_cms * sin_theta_flat * Br_flat) / Vr_flat
    data_Bp = np.reshape(Bp_flat, (dim_siz1, dim_siz2), order='F')

    Br_sign = np.sign(Br_flat)
    Br_adj_mag = np.sqrt(np.maximum(Br_flat**2 - Bp_flat**2, 0.0))
    Br_flat = Br_sign * Br_adj_mag
    data_Br = np.reshape(Br_flat, (dim_siz1, dim_siz2), order='F')
    Bt_flat = Vr_flat * 0.0

    data_final = np.concatenate((rho_flat, Vr_flat, Vp_flat, Vt_flat, P_flat, Br_flat, Bp_flat, Bt_flat))
    dtheta = np.full(dim_siz2, np.pi / dim_siz2)

    with h5py.File(out_file, "w") as data_file:
        data_file.create_dataset("data0", data=data_final)
        data_file["data0"].attrs["time"] = 0
        data_file.attrs['domain'] = [dim_siz1, dim_siz2]
        data_file.attrs['num_components'] = num_components
        data_file.attrs['num_datasets'] = 1
        data_file.attrs['time'] = time
        data_file.attrs['r0'] = r0
        times = 0
        data_file.create_dataset("datasets_time", data=times)
        grp = data_file.create_group("geometry")
        grp.attrs['phys_domain'] = phys_domain
        grp.attrs['step_const'] = step_const
        grp.create_dataset("dtheta", data=dtheta)
        grp.create_dataset("theta", data=theta)
        grp.create_dataset("phi", data=phi)

    if save_png:
        save_frame_png(png_output_dir, data_rho, data_Vr, data_P, data_Br, data_Bp)

    return {
        "out_file": out_file,
        "png_output_dir": png_output_dir,
        "data_shape": data_final.shape,
        "domain": [dim_siz1, dim_siz2],
        "dtype": data_final.dtype,
        "num_components": num_components,
        "phys_domain": phys_domain,
        "time": time,
    }

def save_frame_png(png_output_dir, data_rho, data_Vr, data_P, data_Br, data_Bp):
    if data_Br is not None:
        fig, axs = plt.subplots(2, 3, figsize=(28, 10))
    else:
        fig, axs = plt.subplots(2, 2, figsize=(22, 10))

    pcm0 = axs[0, 0].pcolormesh(data_rho.T, shading='auto')
    axs[0, 0].set_title(f"Density Frame")
    axs[0, 0].invert_yaxis()
    fig.colorbar(pcm0, ax=axs[0, 0])

    pcm1 = axs[0, 1].pcolormesh(data_Vr.T, shading='auto')
    axs[0, 1].set_title(f"Radial Velocity Frame")
    axs[0, 1].invert_yaxis()
    fig.colorbar(pcm1, ax=axs[0, 1])

    if data_Br is not None:
        pcm2 = axs[0, 2].pcolormesh(data_P.T, shading='auto')
        axs[0, 2].set_title(f"Pressure Frame")
        axs[0, 2].invert_yaxis()
        fig.colorbar(pcm2, ax=axs[0, 2])

        pcm3 = axs[1, 0].pcolormesh(data_Br.T, shading='auto')
        axs[1, 0].set_title(f"Radial Magnetic Field Frame")
        axs[1, 0].invert_yaxis()
        fig.colorbar(pcm3, ax=axs[1, 0])

        pcm4 = axs[1, 1].pcolormesh(data_Bp.T, shading='auto')
        axs[1, 1].set_title(f"Azimuthal Magnetic Field Frame")
        axs[1, 1].invert_yaxis()
        fig.colorbar(pcm4, ax=axs[1, 1])

        axs[1, 2].axis('off')
    else:
        pcm2 = axs[1, 1].pcolormesh(data_P.T, shading='auto')
        axs[1, 1].set_title(f"Pressure Frame")
        axs[1, 1].invert_yaxis()
        fig.colorbar(pcm2, ax=axs[1, 1])

    plt.tight_layout()
    plt.savefig(os.path.join(png_output_dir, f"frame_{0:04d}.png"))
    plt.close()

def convert_tree(root_directory, save_png=True):
    files = sorted(glob(os.path.join(root_directory, "**", "vr_rho_p_br_r0.h5"), recursive=True))
    print("Number of br boundary files:", len(files))
    results = []

    for i, file in enumerate(files, start=1):
        print(f"[{i}/{len(files)}] Converting: {file}")
        result = convert_file(file, save_png=save_png)
        results.append(result)
        print("  Output file:", result["out_file"])

    return results

if __name__ == "__main__":
    save_png = True
    batch_root = '/Users/talwindersingh/Library/CloudStorage/Dropbox-GSUDropbox/Talwinder Singh/Reza_project'
    results = convert_tree(batch_root, save_png=save_png)
    if results:
        last = results[-1]
        print("Time in fraction of year:", last["time"])
        print("Conversion completed.")
        print("Last output file:", os.path.abspath(last["out_file"]))
        print("Data shape:", last["data_shape"])
        print("Data domain:", last["domain"])
        print("Data type:", last["dtype"])
        print("Number of components:", last["num_components"])
        print("Physical domain:", last["phys_domain"])
