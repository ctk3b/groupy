from groupy.general import *
from groupy.monolayers import *


# --- user specifications ---
file_name = ('data_files/mpc_0.2_35_shortened.lammpstrj')
n_planes = 100  # number of flux planes

# --- system info ---
# atom indices of regions
botsub_ids = range(3200)
botfilm_ids = range(3200, 4180)
water_ids = range(4180, 6880) + range(11060, 13760)
topsub_ids = range(6880, 10080)
topfilm_ids = range(10080, 11060)

system_info = {'botsub': botsub_ids,
               'botfilm': botfilm_ids,
               'water': water_ids,
               'topsub': topsub_ids,
               'topfilm': topfilm_ids}

# --- main ---
# read first frame to get system dimensions
with open(file_name, 'r') as trj:
    xyz, _, _, _ = read_frame_lammpstrj(trj)

# calculate surface area
x_max = xyz[system_info['botsub']][:, 0].max()
x_min = xyz[system_info['botsub']][:, 0].min()
y_max = xyz[system_info['botsub']][:, 1].max()
y_min = xyz[system_info['botsub']][:, 1].min()
area = (x_max - x_min) * (y_max - y_min)

# determine bounds on planes based on substrate surface coords
z_min = xyz[system_info['botsub']][:, 2].max()
z_max = xyz[system_info['topsub']][:, 2].min()
planes = np.linspace(z_min, z_max, n_planes)

# where the magic happens
fluxes, steps = calc_flux(file_name,
        system_info,
        planes,
        area,
        max_time=1e6)
# convert atoms / (angstrom^2 * fs) to h2o_molecules / (nm^2 * ps)
fluxes = fluxes * (100 * 1000 / 3)

# calculate time averaged flux across each plane
avg_fluxes = np.empty(shape=(n_planes))
for i in range(fluxes.shape[1]):
    p = np.polyfit(steps, fluxes[:, i], 1)
    avg_fluxes[i] = p[1]

# show fluxes at respective z-coords
fig = plt.figure()
plt.plot(avg_fluxes, planes)
plt.xlabel(ur'Flux (molecules/nm$\mathregular{^2\cdot}$ ps)')
plt.ylabel(ur'z (\u00c5)')
fig.savefig('flux_vs_z.pdf', bbox_inches='tight')
