import numpy as np
import rasterio
import matplotlib.pyplot as plt
import pyvista as pv
import pyspaceaware as ps
import seaborn as sns
import json

station_lat_rad = 0.574213323906134
station_lon_rad = -1.84190413727077
station_alt_m = 2.22504
station_height_above_terrain_km = 3/1e3

station_ecef = ps.lla_to_itrf(station_lat_rad, station_lon_rad, station_alt_m)

s_corner_lat, w_corner_lon = 32, -106
n_corner_lat, e_corner_lon = 33, -105
with rasterio.open(f'data/{"n" if s_corner_lat > 0 else "s"}{abs(s_corner_lat)}_{"w" if w_corner_lon < 0 else "e"}{abs(w_corner_lon)}_1arc_v3.tif') as src:
    elev_grid = src.read().squeeze().T / 1e3
    lat = n_corner_lat - (n_corner_lat-s_corner_lat) * np.linspace(0,1,elev_grid.shape[0])
    lon = w_corner_lon - (w_corner_lon-e_corner_lon) * np.linspace(0,1,elev_grid.shape[1])
    lat_grid, lon_grid = np.meshgrid(lat, lon)


ecef = ps.lla_to_itrf(np.deg2rad(lat_grid.flatten()), np.deg2rad(lon_grid.flatten()), elev_grid.flatten())
station_terrain_alt_km = elev_grid.flatten()[np.argmin(ps.vecnorm(station_ecef - ecef))]
station_ecef = ps.lla_to_itrf(station_lat_rad, station_lon_rad, station_terrain_alt_km + station_height_above_terrain_km)

ecef_to_enu_rotm = ps.ecef_to_enu(station_ecef)
station_enu = (ecef_to_enu_rotm @ station_ecef.T).T
enu = (ecef_to_enu_rotm @ ecef.T).T

enu_los = enu - station_enu
x_grid = enu_los[:,0].reshape(lat_grid.shape, order="F")
y_grid = enu_los[:,1].reshape(lat_grid.shape, order="F")
z_grid = enu_los[:,2].reshape(lat_grid.shape, order="F")


ps.tic()
dem = pv.StructuredGrid(x_grid, y_grid, z_grid).extract_surface().triangulate()
dem["Elevation [km]"] = elev_grid.flatten(order="F")
is_local = ps.vecnorm(dem.points).flatten() < 1e0
is_above = dem["Elevation [km]"] > station_terrain_alt_km - station_height_above_terrain_km
dem = dem.extract_points(is_above | is_local).extract_surface().triangulate()
dem.compute_normals(cell_normals=False, point_normals=True, inplace=True)
normals = dem["Normals"]
is_facing_station = ps.dot(ps.hat(dem.points), normals).flatten() < 0
is_ridge = ps.dot(ps.hat(dem.points), normals).flatten() < 0.1
is_close = ps.vecnorm(dem.points).flatten() < 1e1
is_local = ps.vecnorm(dem.points).flatten() < 1e0
ids = np.argwhere((is_ridge & is_facing_station & is_close) | is_local)
dem = dem.extract_points(ids).extract_surface().triangulate()
ps.toc()


naz_mask = 1000
d_el = np.deg2rad(1)/10
az_ray = np.linspace(0, 2*np.pi, naz_mask)
r_ray = np.max(ps.vecnorm(dem.points)) * np.ones_like(az_ray)
origins = np.zeros((naz_mask, 3))
el_horizon = np.array([-1.0 for x in range(naz_mask)], dtype=np.float64)
_, el, _ = ps.cart_to_sph(dem.points[:,0], dem.points[:,1], dem.points[:,2])
el_ray = np.ones_like(az_ray) * np.max(el)
max_ring_enu = np.vstack(ps.sph_to_cart(az_ray, el_ray)).T


hat_los_enu = ps.hat(enu_los)
for i, mre in enumerate(max_ring_enu):
    # dp = hat_los_enu[:,0]*mre[0] + hat_los_enu[:,1]*mre[1] + hat_los_enu[:,2]*mre[2]
    dp = ps.dot(hat_los_enu, mre)
    min_ang_to_dem = np.min(np.arccos(dp))
    el_ray[i] = el_ray[i] - min_ang_to_dem + d_el/2

ray_prev_collided = np.array([False for _ in range(naz_mask)])
i = 0
while -1.0 in el_horizon:
    ps.tic()
    ray_searching = np.array([True if n == -1.0 else False for n in el_horizon])
    ray_inds_searching = np.argwhere(ray_searching)
    dests = np.vstack(ps.sph_to_cart(az_ray, el_ray, r_ray)).T
    (_, iray_inds, _) = dem.multi_ray_trace(origins[ray_searching,:], dests[ray_searching,:], first_point=True)
    irays = ray_inds_searching[iray_inds]
    ray_collided = np.array([True if n in irays else False for n in range(naz_mask)])
    ray_below_horizon = el_ray <= 0
    ray_found_horizon = (ray_searching & ray_prev_collided & ~ray_collided) \
                      | (ray_searching & ~ray_prev_collided & ray_collided)
    el_horizon[ray_found_horizon] = el_ray[ray_found_horizon]
    el_horizon[ray_below_horizon] = 0.0
    print(f'Ray tracing iteration {i}, {np.sum([True if n != -1.0 else False for n in el_horizon])}/{naz_mask} horizon points found!')
    i += 1
    ray_prev_collided = ray_collided
    el_ray[ray_collided & ray_searching] += d_el # Ones that collided need to go up
    el_ray[~ray_collided & ray_searching & ~ray_below_horizon] -= d_el # Ones that didn't collide need to go down
    ps.toc()


mask_pts = np.vstack(ps.sph_to_cart(az_ray, el_horizon)).T

mask_info = {'az': az_ray.tolist(), 'el': el_horizon.tolist(),
             'lat': station_lat_rad, 'lon': station_lon_rad,
             'station_terrain_height_km': station_height_above_terrain_km,
             'station_name': 'Purdue Optical Ground Station'}
with open('pogs_2.aem', 'w') as f:
    json.dump(mask_info, f)


pl = pv.Plotter()
pl.add_mesh(dem, opacity=0.5)
ps.scatter3(pl, mask_pts)
pl.camera.focal_point = [0.0,0.0,0.0]
# ps.plot_basis(pl, np.eye(3), labels=["E", "N", "U"])
pl.show()

sns.scatterplot(x=az_ray, y=el_ray)
plt.show()