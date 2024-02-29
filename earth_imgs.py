import datetime
import pyspaceaware as ps
import numpy as np

import matplotlib.pyplot as plt

from ftplib import FTP
import os
import threading
from PIL import Image
from matplotlib.animation import FuncAnimation
import time
import pyvista as pv

_PLANET_RES = 500
def textured_ellipsoid(
    a: float, f: float, rotm: np.ndarray = np.eye(3)
) -> pv.ParametricEllipsoid:
    """Returns a parametric ellipsoid based on the given semimajor axis and flattening ratio, ready for texturing

    :param a: Semimajor axis of the body
    :type a: float
    :param f: Flattening ratio of the body
    :type f: float
    :param rotm: Rotation matrix to apply to the body points, defaults to np.eye(3)
    :type rotm: np.ndarray [3x3]
    :return: Ellipsoid parametrized by ``a``, ``f``
    :rtype: pv.ParametricEllipsoid
    """
    b = a * (1 - f)  # Semi-minor axis
    body = pv.ParametricEllipsoid(
        a,
        a,
        b,
        u_res=_PLANET_RES,
        v_res=_PLANET_RES,
        w_res=_PLANET_RES,
        min_u=-np.pi / 2 + 0.00001,
        max_u=3 * np.pi / 2 - 0.00001,
    )
    sph_hat_pts = ps.hat(body.points)
    body.active_t_coords = np.zeros((body.points.shape[0], 2))
    body.active_t_coords[:, 0] = 0.5 + np.arctan2(
        -sph_hat_pts[:, 0], sph_hat_pts[:, 1]
    ) / (2 * np.pi)
    body.active_t_coords[:, 1] = 0.5 + np.arcsin(sph_hat_pts[:, 2]) / np.pi
    body.points = (rotm @ body.points.T).T
    return body


_IMG_DIR = os.path.join(os.curdir, 'earth_daily_color')
if not os.path.exists(_IMG_DIR):
    os.makedirs(_IMG_DIR)

def dates_to_udates(dates: np.ndarray[datetime.datetime]) -> np.ndarray[datetime.datetime]:
    dates = [datetime.datetime(year=d.year, month=d.month, day=d.day, tzinfo=datetime.timezone.utc) for d in dates]
    return np.array(sorted(list(set(dates))))

def download_date_range_color(dates: np.ndarray[datetime.datetime]):
    udates = dates_to_udates(dates)
    threads = []

    for date in udates:
        thread = threading.Thread(target=download_day_color, args=(date,))
        threads.append(thread)

    # Start all the threads
    for thread in threads:
        thread.start()

    # Wait for all threads to finish
    for thread in threads:
        thread.join()


def download_day_color(date: datetime.datetime):
    fname = f'TRUE.daily.{date.strftime("%Y%m%d")}.color.png'
    fpath = os.path.join(_IMG_DIR, fname)
    if not os.path.exists(fpath):
        time.sleep(np.random.rand())
        ftp_server = 'ftp.nnvl.noaa.gov'
        ftp = FTP(host=ftp_server)
        ftp.login()  # Logs in as anonymous/guest
        ftp.cwd('/View/TRUE/Images/Color/Daily/')
        with open(fpath, 'wb') as local_file:
            ftp.retrbinary('RETR ' + fname, local_file.write)
        ftp.quit()

if __name__ == "__main__":
    date = datetime.datetime(2023, 6, 10, tzinfo=datetime.timezone.utc)
    dates, _ = ps.date_linspace(date - ps.days(200), date, 10_00)
    udates = dates_to_udates(dates)

    # ps.tic()
    # download_date_range_color(dates)
    # ps.toc()

    fnames = sorted([x for x in os.listdir(_IMG_DIR) if x.endswith('.color.png')])
    # fig, ax = plt.subplots()
    # image_path = os.path.join(_IMG_DIR, fnames[0])
    # im = ax.imshow(Image.open(image_path))

    # def update(i):
    #     image_path = os.path.join(_IMG_DIR, fnames[i])
    #     img = Image.open(image_path)
    #     im.set_data(img)
    #     ax.set_title(udates[i].strftime('%Y-%m-%d'))
    #     print(fnames[i])
    #     return im,

    # animation = FuncAnimation(fig, update, frames=len(fnames), interval=50)

    # # Save the animation as a GIF
    # animation_filename = 'animation.gif'
    # animation.save(animation_filename, writer='pillow')

    pl = pv.Plotter()
    pl.enable_anti_aliasing('ssaa')
    body = textured_ellipsoid(ps.AstroConstants.earth_r_eq, ps.AstroConstants.earth_f)
    pl.set_background('black')
    pl.open_gif("wave_day.gif", fps=10)

    dates, _ = ps.date_linspace(udates[0], udates[4], 40)

    for date in dates:
        today = datetime.datetime(year=date.year, month=date.month, day=date.day, tzinfo=datetime.timezone.utc)
        today_ind = np.argmin(np.abs([(d - today).total_seconds() for d in udates]))
        print(date, today_ind)
        date_frac = (date - today).total_seconds() / 86400
        prev_file = os.path.join(_IMG_DIR, fnames[today_ind])
        next_file = os.path.join(_IMG_DIR, fnames[today_ind+1])
        prev_im = np.array(Image.open(prev_file))
        next_im = np.array(Image.open(next_file))
        prev_im[:,:,-1] = 255
        next_im[:,:,-1] = 255
        im = prev_im * (1 - date_frac) + date_frac * next_im

        tex = pv.numpy_to_texture(im)
        pl.clear()
        pl.add_mesh(body, texture=tex)
        pl.write_frame()

    pl.close()