import datetime

import matplotlib.pyplot as plt
from terrainman import TsiDataHandler
import mirage as mr
import numpy as np


import terrainman as tm
tm.purge_data()

if __name__ == "__main__":
    date = mr.utc(2000, 12, 9)
    dates = mr.date_linspace(date, mr.now(), 10_000)
    tsi_dh = TsiDataHandler()
    mr.tic()
    sc = tsi_dh.eval(dates)
    mr.toc()
    print(sc)


    plt.scatter(mr.date_to_jd(dates), sc, s=1)
    plt.title("Total Solar Irradiance")
    plt.xlabel("Julian date [days]")
    plt.ylabel(r"Irradiance at 1 AU $\left[\frac{W}{m^2}\right]$")
    plt.show()