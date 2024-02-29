import requests
from bs4 import BeautifulSoup

# list all linked files here
# https://www.ncei.noaa.gov/data/total-solar-irradiance/access/daily/

import mirage as mr

mr.tic()
page = requests.get("https://www.ncei.noaa.gov/data/total-solar-irradiance/access/daily/")
soup = BeautifulSoup(page.content, 'html.parser')
# get all hrefs
for link in soup.find_all('a'):
    if 'tsi_' in link.get('href'):
        # print(link.get('href'))
        pass
mr.toc()