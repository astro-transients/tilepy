import warnings

import matplotlib.cbook
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap

warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)


fig = plt.figure(num=None, figsize=(12, 8))
m = Basemap(projection="moll", lon_0=0, resolution="c")
m.drawcoastlines()
# m.fillcontinents(color='tan',lake_color='lightblue')
# draw parallels and meridians.
m.drawparallels(
    np.arange(-90.0, 91.0, 30.0), labels=[True, True, False, False], dashes=[2, 2]
)
m.drawmeridians(
    np.arange(-180.0, 181.0, 60.0), labels=[False, False, False, False], dashes=[2, 2]
)
m.drawmapboundary(fill_color="white")
plt.show()
