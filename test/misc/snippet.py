import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def main():
    bm = Basemap(projection="rotpole",
                 o_lat_p=36.0,
                 o_lon_p=180.0,
                 llcrnrlat=-10.590603,
                 urcrnrlat=46.591976,
                 llcrnrlon=-139.08585,
                 urcrnrlon=22.661009,
                 lon_0=-106.0,
                 rsphere=6370000,
                 resolution='l')

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    bm.drawcoastlines(linewidth=.5)

    print(bm.proj4string)

    plt.savefig("basemap_map.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
