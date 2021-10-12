"""
# script for downloading IS2 data with ICEPYX
# built from example at: https://icepyx.readthedocs.io/en/latest/getting_started/example_link.html
"""

import icepyx as ipx

def download_atl03(short_name, spatial_extent, date_range, earthdata_uid, earthdata_email, page_size, page_num, out_path):

    """ Function to search and download IS2"""
    region_a = ipx.Query(short_name, spatial_extent, date_range)

    print('\n print avail granules = ', region_a.avail_granules(), '\n')

    cont = input("Continue...? (yes/no): ")
    if cont == "yes":
        session=region_a.earthdata_login(earthdata_uid, earthdata_email)
        region_a.download_granules(session, out_path, page_size, page_num)
    elif cont == "no":
        print("answer 'yes' to continue")
    else:
        print('Answer "yes" or "no"')

def main():
    """ run functions """
    #short_name = the dataset of interest, known as its "short name".
    #See https://nsidc.org/data/icesat-2/data-sets for a list of the available datasets.
    #spatial extent = a region of interest to search within. This can be entered as a bounding box,
    #polygon vertex coordinate pairs, or a polygon geospatial file (currently shp, kml,
    #and gpkg are supported).
    #bounding box: Given in decimal degrees for the lower left longitude,
    #lower left latitude, upper right longitude, and upper right latitude
    #polygon vertices: Given as longitude, latitude coordinate pairs
    #of decimal degrees with the last entry a repeat of the first.
    #polygon file: A string containing the full file path and name.
    #date_range = the date range for which you would like to search for results.
    #Must be formatted as a set of 'YYYY-MM-DD' strings.

    short_name = 'ATL03'
    spatial_extent = [-82.14258, 21.60278, -81.47461, 22.18276]
    #spatial_extent = './supporting_files/data-access_PineIsland/glims_polygons.kml'
    date_range = ['2021-02-05','2021-02-07']
    earthdata_uid = 'bhylee'
    earthdata_email = 'brianlee52@ucsb.edu'
    page_size = 10
    page_num = 3
    out_path = '/Users/brian.h.lee/Desktop/icesat/data/raw'

    download_atl03(short_name, spatial_extent, date_range, earthdata_uid, earthdata_email, page_size, page_num, out_path)

if __name__=='__main__':
    main()
