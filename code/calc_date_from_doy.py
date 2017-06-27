from astropy.time import Time
def calc_date_from_doy(doy, year):
    beginning = Time('{}-12-31 00:00:00'.format(year-1), format='iso')
    new_date = Time(beginning.mjd+doy, format='mjd')
    return new_date.iso