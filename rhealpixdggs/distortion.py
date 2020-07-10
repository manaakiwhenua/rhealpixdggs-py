"""
This Python 3.3 module computes linear and areal distortion statistics of 
map projections.
It is used for analysis only and so is not essential to manipulating 
the rHEALPix discrete global grid system.

CHANGELOG:

- Alexander Raichev (AR), 2012-10-01: Initial version.
- AR, 2013-01-21: Wasn't working in degrees mode. Fixed that bug.
- AR, 2013-07-23: Ported to Python 3.3.

"""
# *****************************************************************************
#       Copyright (C) 2012 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http: //www.gnu.org/licenses/
# *****************************************************************************
# Import third-party modules.
from numpy import pi, floor, sqrt, sin, cos, arcsin, deg2rad, rad2deg, array, median

# Import my modules.
from rhealpixdggs.utils import my_round, wrap_longitude, wrap_latitude


def fff_coeffs(T, lam, phi):
    """
    Return numerical approximations of the `first fundamental form <http://en.wikipedia.org/wiki/First_fundamental_form>`_ coefficients
    E, F, and G (in that order) of map projection `T` 
    (an projection_tools.Proj or projection_tools.Proj4 instance) at longitude `lam` and 
    latitude `phi`.
    
    EXAMPLES::
    
        >>> from rhealpixdggs.projection_wrapper import Proj
        >>> from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID_RADIANS
        >>> f = Proj(ellipsoid=WGS84_ELLIPSOID_RADIANS, proj='healpix')
        >>> print(my_round(fff_coeffs(f, 0, pi/6), 15))
        (40635288880650.484, 0.0, 42251277118793.328)

    """
    dellam_rad = 1e-6
    delphi_rad = 1e-6
    radians = T.ellipsoid.radians
    if not radians:
        # Convert to degrees.
        dellam, delphi = rad2deg([dellam_rad, delphi_rad])
    else:
        dellam, delphi = dellam_rad, delphi_rad
    # Standardize new lon and lat to avoid out of bounds errors.
    lam2 = wrap_longitude(lam + dellam, radians=radians)
    phi2 = wrap_latitude(phi + delphi, radians=radians)
    delxdellam, delydellam = (array(T(lam2, phi)) - array(T(lam, phi))) / dellam_rad
    delxdelphi, delydelphi = (array(T(lam, phi2)) - array(T(lam, phi))) / delphi_rad
    E = delxdellam ** 2 + delydellam ** 2
    F = delxdellam * delxdelphi + delydellam * delydelphi
    G = delxdelphi ** 2 + delydelphi ** 2
    # if G == 0:
    #     print 'The projection is not injective near:', lam, phi
    # if E == 0:
    #     print 'E = 0 at', lam, phi
    return E, F, G


def scale_factors(T, lam, phi):
    """
    Return numerical approximations of the local scale factors s_M, s_P, and 
    s_A of the map projection `T` 
    (an projection_tools.Proj or projection_tools.Proj4 instance)
    at longitude `lam` and  latitude `phi`, where 
    s_M is the local linear scale along meridians, 
    s_P is the local linear scale along parallels, and
    s_A is the local area scale. 
    Also return theta (in radians), the angle between the vectors (delxdellam, delydellam)
    and (delxdelphi, delydelphi).

    OUTPUT:

    (s_M, s_P, s_A, theta)
    
    EXAMPLES::
    
        >>> from rhealpixdggs.projection_wrapper import Proj
        >>> from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID_RADIANS
        >>> f = Proj(ellipsoid=WGS84_ELLIPSOID_RADIANS, proj='healpix')
        >>> print(my_round(scale_factors(f, 0, pi/6), 15))
        (1.0212575853790069, 1.1535746974071359, 1.1780969100283301, 1.5707963267948959)


    """
    # Set length of a radian of longitude on the sphere/ellipsoid and
    # length of a radian of latitude on the sphere/ellipsoid.
    ellps = T.ellipsoid
    if not ellps.radians:
        # Convert to radians.
        phi_rad = deg2rad(phi)
    else:
        phi_rad = phi
    if ellps.sphere:
        R = ellps.R
        L_lam = R * cos(phi_rad)
        L_phi = R
    else:
        a = ellps.a
        e = ellps.e
        L_lam = a * cos(phi_rad) / sqrt(1 - e ** 2 * sin(phi_rad) ** 2)
        L_phi = a * (1 - e ** 2) / sqrt(1 - e ** 2 * sin(phi_rad) ** 2) ** 3
    E, F, G = fff_coeffs(T, lam, phi)
    H = sqrt(E * G - F ** 2)
    s_M = sqrt(G) / L_phi
    s_P = sqrt(E) / L_lam
    s_A = H / (L_lam * L_phi)
    theta = arcsin(H / sqrt(E * G))
    return s_M, s_P, s_A, theta


def distortion(T, lam, phi):
    """
    Return ((x, y), mad, ld, ad)
    Here `(x, y) = T(lam, phi)`, the image under the map projection `T` 
    (an projection_tools.Proj or projection_tools.Proj4 instance)
    of the longitude-latitude point `(lam, phi)`;
    mad = maximum angular distortion at (x, y) = 2*arcsin((A - B)/(A + B))
    ld = linear distortion at (x, y) = A/B;
    ad = areal distortion at (x, y) = AB;
    A and B are the major and minor radii, respectively, of the Tissot ellipse 
    at (x, y).
    
    EXAMPLES::
    
        >>> from rhealpixdggs.projection_wrapper import Proj
        >>> from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID_RADIANS
        >>> f = Proj(ellipsoid=WGS84_ELLIPSOID_RADIANS, proj='healpix')
        >>> print(my_round(distortion(f, 0, pi/6), 15))
        ((0.0, 3748655.1150495014), 0.121755482930707, 1.1295629172526771, 1.1780969100283301)

    """
    x, y = T(lam, phi)
    s_M, s_P, s_A, theta = scale_factors(T, lam, phi)
    AA = sqrt(s_M ** 2 + s_P ** 2 + 2 * s_M * s_P * sin(theta))
    if s_M ** 2 + s_P ** 2 - 2 * s_M * s_P * sin(theta) < 0:
        # This case is analytically impossible but numerically possible,
        # because of my_rounding errors.
        BB = 0
    else:
        BB = sqrt(s_M ** 2 + s_P ** 2 - 2 * s_M * s_P * sin(theta))
    A = (AA + BB) / 2.0
    B = (AA - BB) / 2.0
    mad = 2 * arcsin((A - B) / (A + B))
    if not T.ellipsoid.radians:
        # Convert to degrees.
        mad = rad2deg(mad)
    ld = A / B
    ad = A * B
    return (x, y), mad, ld, ad


def distortion_stats(T, sample, my_round_numbers=3):
    """
    Return the sample minimum, sample maximum, sample median, sample mean, 
    and sample standard deviation of the maximum angular distortion, linear 
    distortion, and area distortion functions
    for the map projection `T` 
    (an projection_tools.Proj or projection_tools.Proj4 instance)
    of the list `sample` of longitude-latitude 
    points chosen from the surface `T.ellipsoid`.
    Most likely you will want sample to comprise points sampled uniformly 
    at random from the surface of `T.ellipsoid` (and not simply sampled 
    uniformly at random from the rectangle (-pi, pi) x (-pi/2, pi/2)).
        
    OUTPUT:
    
    (distortions, stats), where
    distortions is a list of distortion() outputs for each longitude-
    latitude point sampled;
    stats is the list of lists [maximum angular distortion stats, 
    linear distortion stats, area distortion stats], where each stats sublist 
    is of the form [sample mean, sample standard deviation, sample minimum,
    sample maximum, sample median].
    
    EXAMPLES::
    
        >>> from rhealpixdggs.projection_wrapper import Proj
        >>> from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID_RADIANS
        >>> E = WGS84_ELLIPSOID_RADIANS
        >>> f = Proj(ellipsoid=E, proj='healpix')
        >>> sample = [E.random_point() for i in range(100)]
        >>> print(distortion_stats(f, sample)[1])  # doctest: +SKIP 
        [[0.309, 0.238, 0.001, 0.838, 0.178], [1.41, 0.375, 1.001, 2.372, 1.195], [1.178, 0.0, 1.178, 1.178, 1.178]]

    """
    distortions = []
    distortions_slice = []
    k = len(distortion(T, 0, 0)) - 1
    sample_sum = array([0.0 for i in range(k)])
    sample_min = array([float("inf") for i in range(k)])
    sample_max = array([float("-inf") for i in range(k)])
    N = len(sample)
    for (lam, phi) in sample:
        d = distortion(T, lam, phi)
        distortions.append(d)
        distortions_slice.append(d[1:])
        for i in range(k):
            if d[1 + i] < sample_min[i]:
                sample_min[i] = d[1 + i]
            if d[1 + i] > sample_max[i]:
                sample_max[i] = d[1 + i]
        sample_sum += array(d[1:])
    distortions_slice = array(distortions_slice)
    sample_median = median(distortions_slice, axis=0)
    sample_mean = sample_sum / N
    sample_var = sum(((d - sample_mean) ** 2 for d in distortions_slice)) / (N - 1)
    sample_std = sqrt(sample_var)
    # sample_mean_error = sample_std/sqrt(N)
    # Zero until i find an appropriate estimate:
    # sample_std_error = [0 for i in range(k)]
    # Compile stats.
    stats = [sample_mean, sample_std, sample_min, sample_max, sample_median]
    # errors = [sample_mean_error, sample_std_error]
    # Reorder to have ld stuff followed by ad stuff.
    stats = list(zip(*stats))
    # errors = zip(*errors)
    ## Compile result
    # stats = [stats[0], errors[0], stats[1], errors[1]]
    # Round stats if desired.
    if my_round_numbers:
        new_stats = []
        for x in stats:
            new_stats.append([my_round(xx, my_round_numbers) for xx in x])
        stats = new_stats
    return distortions, stats


def utm_zone(lam, phi):
    """
    Return the Universal Transverse Mercator zone and hemisphere for longitude
    `lam` and latitude `phi` given in radians.
    Based on the WGS84 ellipsoid.
    Return None if `phi` is out of bounds, that is, if `phi` is greater than
    84 degrees or less than 80 degrees.
    
    EXAMPLES::
    
        >>> print(utm_zone(0, 84*pi/180))
        31
        >>> print(utm_zone(0, 85*pi/180)) 
        None

    """
    if phi > 84 * pi / 180 or phi < -80 * pi / 180:
        return None
    dlam = 6 * pi / 180
    zone = int(floor((lam + pi) / dlam) + 1)
    return zone
