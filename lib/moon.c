/**
 * Location of the Moon.
 *
 * Calculation of the (Earth-fixed) latitude and longitude of the Moon.
 *
 * Equations are from Jean Meeus, Astronomical Algorithms, 1st ed.
 *
 * @author magne
 */
#include <stdio.h>
#include <math.h>

double deg2rad (const double d)
{
  return (d * M_PI) / 180.0;
}

/**
 * Calculate Julian century.
 *
 * @see AA eq. 22.1
 *
 * @param jde Julian Ephemeris Day
 * @return Julian century
 */
double  get_julian_century (const double jde)
{
  return (jde - 2451545) / 36525.0;
}

/**
 * Mean anomaly of the Sun.
 *
 * @see AA, page 308
 *
 * @param T time in Julian centuries
 * @return mean anomaly of the Sun
 */
double mean_anomaly_of_sun (const double T)
{
  double M;
  M = 
    357.5291092
    + (35999.0502909 * T)
    - (0.0001536 * T * T)
    - ((T * T * T) / 24490000.0);

  return M;
}

/**
 * Argument of latitude for the Moon.
 *
 * @param T time in julian centuries
 * @return Moons argument of latitude
 */
double moons_argument_of_latitude (const double T)
{
  double F;
  F = 
    93.2720993
    + (483202.0175273 * T)
    - (0.0034029 * T * T)
    - ((T * T * T) / 3526000.0)
    + ((T * T * T * T) / 863310000.0);

  return F;
}

/**
 * Mean anomaly of the Moon.
 *
 * @see AA, page 121
 *
 * @param T time in Julian centuries
 * @return mean anomaly of the Moon
 */
double mean_anomaly_of_moon (const double T)
{
  double Mprime;
  Mprime =
    134.96298
    + (477198.867398 * T)
    + (0.0086972 * T * T)
    + ((T * T * T) / 56250.0);

  return Mprime;
  
}

/**
 * Mean elongation of the Moon from the Sun.
 *
 * @see AA, page 308
 *
 * @param T time in Julian centuries
 * @return mean elongation
 */
double mean_elongation_of_moon_from_sun (const double T)
{
  double D;
  D =
    297.85022042
    + (445267.1115168 * T)
    - (0.0016300 * T * T)
    + ((T * T * T) / 545868.0)
    - ((T * T * T * T) / 113065000.0);

  return D;
}


double mean_longitude_of_moon (const double T)
{
  double T2 = (T * T);
  double T3 = (T * T2);
  double T4 = (T2 * T2);
  double Lp = 
      218.3164591 
    + 481267.88134236 * T
    - 0.0013268 * T2
    + (T3 / 538841.0)
    - (T4 / 65194000);

  return Lp;
}

double getEcorrection (const double T)
{
  return 1.0 - 0.002516 * T - 0.0000074 * T * T;
}


/**
 * Periodic terms for the longitude of the Moon.
 *
 * @see AA, page 308-311
 * @param T time, in Julian centuries
 * @param D mean elongation of the Moon
 * @param M Sun's mean anomaly
 * @param Mp Moon's mean anomaly
 * @param F  Moon's argument of latitude
 * @return longitude contribution by the periodic terms
 */
double periodic_terms_for_the_longitude (double T, double D, double M, double Mp, double F)
{
  double m  = deg2rad(M);
  double mp = deg2rad(Mp);
  double d  = deg2rad(D);
  double f  = deg2rad(F);

  /* 
     E is a correction to the M terms due to the decreasing eccentricity
     of the Earth's orbit. Each M term is multiplied by E, and each
     2M term is multiplied  by E^2
  */
  double E = 1.0 - 0.002516 * T - 0.0000074 * T * T;
  double E2 = E * E;

  /* periodic terms from AA, table 45.A, page 309 */
  double suml =
      6.288774 * sin (mp)
    + 1.274027 * sin (2.0 * d - mp)
    + 0.658314 * sin (2.0 * d)
    + 0.213618 * sin (2.0 * mp)
    - 0.185116 * sin (m) * E
    - 0.114332 * sin (2.0 * f)
    + 0.058793 * sin (2.0 * d - 2.0 * mp)
    + 0.057066 * sin (2.0 * d - m - mp) * E
    + 0.053322 * sin (2.0 * d + mp)
    + 0.045758 * sin (2.0 * d - m) * E
    - 0.040923 * sin (2.0 * m - mp) * E
    - 0.034720 * sin (2.0 * d)
    - 0.030383 * sin (m - mp) * E
    + 0.015327 * sin (2.0 * d - 2.0 * f)
    - 0.012528 * sin (mp + 2.0 * f)
    + 0.010980 * sin (mp - 2.0 * f)
    + 0.010675 * sin (4.0 * d - mp)
    + 0.010034 * sin (3.0 * mp)
    + 0.008548 * sin (4.0 * d - 2.0 * mp)
    - 0.007888 * sin (2.0 * d + m - mp) * E
    - 0.006766 * sin (2.0 * d + m) * E
    - 0.005163 * sin (d - mp)
    + 0.004987 * sin (d + m) * E
    + 0.003994 * sin (2.0 * d + 2.0 * mp)
    + 0.003861 * sin (4.0 * d)
    + 0.003665 * sin (2.0 * d - 3.0 * mp)
    - 0.002689 * sin (m - 2.0 * mp) * E
    - 0.002602 * sin (2.0 * d - mp)
    + 0.002390 * sin (2.0 * d - m - mp + 2.0 * f) * E
    - 0.002348 * sin (d - mp)
    + 0.002236 * sin (2.0 * d - 2.0 * m) * E2
    - 0.002120 * sin (m + 2.0 * mp) * E
    - 0.002069 * sin (2.0 * m) * E2
    + 0.002048 * sin (2.0 * d - 2.0 * m - mp) * E2
    - 0.001773 * sin (2.0 * d + m - 2.0 * f)
    - 0.001595 * sin (2.0 * d + 2.0 * f)
    + 0.001215 * sin (4.0 * d - m - mp) * E
    - 0.001110 * sin (2.0 * mp + 2.0 * f)
    - 0.000892 * sin (3.0 * d - mp)
    - 0.000810 * sin (2.0 + m + mp) * E
    + 0.000759 * sin (4.0 * d - m - 2.0 * mp) * E
    - 0.000713 * sin (2.0 * m - mp) * E2
    - 0.000700 * sin (2.0 * d + 2.0 * m - mp) * E2
    + 0.000691 * sin (2.0 * d + m - 2.0 * mp) * E
    + 0.000596 * sin (2.0 * d - m - 2.0 * f) * E
    + 0.000549 * sin (4.0 * d + mp)
    + 0.000537 * sin (4.0 * mp)
    + 0.000520 * sin (4.0 * d - m) * E
    - 0.000487 * sin (d - 2.0 * mp)
    - 0.000399 * sin (2.0 * d + m - 2.0 * f) * E
    - 0.000381 * sin (2.0 * mp  - 2.0 * f)
    + 0.000351 * sin (d + m + mp) * E
    - 0.000340 * sin (3.0 * d - 2.0 * mp)
    + 0.000330 * sin (4.0 * d - 3.0 * mp)
    + 0.000327 * sin (2.0 * d - m + 2.0 * mp) * E
    - 0.000323 * sin (2.0 * m + mp) * E2
    + 0.000299 * sin (d + m - mp) * E
    + 0.000294 * sin (2.0 * d - mp - 3.0 * f)
    ;

  return suml;
}

/**
 * Periodic terms for the latitude of the Moon.
 *
 * @see AA, page 308-311
 * @param T time, in Julian centuries
 * @param D mean elongation of the Moon
 * @param M Sun's mean anomaly
 * @param Mp Moon's mean anomaly
 * @param F  Moon's argument of latitude
 * @return latitude contribution by the periodic terms
 */
double periodic_terms_for_the_latitude (double T, double D, double M, double Mp, double F)
{
  double m  = deg2rad(M);
  double mp = deg2rad(Mp);
  double d  = deg2rad(D);
  double f  = deg2rad(F);

  /* 
     E is a correction to the M terms due to the decreasing eccentricity
     of the Earth's orbit. Each M term is multiplied by E, and each
     2M term is multiplied  by E^2
  */
  double E = 1.0 - 0.002516 * T - 0.0000074 * T * T;
  double E2 = E * E;

  /* periodic terms from AA, table 45.B, page 311 */
  double sumb =
      5.128122 * sin (f)
    + 0.280602 * sin (mp + f)
    + 0.277693 * sin (mp - f)
    + 0.173237 * sin (2.0 * d - f)
    + 0.055413 * sin (2.0 * d - mp + f)
    + 0.046271 * sin (2.0 * d - mp - f)
    + 0.032573 * sin (2.0 * d + f)
    + 0.017198 * sin (2.0 * mp + f)
    + 0.009266 * sin (2.0 * d + mp - f)
    + 0.008822 * sin (2.0 * mp - f)
    + 0.008216 * sin (2.0 * d - m - f) * E
    + 0.004324 * sin (2.0 * d - 2.0 * mp - f)
    /* TODO(magne): Add remaining terms */
    ;
  return sumb;
}

double getA1correction (const double T)
{
  return 119.75 + 131.849 * T;
}

double getA2correction (const double T)
{
  return 53.09 + 479264.290 * T;
}

double getA3correction (const double T)
{
  return 313.45 + 481266.484 * T;
}

/**
 * Latitude of the Moon.
 *
 * @see AA, page 312
 * @param T time, in Julian centuries
 * @return latitude of the Moon
 */
double latitude_of_moon (const double T)
{
  double Lp = mean_longitude_of_moon (T); 
  double D  = mean_elongation_of_moon_from_sun (T);
  double M  = mean_anomaly_of_sun (T);
  double Mp = mean_anomaly_of_moon (T);
  double F  = moons_argument_of_latitude (T);
  double A1 = getA1correction (T);
  double A2 = getA2correction (T);

  double suml = periodic_terms_for_the_latitude (T, D, M, Mp, F);
  /* Additives due to Venus, Jupiter and the flattening of Earth */
  suml +=
      0.003958 * sin (deg2rad(A1))     /* Venus */
    + 0.001962 * sin (deg2rad(Lp - F)) /* Earth's flattening */
    + 0.000318 * sin (deg2rad(A2))     /* Jupiter */
    ;

  double lat = Lp + suml;
  return lat;
}


/**
 * Latitude of the Moon.
 *
 * @see AA, page 312
 * @param T time, in Julian centuries
 * @return latitude of the Moon
 */
double longitude_of_moon (const double T)
{
  double Lp = mean_longitude_of_moon (T); 
  double D  = mean_elongation_of_moon_from_sun (T);
  double M  = mean_anomaly_of_sun (T);
  double Mp = mean_anomaly_of_moon (T);
  double F  = moons_argument_of_latitude (T);

  double A1 = getA1correction (T);
  double A2 = getA2correction (T);
  double A3 = getA3correction (T);

  double sumb = periodic_terms_for_the_longitude (T, D, M, Mp, F);
  /* Additives due to Venus, Jupiter and the flattening of Earth */
  sumb +=
    - 0.002235 * sin (deg2rad(Lp))
    + 0.000382 * sin (deg2rad(A3))
    + 0.000175 * sin (deg2rad(A1 - F))
    + 0.000175 * sin (deg2rad(A2 + F))
    + 0.000127 * sin (deg2rad(Lp - Mp))
    - 0.000115 * sin (deg2rad(Lp + Mp));

  return sumb;
}

/**
 * Apparent longitude of the Moon.
 *
 * The apparent longitude of the Moon is the real
 * longitude corrected for the nutation in longitude.
 *
 * @see AA, page 312
 * @param T time, in Julian centuries
 * @return apparent longitude of moon
 */


double angle (double a)
{
  while (a < 0.0)
    a += 360.0;

  while (a >= 360.0)
    a -= 360.0;

  return a;
}


static void test1 ()
{
  const double jde = 2448724.5;

  const double T = get_julian_century (jde);
  
  printf("JDE  = %lf\n", jde);
  printf("T    = %lf\n", T);
  printf("L'   = %lf\n", angle(mean_longitude_of_moon (T)));
  printf("D    = %lf\n", angle(mean_elongation_of_moon_from_sun (T)));
  printf("M    = %lf\n", angle(mean_anomaly_of_sun (T)));
  printf("M'   = %lf\n", angle(mean_anomaly_of_moon (T)));
  printf("F    = %lf\n", angle(moons_argument_of_latitude (T)));
  printf("A1   = %lf\n", angle(getA1correction (T)));
  printf("A2   = %lf\n", angle(getA2correction (T)));
  printf("A3   = %lf\n", angle(getA3correction (T)));
  printf("E    = %lf\n", getEcorrection (T));
  printf("\n");
  printf("lat  = %lf\n", angle(latitude_of_moon (T)));
  printf("lon  = %lf\n", angle(longitude_of_moon (T)));
  return;
}

int main()
{
  test1();
}
