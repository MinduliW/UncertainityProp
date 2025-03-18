#ifndef dynamics_H
#define dynamics_H

template <typename U>
U Density_HP(U heightm)
{

  U h[50] = {100.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0, 720.0, 740.0, 760.0, 780.0, 800.0, 840.0, 880.0, 920.0, 960.0, 1000.0};

  U c_min[50] = {4.974e+05, 2.490e+04, 8.377e+03, 3.899e+03, 2.122e+03, 1.263e+03,
                 8.008e+02, 5.283e+02, 3.617e+02, 2.557e+02, 1.839e+02, 1.341e+02,
                 9.949e+01, 7.488e+01, 5.709e+01, 4.403e+01, 3.430e+01, 2.697e+01,
                 2.139e+01, 1.708e+01, 1.099e+01, 7.214e+00, 4.824e+00, 3.274e+00,
                 2.249e+00, 1.558e+00, 1.091e+00, 7.701e-01, 5.474e-01, 3.916e-01,
                 2.819e-01, 2.042e-01, 1.488e-01, 1.092e-01, 8.070e-02, 6.012e-02,
                 4.519e-02, 3.430e-02, 2.632e-02, 2.043e-02, 1.607e-02, 1.281e-02,
                 1.036e-02, 8.496e-03, 7.069e-03, 4.680e-03, 3.200e-03, 2.210e-03,
                 1.560e-03, 1.150e-03};

  heightm = heightm / 1000.0;

  U density = 0;

  if (heightm > 1000.0)
  {
    density = 1.150e-03;
  }
  else if (heightm < 100.0)
  {
    density = 4.974e+05;
  }
  else
  {

    for (int i = 0; i < 50; i++)
    {
      if (h[i] <= heightm && h[i + 1] >= heightm)
      {
        U diffh = heightm - h[i];
        U diffh2 = h[i + 1] - h[i];

        density = c_min[i] + (c_min[i + 1] - c_min[i]) * diffh / diffh2;
      }
    }
  }

  return density * 1e-12;
}

template <typename U>
U Smoothed_eclipse(U cs, U ct, DACE::AlgebraicVector<double> x_sun, DACE::AlgebraicVector<double> x_sat, U Rs, U Re)
{

  using DACE::cos;
  using DACE::pow;
  using DACE::sin;
  using DACE::sqrt;
  using std::cos;
  using std::pow;
  using std::sin;
  using std::sqrt;

  // cout << "Rs =" << Rs << endl;

  U a_SR = asin(Rs / x_sun.vnorm());
  U a_BR = asin(Re / x_sat.vnorm());

  U a_D = acos(x_sat.dot(x_sun) / x_sat.vnorm() / x_sun.vnorm());

  U gamma_L = 1 / (1 + exp(-cs * (a_D - ct * (a_SR + a_BR))));

  return gamma_L;
}

template <typename U>
DACE::AlgebraicVector<U> getSunPosVec(U JD, U LU)
{

  U pi = 4 * atan(1.0);

  U T_UT1 = (JD - 2451545.0) / 36525;

  U lambdaMSun = 280.460 + 36000.771 * T_UT1;

  U MSun = 357.5291092 + 35999.05034 * T_UT1;

  U lambdaEcliptic = lambdaMSun + 1.914666471 * sin(MSun * pi / 180.0) + 0.019994643 * sin(2 * MSun * pi / 180.0);

  U r_Sun = 1.000140612 - 0.016708617 * cos(MSun * pi / 180.0) - 0.000139589 * cos(2 * MSun * pi / 180.0);

  U epsilon = 23.439291 - 0.0130042 * T_UT1;

  DACE::AlgebraicVector<U> SunVec(3);
  SunVec[0] = r_Sun * cos(lambdaEcliptic * pi / 180.0);
  SunVec[1] = r_Sun * cos(epsilon * pi / 180.0) * sin(lambdaEcliptic * pi / 180.0);
  SunVec[2] = r_Sun * sin(epsilon * pi / 180.0) * sin(lambdaEcliptic * pi / 180.0);

  SunVec = 1.495978707e11 * SunVec / LU;

  // cout << SunVec << endl;
  return SunVec;
}

template <typename T, typename U>
double getEclipse(U t, DACE::AlgebraicVector<T> x, U cs, U ct, U LU, U Rs, U Re)
{

  DACE::AlgebraicVector<double> r_Earth2Sun(3);

  r_Earth2Sun = getSunPosVec(t + 2400000.5 + 60000, LU);

  // convert dace to double.
  DACE::AlgebraicVector<double> x_double(6);
  for (unsigned int i = 0; i < 6; i++)
  {
    x_double[i] = DACE::cons(x[i]);
  }

  DACE::AlgebraicVector<double> CART = x_double;

  DACE::AlgebraicVector<double> rSat2Earth(3);
  DACE::AlgebraicVector<double> rSat2Sun(3);

  for (unsigned int i = 0; i < 3; i++)
  {
    rSat2Earth[i] = -CART[i];
    rSat2Sun[i] = rSat2Earth[i] + r_Earth2Sun[i];
  }

  double numinus1 = Smoothed_eclipse(cs, ct, rSat2Sun, rSat2Earth, Rs, Re);

  return 1.0 - numinus1;
}

template <typename T, typename U>
DACE::AlgebraicVector<T> DynamicsCart(DACE::AlgebraicVector<T> x, U t, const double *cst)
{

  using DACE::cos;
  using DACE::pow;
  using DACE::sin;
  using DACE::sqrt;
  using std::cos;
  using std::pow;
  using std::sin;
  using std::sqrt;

  const double mu = cst[5];
  const double J2 = cst[6];
  const double Re = cst[7];
  const double LU = cst[9];
  const double Area = cst[10];
  const double Cd = cst[11];
  const double m0 = cst[4];
  const double massunit = cst[12];

  DACE::AlgebraicVector<T> dxdt(6);

  T r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  // cout << DACE::cons(term) << endl;

  if (DACE::cons(r2) <= 0)
  {
    cout << "So the iter.tex has a mistake, its giving all states as zero, likely after infeasible run" << endl;
  }

  T r = sqrt(r2);

  T r3 = r * r * r;

  for (unsigned int i = 0; i < 3; i++)
  {
    dxdt[i] = x[3 + i];
    dxdt[i + 3] = -mu * x[i] / (r3);
  }

   // cout << J2 << endl;

  dxdt[3] = dxdt[3] - mu * x[0] / (r3)*1.5 * J2 * (Re / r) * (Re / r) * (1 - 5 * x[2] * x[2] / r2);
  dxdt[4] = dxdt[4] - mu * x[1] / (r3)*1.5 * J2 * (Re / r) * (Re / r) * (1 - 5 * x[2] * x[2] / r2);
  dxdt[5] = dxdt[5] - mu * x[2] / (r3)*1.5 * J2 * (Re / r) * (Re / r) * (3 - 5 * x[2] * x[2] / r2);

  // ADD DRAG

//     if (cst[13] > 0)
//     {
// 
//         double rho = Density_HP(DACE::cons(r * LU - Re * LU))*LU*LU*LU/massunit;
// 
//         T vmag = sqrt(x[3] * x[3] + x[4] * x[4] + x[5] * x[5]);
//         for (unsigned int i = 0; i < 3; i++)
//         {
//             dxdt[3 + i] += 0.5 * rho * Cd * Area * vmag * vmag / m0 * -x[3 + i] / vmag;
//         }
//     }
// 
//   if (cst[8] < 2)
//   {
//     DACE::AlgebraicVector<T> u(3);
//     for (unsigned int i = 0; i < 3; i++)
//     {
//       // cout << params.u0[i] << endl;
//       u[i] = cst[i] + DACE::DA(i + 7); //
//     }
// 
//     dxdt[3] += u[0];
//     dxdt[4] += u[1];
//     dxdt[5] += u[2];
// 
//     // cout << u[2] << endl;
//   }

  return dxdt;
}

#endif