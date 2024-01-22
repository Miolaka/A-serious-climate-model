// ##########################################################
// ############ functions ###################################
// ##########################################################
void sort_vertically(double Theta[], int nlay, int unstable)
{
    unstable = 1;
    while (unstable)
    {
        unstable = 0;
        for (int i = nlay - 2; i >= 0; i--)
        { // 9,8,7,6,5,...,0
            if (Theta[i] < Theta[i + 1])
            {
                unstable = 1;
                double tmp = Theta[i + 1];
                Theta[i + 1] = Theta[i];
                Theta[i] = tmp;
            }
        }
    }
}

void get_Theta_from_T(double Theta[], double T[], int nlay, double p[])
{
    for (int i = nlay - 1; i >= 0; i--)
    {
        Theta[i] = T[i] * pow((100000.0 / ((p[i] + p[i + 1]) / 2)), 2.0 / 7.0);
    }
}

void get_T_from_Theta(double Theta[], double T[], int nlay, double p[])
{
    for (int i = nlay - 1; i >= 0; i--)
    {
        T[i] = Theta[i] * pow((((p[i] + p[i + 1]) / 2.0) / 100000.0), 2.0 / 7.0);
    }
}

double calc_epsilon(double tau, double mu)
{
    return 1.0 - exp(-tau / mu);
}

double calc_B(double T, double lambda, double wavelength)
{
    double h;
    double c;
    double k_B;
    wavelength = wavelength * 1e-9;                                                                                            /* convert from nm to m */
    return (2.0 * h * c * c / (lambda * lambda * lambda * lambda * lambda)) / (exp(h * c / (lambda * k_B * T)) - 1.0) / 1.0e9; /* W/(m2 nm sterad) */
}



void eddington_v2(double dtau, double g, double omega0, double mu0,
                  double *t, double *r, double *rdir, double *sdir, double *tdir)
{
    double alpha1 = 0, alpha2 = 0, alpha3 = 0, alpha4 = 0, alpha5 = 0, alpha6 = 0;
    double a11 = 0, a12 = 0, a13 = 0, a23 = 0, a33 = 0;
    double lambda = 0, b = 0, A = 0;
    double denom = 0;

    /* first, avoid omega0=1 because of instability */
    if (omega0 > 0.999999)
        omega0 = 0.999999;

    alpha1 = (1.0 - omega0) + 0.75 * (1.0 - omega0 * g);
    alpha2 = -(1.0 - omega0) + 0.75 * (1.0 - omega0 * g);

    lambda = sqrt(alpha1 * alpha1 - alpha2 * alpha2);

    A = 1.0 / (alpha2 / (alpha1 - lambda) * exp(lambda * dtau) - alpha2 / (alpha1 + lambda) * exp(-lambda * dtau));

    a11 = A * 2.0 * lambda / alpha2;
    a12 = A * (exp(lambda * dtau) - exp(-lambda * dtau));

    b = 0.5 - 0.75 * g * mu0;
    alpha3 = -omega0 * b;
    alpha4 = omega0 * (1 - b);

    denom = (1.0 / mu0 / mu0 - lambda * lambda);
    alpha5 = ((alpha1 - 1.0 / mu0) * alpha3 - alpha2 * alpha4) / denom;
    alpha6 = (alpha2 * alpha3 - (alpha1 + 1.0 / mu0) * alpha4) / denom;

    a33 = exp(-dtau / mu0);

    a13 = alpha5 * (1.0 - (a11) * (a33)) - alpha6 * (a12);
    a23 = -(a12)*alpha5 * (a33) + alpha6 * ((a33) - (a11));

    *t = a11;
    *r = a12;
    *tdir = a33;
    *rdir = a13 / mu0;
    *sdir = a23 / mu0;
}
