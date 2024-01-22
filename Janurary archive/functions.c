// ##########################################################
// ############ functions ###################################
// ##########################################################

void sort_vertically(double Theta[], int nlay)
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

double get_Theta_from_T(double Theta[], double T[])
{
    for (int i = nlay - 1; i >= 0; i--)
    {
        Theta[i] = T[i] * pow((100000.0 / ((p[i] + p[i + 1]) / 2)), 2.0 / 7.0);
        return Theta[i];
    }
}

double get_T_from_Theta(double Theta[], double T[])
{
    for (int i = nlay - 1; i >= 0; i--)
    {
        T[i] = Theta[i] * pow((((p[i] + p[i + 1]) / 2.0) / 100000.0), 2.0 / 7.0);
        return T[i];
    }
}

inline double calc_epsilon(double tau, double mu)
{
    return 1.0 - exp(-tau / mu);
}

double calc_B(double T, double lambda)
{

    wavelength *= 1e-9; /* convert from nm to m */

    return (2.0 * h * c * c / (lambda * lambda * lambda * lambda * lambda)) / (exp(h * c / (lambda * k_B * T)) - 1.0) / 1.0e9; /* W/(m2 nm sterad) */
}
