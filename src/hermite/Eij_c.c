int E(int i, int j, int t, int r, double alpha, double beta)
{
    // Hermite Coefficient
    double p, mu

    p  = alpha + beta
    mu = alpha * beta / p

    if ((t < 0) || (t > (i + j)))
        {
        return 0.0;
        }
    else if (i == j) && (j == t) && (j == 0)
        {
        // KAB
        return exp(-mu * r * r);
        }
    else if (j == 0)
        {
        // decrement index i
        return (
            (1 / (2 * p)) * E(i - 1, j, t - 1, r, alpha, beta)
            - (mu * r / alpha) * E(i - 1, j, t, r, alpha, beta)
            + (t + 1) * E(i - 1, j, t + 1, r, alpha, beta)
        );
        }
    else
        {
        // decrement index j
        return (
            (1 / (2 * p)) * E(i, j - 1, t - 1, r, alpha, beta)
            + (mu * r / beta) * E(i, j - 1, t, r, alpha, beta)
            + (t + 1) * E(i, j - 1, t + 1, r, alpha, beta)
        );
        }
}