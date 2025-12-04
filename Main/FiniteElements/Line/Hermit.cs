using Real = double;

using TelmaCore;
using static Quadrature.Gauss;

namespace FiniteElements.Line.Hermit;

public static class Cubic
{
    public static readonly Func<Real, Real>[] BasisTemplate =
    {
        a => 1 - 3*a*a + 2*a*a*a,
        a => a - 2*a*a + a*a*a,
        a => 3*a*a - 2*a*a*a,
        a => -a*a + a*a*a
    };
    
    public static Real BasisConverted(int i, Real p0, Real p1, Real p)
    {
        Real h = p1 - p0;
        Real p01 = (p - p0)/h;
        // см. кирпич с.151
        Real[] coeffs = [1, h, 1, h];

        return coeffs[i] * BasisTemplate[i](p01);
    }
    
    public static readonly Func<Real, Real>[] BasisGradTemplate =
    {
        a => -6*a + 6*a*a,
        a => 1 - 4*a + 3*a*a,
        a => 6*a - 6*a*a,
        a => -2*a + 3*a*a
    };
    
    public static Real BasisGradConverted(int i, Real p0, Real p1, Real p)
    {
        Real h = p1 - p0;
        Real p01 = (p - p0)/h;
        // см. кирпич с.151
        Real[] coeffs = [1, h, 1, h];

        return coeffs[i] * BasisGradTemplate[i](p01);
    }
    
    public static readonly Func<Real, Real>[] BasisGradGradTemplate =
    {
        a => -6 + 12*a,
        a => -4 + 6*a,
        a => 6 - 12*a,
        a => -2 + 6*a
    };
    
    public static Real BasisGradGradConverted(int i, Real p0, Real p1, Real p)
    {
        Real h = p1 - p0;
        Real p01 = (p - p0)/h;
        // см. кирпич с.151
        Real[] coeffs = [1, h, 1, h];

        return coeffs[i] * BasisGradGradTemplate[i](p01);
    }
}
