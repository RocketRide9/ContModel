using Real = double;

using TelmaCore;
using static Quadrature.Gauss;

namespace FiniteElements.Line.Hermit;

public static class Cubic
{
    public static readonly Func<Real, Real>[] Basis =
    {
        a => 1 - 3*a*a + 2*a*a*a,
        a => a - 2*a*a + a*a*a,
        a => 3*a*a - 2*a*a*a,
        a => -a*a + a*a*a
    };
    
    public static readonly Func<Real, Real>[] BasisGrad =
    {
        a => -6*a + 6*a*a,
        a => 1 - 4*a + 3*a*a,
        a => 6*a - 6*a*a,
        a => -2*a + 3*a*a
    };
}
