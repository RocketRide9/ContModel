using Real = double;

namespace CoordSystem.Dim2;

public class RZ : ICoordSystem
{
    public static Real Jacobian(Real r, Real z)
    {
        return r;
    }
}
