using Real = double;

namespace CoordSystem.Dim2;

public class XY : ICoordSystem
{
    public static Real Jacobian(Real x, Real y)
    {
        return 1;
    }
}
