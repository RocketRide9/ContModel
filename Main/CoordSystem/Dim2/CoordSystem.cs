using Real = double;

namespace CoordSystem.Dim2;

public interface ICoordSystem
{
    static abstract Real Jacobian(double a, double b);
}
