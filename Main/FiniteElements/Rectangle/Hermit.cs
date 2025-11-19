using Real = double;

using TelmaCore;
using static Quadrature.Gauss;
using CoordSystem.Dim2;

namespace FiniteElements.Rectangle.Hermit;
using Dim1 = FiniteElements.Line.Hermit.Cubic;

public static class Cubic
{
    // 16 функций
    public static Func<PairF64, Real> Basis(int i)
    {
        int mu = 2*(i/4%2) + i%2;
        int nu = 2*(i/8) + i/2%2;

        return pair
            => Dim1.Basis[mu](pair.X) * Dim1.Basis[nu](pair.Y);
    }
    
    /// xy==0 => diff x
    /// xy==1 => diff y
    public static Func<PairF64, Real> BasisGrad(int i, int xy)
    {
        int mu = 2*(i/4%2) + i%2;
        int nu = 2*(i/8) + i/2%2;

        if (xy == 0)
        {
            return pair
                => Dim1.BasisGrad[mu](pair.X) * Dim1.Basis[nu](pair.Y);
        } else if (xy == 1) {
            return pair
                => Dim1.Basis[mu](pair.X) * Dim1.BasisGrad[nu](pair.Y);
        } else {
            throw new ArgumentException("Invalid xy argument");
        }
    }
    
    public static Real[,] ComputeLocal<Tc>(TaskFuncs funcs, PairF64 p0, PairF64 p1, int subDom)
    where Tc : ICoordSystem
    {
        // side
        var sd = 16;
        var values = new Real[sd, sd];

        for (int i = 0; i < sd; i++)
        {
            for (int j = 0; j < sd; j++)
            {
                var ph = p1 - p0;
                
                var funcMass = (PairF64 point) => {
                    // в координатах шаблонного базиса, [0;1]
                    var p01 = new PairF64 (
                        (point.X - p0.X) / ph.X,
                        (point.Y - p0.Y) / ph.Y
                    );
                    return funcs.Gamma(subDom, point.X, point.Y)
                        * Basis(i)(p01)
                        * Basis(j)(p01)
                        * Tc.Jacobian(point.X, point.Y);
                };
                
                values[i, j] = Integrate2DOrder5(p0, p1, funcMass);
                
                var funcStiffness = (PairF64 point) => {
                    // в координатах шаблонного базиса - [0;1]
                    var p01 = new PairF64 (
                        (point.X - p0.X) / ph.X,
                        (point.Y - p0.Y) / ph.Y
                    );
                    return  funcs.Lambda(subDom, point.X, point.Y)
                        *
                        (
                            BasisGrad(i, 0)(p01)
                            * BasisGrad(j, 0)(p01) / ph.X / ph.X
                        +
                            BasisGrad(i, 1)(p01)
                            * BasisGrad(j, 1)(p01) / ph.Y / ph.Y
                        )
                        * Tc.Jacobian(point.X, point.Y);
                };
                
                values[i, j] += Integrate2DOrder5(p0, p1, funcStiffness);
            }
        }

        return values;
    }
    
    public static Real[] ComputeLocalB<Tc>(TaskFuncs funcs, PairF64 p0, PairF64 p1, int subDom)
    where Tc : ICoordSystem
    {
        // side
        var sd = 16;
        var ph = p1 - p0;
        var res = new Real[sd];
        
        for (int i = 0; i < sd; i++)
        {
            var func = (PairF64 point) =>
            {
                // в координатах шаблонного базиса - [0;1]
                var p01 = new PairF64 (
                    (point.X - p0.X) / ph.X,
                    (point.Y - p0.Y) / ph.Y
                );
                return funcs.F(subDom, point.X, point.Y)
                    * Basis(i)(p01)
                    * Tc.Jacobian(point.X, point.Y);
            };
            res[i] = Integrate2DOrder5(p0, p1, func);
        }
        
        return res;
    }
}
