using TelmaCore;

namespace FiniteElements.Rectangular;

public static class BiLinear
{
    public static readonly Func<PairF64, double>[] LagrangeBasis =
    {
        (PairF64 vert) => (1 - vert.X) * (1 - vert.Y),
        (PairF64 vert) => vert.X       * (1 - vert.Y),
        (PairF64 vert) => (1 - vert.X) * vert.Y,
        (PairF64 vert) => vert.X       * vert.Y
    };
    
    public static readonly Func<PairF64, double>[,] LagrangeBasisGrad =
    {
        {
            (PairF64 vert) => -(1 - vert.Y),
            (PairF64 vert) => -(1 - vert.X),
        },
        {
            (PairF64 vert) => (1 - vert.Y),
            (PairF64 vert) => -vert.X,
        },
        {
            (PairF64 vert) => -vert.Y,
            (PairF64 vert) => (1 - vert.X),
        },
        {
            (PairF64 vert) => vert.Y,
            (PairF64 vert) => vert.X
        }
    };
}
