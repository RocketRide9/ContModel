using TelmaCore;

namespace Quadrature;

public static class Gauss
{
    /// p0 - нижний левый угол прямоугольной области
    /// p1 - верхний правый угол
    /// func - на промежутке [-1:1]
    public static double Integrate2DOrder5(
        PairF64 p0, PairF64 p1,
        Func<PairF64, double> func
    ) {
        var quad = Get2DOrder5();

        var res = 0.0;
        foreach (var node in quad.Nodes)
        {
            var point = node.Point;
            var weight = node.Weight;
            res += func(point) * weight;
        }

        // с Якобианом
        return res * (p1.X - p0.X) * (p1.Y - p0.Y) / 4.0;
    }
    
    static Quadrature<PairF64> Get2DOrder5()
    {
        var dim1 = Get1DOrder5();
        return Make2D(dim1);
    }
    
    // https://w.wiki/6kqB
    static Quadrature<double> Get1DOrder5()
    {
        double[] points = {
            -Math.Sqrt(0.6),
             0.0,
             Math.Sqrt(0.6)
        };
        double[] weights = {
            5 / 9.0,
            8 / 9.0,
            5 / 9.0
        };

        var res = new Node<double>[3];

        for (int i = 0; i < 3; i++)
        {
            // есть кусочек от Якобиана и перехода на мастер-элемент
            res[i] = new Node<double>(points[i], weights[i]);
        }

        return new Quadrature<double>(res);
    }
    
    static Quadrature<PairF64> Make2D(Quadrature<double> dim1)
    {
        var len = dim1.Nodes.Length;
        var res = new Node<PairF64>[len * len];

        int i = 0;
        foreach (var node1 in dim1.Nodes)
        {
            foreach (var node2 in dim1.Nodes)
            {
                res[i] = new Node<PairF64> (
                    new(node1.Point, node2.Point),
                    node1.Weight * node2.Weight
                );
                i++;
            }
        }

        return new Quadrature<PairF64>(res);
    }
}
