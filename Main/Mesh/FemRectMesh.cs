using Real = double;

// TODO: можно ли обойтись без наследования?
public class FemRectMesh : RectMesh
{
    Subdomain[] _subDomains;
    public Subdomain[] SubDomains { get => _subDomains; }

    BoundaryCondition[] _boundaryConditions;
    public BoundaryCondition[] BoundaryConditions { get => _boundaryConditions; }

    public FemRectMesh(
        Real[] xAxis, Real[] yAxis,
        Subdomain[] subDomains,
        BoundaryCondition[] boundaryConditions
    ) : base(xAxis, yAxis) {
        _subDomains = subDomains;
        _boundaryConditions = boundaryConditions;
    }    

    public int? GetSubdomNumAtElCoords (int x1, int y1)
    {
        foreach (var a in SubDomains)
        {
            if (x1 >= IXw[a.X1] && x1 < IXw[a.X2] &&
                y1 >= IYw[a.Y1] && y1 < IYw[a.Y2]
            ) {
                return a.Num;
            }
        }

        return null;
    }

    public int? GetSubdomNumAtPoint (Real x1, Real y1)
    {
        foreach (var a in SubDomains)
        {
            if (x1 >= Xw[a.X1] && x1 <= Xw[a.X2] &&
                y1 >= Yw[a.Y1] && y1 <= Yw[a.Y2]
            ) {
                return a.Num;
            }
        }

        return null;
    }
}
