using Real = double;

using System.Text.Json.Serialization;

[JsonSerializable(typeof(SolverParams))]
[JsonSerializable(typeof(Real))]
[JsonSerializable(typeof(int))]        
internal partial class SolverParamsSourceGenerationContext : JsonSerializerContext
{
}

public struct SolverParams
{
    public Real eps { get; set; }
    public int maxIter { get; set; }
}

public struct Subdomain
{
    // номера подобласти начинаются с 0
    public int Num;
    public int X1;
    public int X2;
    public int Y1;
    public int Y2;
}

[JsonSerializable(typeof(RefineParams))]
[JsonSerializable(typeof(int[]))]
[JsonSerializable(typeof(Real[]))]        
internal partial class RefineParamsSourceGenerationContext : JsonSerializerContext
{
}

public struct RefineParams
{
    public int[] XSplitCount { get; set; }
    public Real[] XStretchRatio { get; set; }

    public int[] YSplitCount { get; set; }
    public Real[] YStretchRatio { get; set; }
}

[JsonSerializable(typeof(int))]
[JsonSerializable(typeof(BoundaryCondition[]))]        
[JsonSerializable(typeof(BoundaryConditionsFile))]        
internal partial class BoundaryConditionsFileSourceGenerationContext
    : JsonSerializerContext
{}

public struct BoundaryConditionsFile
{
    public BoundaryCondition[] BoundaryConditions { get; set; }
}

public struct BoundaryCondition
{
    // в файлах нумерация с 1
    // в программе - с 0
    public int Num { get; set; }
    // тип краевого условия (первый, второй ...) нумеруется с 1
    public int Type { get; set; }
    public int X1 { get; set; }
    public int X2 { get; set; }
    public int Y1 { get; set; }
    public int Y2 { get; set; }
}

#if false
public struct MsrMatrix
{
    public Real[] Elems;
    public Real[] Di;
    public int[] Ia;
    public int[] Ja;

    public SparkAlgos.Types.MsrMatrixRef AsRef() => new() {
        Elems = Elems,
        Di = Di,
        Ia = Ia,
        Ja = Ja,
    };

    static void ArraySerialize(Real[] arr, string fileName)
    {
        var stream = new StreamWriter(fileName);
        stream.Write(
            string.Join(
                "\n",
                arr.Select(
                    e =>
                    {
                        var bytes = BitConverter.GetBytes(e);
                        var i = BitConverter.ToInt64(bytes, 0);
                        return "0x" + i.ToString("X");
                    }
                )
            )
        );
        stream.Close();
    }

    static void ArraySerialize(int[] arr, string fileName)
    {
        var stream = new StreamWriter(fileName);
        stream.Write(
            string.Join(
                "\n",
                arr
            )
        );
        stream.Close();
    }

    public void Serialize()
    {
        ArraySerialize(Elems, "mat.txt");
        ArraySerialize(Di, "di.txt");
        ArraySerialize(Ia, "ia.txt");
        ArraySerialize(Ja, "ja.txt");
    }
}
#endif

public struct MeshAxes
{
    public Real[] xAxis;
    public Real[] yAxis;
}

interface IElement
{
    int NumberOfDofs { get; }
    Real[,] LocalGMatrix (Real hx, Real hy, Real gamma);
    Real[,] LocalMMatrix (Real hx, Real hy, Real gamma);
}

public struct SplineParams
{
    public Real Alpha { get; set; }
    public Real Beta { get; set; }
    public Real W { get; set; }
}
