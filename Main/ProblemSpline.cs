#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using System.Text.Json;
using System.Text.Unicode;
using System.Diagnostics;

using Types;
using TelmaCore;
using static FiniteElements.Rectangle.Lagrange.BiLinear;
using SplineSlaeBuilder;

class ProblemSpline
{
    RectMesh _inMesh;
    Real[] _values;
    RectMesh _mesh;
    MeshAxes _meshAxes;
    SplineParams _splineParams;
    SolverParams _solverParams;
    
    public IMatrix matrix;
    public Real[] b;
    
    public ProblemSpline(
        string taskFolder,
        Real[] values, RectMesh inMesh
    ) {
        var json = File.ReadAllText(Path.Combine(taskFolder, "SplineParams.json"));
        _splineParams = JsonSerializer.Deserialize<SplineParams>(json)!;

        json = File.ReadAllText(Path.Combine(taskFolder, "SolverParams.json"));
        _solverParams = JsonSerializer.Deserialize<SolverParams>(json)!;
        
        _meshAxes = ReadMesh(taskFolder);

        _mesh = new RectMesh(
            _meshAxes.xAxis,
            _meshAxes.yAxis
        );
        
        _values = values;
        _inMesh = inMesh;
    }
    
    public void Build<T>()
    where T: ISplineSlaeBuilder
    {
        var sw = Stopwatch.StartNew();
        var slaeBuilder = T.Construct(_inMesh, _values, _mesh);
        (matrix, b) = slaeBuilder.Build();
        Trace.WriteLine($"ProblemLine.Build total: {sw.ElapsedMilliseconds}");
    }
    
    public (Real[] ans, int iters, Real rr) SolveOCL<T> ()
    where T: SparkAlgos.SlaeSolver.ISlaeSolver
    {
        var sw = Stopwatch.StartNew();
        var x0 = Enumerable.Repeat((Real)0, b.Length).ToArray();
        var solver = T.Construct(
            _solverParams.maxIter,
            _solverParams.eps
        );
        solver.AllocateTemps(x0.Length);
        Trace.WriteLine($"Solver prepare: {sw.ElapsedMilliseconds}ms");
        
        sw.Restart();
        
        var cm = matrix.GetComputeMatrix();
        Trace.WriteLine($"Matrix Host->Device: {sw.ElapsedMilliseconds}ms");
        
        sw.Restart();
        var (rr, iter) = solver.Solve(cm, b, x0);
        Trace.WriteLine($"Solver {typeof(T).FullName}: {sw.ElapsedMilliseconds}ms");

        return (x0, iter, rr);
    }
    
    static MeshAxes ReadMesh(string taskFolder)
    {
        var file = new StreamReader(Path.Combine(taskFolder, "Mesh.txt"));

        MeshAxes res;
        res.xAxis = file.ReadLine()!.Split().Select(Real.Parse).ToArray();
        res.yAxis = file.ReadLine()!.Split().Select(Real.Parse).ToArray();

        return res;
    }
}
