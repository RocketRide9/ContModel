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

class ProblemSpline1D
{
    LineMesh _inMesh;
    Real[] _values;
    public LineMesh _mesh;
    Real[] _meshAxe;
    SplineParams _splineParams;
    SolverParams _solverParams;
    RefineParams1D _refineParams;
    
    public IMatrix matrix;
    public Real[] b;
    
    public ProblemSpline1D(
        string taskFolder,
        Real[] values, LineMesh inMesh
    ) {
        var json = File.ReadAllText(Path.Combine(taskFolder, "SplineParams.json"));
        _splineParams = JsonSerializer.Deserialize<SplineParams>(json)!;

        json = File.ReadAllText(Path.Combine(taskFolder, "SolverParams.json"));
        _solverParams = JsonSerializer.Deserialize<SolverParams>(json)!;
        
        json = File.ReadAllText(Path.Combine(taskFolder, "RefineParams.json"));
        _refineParams = JsonSerializer.Deserialize<RefineParams1D>(json)!;

        // _meshAxe = ReadMesh(taskFolder);
        _meshAxe = [inMesh.X[0], inMesh.X[inMesh.X.Length - 1]];

        _mesh = new LineMesh(_meshAxe);
        _mesh.Refine(_refineParams);

        _values = values;
        _inMesh = inMesh;
    }
    
    public void Build<T>()
    where T: ISplineSlaeBuilder1D
    {
        var sw = Stopwatch.StartNew();
        var slaeBuilder = T.Construct(_inMesh, _values, _mesh, _splineParams);
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
    
    static Real[] ReadMesh(string taskFolder)
    {
        var file = new StreamReader(Path.Combine(taskFolder, "Mesh.txt"));

        return file.ReadLine()!.Split().Select(Real.Parse).ToArray();
    }
}
