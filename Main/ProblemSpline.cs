/*
MathShards
Copyright (C) 2025 Afonin Anton

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

using System.Text.Json;
using System.Text.Unicode;
using System.Diagnostics;

using MathShards.Matrices.Types;
using MathShards.TelmaCore;
using static MathShards.FiniteElements.Rectangle.Lagrange.BiLinear;
using MathShards.SlaeBuilder.Spline;
using MathShards.Mesh.RectMesh;

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
