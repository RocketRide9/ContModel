using Real = double;

using System.Globalization;
using System.Diagnostics;

using SparkCL;
using CoordSystem.Dim2;
using FemSlaeBuilder;
using SparkAlgos.SlaeSolver;
using Matrices;

namespace Main.Tests;

[TestFixture]
class HalfMultiplyTest
{
    const string MAIN_SRC_DIR = "../../../../Main/";
    
    [SetUp]
    public void Setup()
    {
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
        var unixMs = new DateTimeOffset(DateTime.Now).ToUnixTimeMilliseconds();
        Directory.CreateDirectory(MAIN_SRC_DIR + "measurements");
        var measures = new StreamWriter(MAIN_SRC_DIR + "measurements/" + unixMs + ".txt");
        Trace.Listeners.Add(new TextWriterTraceListener(measures));
        Trace.AutoFlush = true;
        
        Core.Init();
    }
    
    [Test]
    public void HalfMultiplies() {
        var task = new TaskRect4x5x3();
        var prob = new ProblemLine(task, MAIN_SRC_DIR + "InputRect4x5");

        prob.Build<DiagSlaeBuilder<XY>>();

        var (ans, iters, rr) = prob.SolveOCL<BicgStabOCL>();

        var matrix = prob.matrix as Diag9Matrix;

        var vec = new Real[ans.Length];
        var vec2 = new Real[ans.Length];
        // L
        ans.AsSpan().CopyTo(vec);

        matrix!.InvLMul(vec);
        matrix.LMul(vec, vec2);

        var err = vec2
            .Zip(ans)
            .Select(a =>
            {
                return Real.Abs(a.First - a.Second);
            })
            .Sum();
        
        Assert.That(err, Is.EqualTo(0).Within(1e-12));
        
        // U
        ans.AsSpan().CopyTo(vec);
        
        matrix!.InvUMul(vec);
        matrix.UMul(vec, vec2);
        
        err = vec2
            .Zip(ans)
            .Select(a =>
            {
                return Real.Abs(a.First - a.Second);
            })
            .Sum();
        
        Assert.That(err, Is.EqualTo(0).Within(1e-12));
    }
}
