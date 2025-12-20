using Real = double;

using System.Globalization;
using System.Diagnostics;

using SparkCL;
using MathShards.CoordSystem.Dim2;
using MathShards.SlaeBuilder.Fem;
using SparkAlgos.SlaeSolver;
using MathShards.Matrices;
using Silk.NET.OpenCL;

namespace Main.Tests.HalvesTests;

[TestFixture]
class MsrHalvesTest
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
    
    [TearDown]
    public void Teardown()
    {
        Trace.Listeners.Clear();
        
        Core.Deinit();
    }

    static (MsrMatrix matrix, Real[] b) SomeSlae() {
        Real[] identity = [1, 1, 1, 1, 1];
        Real[] elems = [0.5, 0.5];
        var matrix = new MsrMatrix {
            Elems = [..elems],
            Ia = [0, 1, 2, 2, 2, 2],
            Ja = [1, 0],
            Di = [..identity],
        };
        Real[] b = [..identity];
        b[0] += 0.5;        
        b[1] += 0.5;

        return (matrix, b);
    } 
    
    [Test]
    public static void Host()
    {
        var (matrix, b) = SomeSlae();
        Common.HalfMultiplies(matrix, b);
    }

    [Test]
    public static void OpenCL()
    {
        var (matrix, b) = SomeSlae();
        Common.HalfMultipliesOpenCL(matrix, b);
    }
}
