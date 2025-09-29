#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using System.Globalization;
using System.Diagnostics;

using SparkCL;
using SlaeBuilder;
using SparkAlgos.SlaeSolver;

class Program
{
    const string SRC_DIR = "../../../";
    static void Main(string[] args)
    {
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
        var unixMs = new DateTimeOffset(DateTime.Now).ToUnixTimeMilliseconds();
        Directory.CreateDirectory(SRC_DIR + "measurements");
        Trace.Listeners.Add(new TextWriterTraceListener(new StreamWriter(SRC_DIR + "measurements/" + unixMs + ".txt")));
        Trace.AutoFlush = true;
        
        var sw = Stopwatch.StartNew();
        Core.Init();
        Trace.WriteLine($"SparkCL Init: {sw.ElapsedMilliseconds}ms");

        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");
        
        prob.MeshRefine(new()
        {
            XSplitCount = [64],
            YSplitCount = [64],
            XStretchRatio = [1],
            YStretchRatio = [1],
        });
        
        prob.buildType = GlobalMatrixImplType.Host;
        prob.Build<DiagSlaeBuilder>();
        SolveOCL<CgmOCL>(prob);
    }
    
    static void SolveOCL<T>(ProblemLine prob)
    where T: SparkAlgos.SlaeSolver.ISlaeSolver
    {
        Trace.WriteLine("SolveOCL");
        Trace.Indent();
#if SPARKCL_COLLECT_TIME
        Core.ResetTime();
#endif
        var sw = Stopwatch.StartNew();
        var (ans, iters, rr) = prob.SolveOCL<T>();
        sw.Stop();
        var err = prob.Lebeg2Err(ans);
#if SPARKCL_COLLECT_TIME
        var (ioTime, kernTime) = Core.MeasureTime();
        ioTime /= (ulong)1e+6;
        kernTime /= (ulong)1e+6;
#endif
        Console.WriteLine($"(err {err}) (iters {iters}) (discrep: {rr})");
        Trace.Unindent();
        Trace.Write($"Solver total: {sw.ElapsedMilliseconds}мс");
#if SPARKCL_COLLECT_TIME
        Trace.Write($": {kernTime}мс + {ioTime}мс");
#endif
        Trace.WriteLine("");
    }
}
