#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using System.Globalization;
using System.Diagnostics;

using SparkCL;
using FemSlaeBuilder;
using SparkAlgos.SlaeSolver;
using SlaeSolver;
using CoordSystem.Dim2;
using TelmaCore;

class Program
{
    const string SRC_DIR = "../../../";
    static void Main(string[] args)
    {
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
        var unixMs = new DateTimeOffset(DateTime.Now).ToUnixTimeMilliseconds();
        Directory.CreateDirectory(SRC_DIR + "measurements");
        var measures = new StreamWriter(SRC_DIR + "measurements/" + unixMs + ".txt");
        Trace.Listeners.Add(new TextWriterTraceListener(measures));
        Trace.AutoFlush = true;
        
        var sw = Stopwatch.StartNew();
        Core.Init();
        Trace.WriteLine($"SparkCL Init: {sw.ElapsedMilliseconds}ms");

        ReverseSigma();
        return;
        ElectroMany();
        Spline();

        Iterate();

        var e = ElectroOnce(7);
        Console.WriteLine(
            string.Join(
                "\n",
                e.Select(
                    (val, idx) => $"{idx}: {val}."
                )
            )
        );
    }
    
    static void Spline()
    {
        var task = new TaskRect4x5XY1();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");

        prob.Build<DiagSlaeBuilder<XY>>();

        var vals = SolveOCL<BicgStabOCL>(prob);

        var splineProb = new ProblemSpline(SRC_DIR + "InputSpline", vals, prob.Mesh);
        splineProb.Build<SplineSlaeBuilder.DiagSlaeBuilderHermit<XY>>();
        
        SplineSolveOCL<BicgStabOCL>(splineProb);
    }
    
    static void SplineSolveOCL<T>(ProblemSpline prob)
    where T: SparkAlgos.SlaeSolver.ISlaeSolver
    {
        var (ans, iters, rr) = prob.SolveOCL<T>();
        
        Console.WriteLine($"(iters {iters}) (discrep: {rr})");
    }
    
    static void ReverseI()
    {
        // TODO: дописать
        // u -> e
        Real[] weights = [1, 1, 1, 1];
        // измеренные напряжённости, с какой-то погрешностью
        PairF64[] e_measured = [
            new(3.710883057100038e-07, -4.652739107054252e-08),
            new(1.046909531845134e-07, -6.388825590065409e-09),
            new(4.878075851331915e-08, -1.811845815180536e-09),
            new(2.858623606945764e-08, -7.38675306794607e-10)
        ];
        // начальное значение проводимости
        Real u_guide = 6.5; // 6
        // истинные напряжённости
        PairF64[] e_true = [
            new(3.710883057100038e-07, -4.652739107054252e-08),
            new(1.046909531845134e-07, -6.388825590065409e-09),
            new(4.878075851331915e-08, -1.811845815180536e-09),
            new(2.858623606945764e-08, -7.38675306794607e-10)
        ];
        // e_measured = e_true;
        Real alpha = 0;
        
        Real j0 = 0;
        var e0 = ElectroOnce(u_guide);
        for (int i = 0; i < 4; i++)
        {
            var w = (e_measured[i] - e0[i]) * weights[i] / e_measured[i].Norm;
            j0 += w*w;
        }
        
        Console.WriteLine($"j0: {j0}");
        
        Real beta = 1;
        Real u0 = u_guide;

        int iter = 0;
        while (true)
        {
            // дифференциал
            var diffU = new PairF64[4];
            var e1 = ElectroOnce(1.05*u0);
            for (int i = 0; i < 4; i++)
            {  
                diffU[i] = (e0[i]-e1[i])/(0.05*u0);
            }
            
            // A
            Real a = alpha;
            for (int i = 0; i < 4; i++)
            {
                var w = weights[i] / e_measured[i].Norm;
                a += w*w * diffU[i]*diffU[i];
            }
            
            // f
            Real f = -alpha * (u0 - u_guide);
            for (int i = 0; i < 4; i++)
            {
                var w = weights[i] / e_measured[i].Norm;
                f -= w*w * (e_measured[i] - e0[i]) * diffU[i];
            }
            
            // новое решение
            var du = f/a;
            u0 += beta * du;
            
            // новое значение функционала
            Real j1 = 0;
            e0 = ElectroOnce(u0);
            for (int i = 0; i < 4; i++)
            {
                var w = (e_measured[i] - e0[i]) * weights[i] / e_measured[i].Norm;
                j1 += w*w;
            }

            iter++;
            if (j0 < j1) {
                beta /= 2.0;
                if (beta < 1.0/4.0) {
                    break;
                }
            }
            
            Console.WriteLine($"sigma: {u0}");
            Console.WriteLine($"j: {j1}");
            
            j0 = j1;
        }
        
        Console.WriteLine($"Iterations: {iter}");
    }
    
    static void ReverseSigma()
    {
        // u -> e
        Real[] weights = [1, 1, 1, 1];
        // измеренные напряжённости, с какой-то погрешностью
        PairF64[] e_measured = [
            new(3.710883057100038e-07, -4.652739107054252e-08),
            new(1.046909531845134e-07, -6.388825590065409e-09),
            new(4.878075851331915e-08, -1.811845815180536e-09),
            new(2.858623606945764e-08, -7.38675306794607e-10)
        ];
        // начальное значение проводимости
        Real u_guide = 6.5; // 6
        // истинные напряжённости
        PairF64[] e_true = [
            new(3.710883057100038e-07, -4.652739107054252e-08),
            new(1.046909531845134e-07, -6.388825590065409e-09),
            new(4.878075851331915e-08, -1.811845815180536e-09),
            new(2.858623606945764e-08, -7.38675306794607e-10)
        ];
        // e_measured = e_true;
        Real alpha = 0;
        
        Real j0 = 0;
        var e0 = ElectroOnce(u_guide);
        for (int i = 0; i < 4; i++)
        {
            var w = (e_measured[i] - e0[i]) * weights[i] / e_measured[i].Norm;
            j0 += w*w;
        }
        
        Console.WriteLine($"j0: {j0}");
        
        Real beta = 1;
        Real u0 = u_guide;

        int iter = 0;
        while (true)
        {
            // дифференциал
            var diffU = new PairF64[4];
            var e1 = ElectroOnce(1.05*u0);
            for (int i = 0; i < 4; i++)
            {  
                diffU[i] = (e0[i]-e1[i])/(0.05*u0);
            }
            
            // A
            Real a = alpha;
            for (int i = 0; i < 4; i++)
            {
                var w = weights[i] / e_measured[i].Norm;
                a += w*w * diffU[i]*diffU[i];
            }
            
            // f
            Real f = -alpha * (u0 - u_guide);
            for (int i = 0; i < 4; i++)
            {
                var w = weights[i] / e_measured[i].Norm;
                f -= w*w * (e_measured[i] - e0[i]) * diffU[i];
            }
            
            // новое решение
            var du = f/a;
            u0 += beta * du;
            
            // новое значение функционала
            Real j1 = 0;
            e0 = ElectroOnce(u0);
            for (int i = 0; i < 4; i++)
            {
                var w = (e_measured[i] - e0[i]) * weights[i] / e_measured[i].Norm;
                j1 += w*w;
            }

            iter++;
            if (j0 < j1) {
                beta /= 2.0;
                if (beta < 1.0/4.0) {
                    break;
                }
            }
            
            Console.WriteLine($"sigma: {u0}");
            Console.WriteLine($"j: {j1}");
            
            j0 = j1;
        }
        
        Console.WriteLine($"Iterations: {iter}");
    }
    
    static PairF64[] ElectroOnce(Real sigma)
    {
        var task = new TaskElectro();
        task.Sigma = sigma;
        var prob = new ProblemLine(task, SRC_DIR + "InputElectro");
        
        // prob.MeshRefine(new()
        // {
        //     XSplitCount   = [  8,   8,   8,   8, 120],
        //     XStretchRatio = [1.0, 1.0, 1.0, 1.0, 1.0],
        //     YSplitCount   = [120,   8],
        //     YStretchRatio = [1.0, 1.0],
        // });
        
        prob.buildType = GlobalMatrixImplType.Host;
        prob.Build<DiagSlaeBuilder<RZ>>();
        // ElectroSolveHost<CgmHost>(prob);
        return ElectroSolveOCL<BicgStabOCL>(prob);
    }
    
    static void ElectroMany()
    {
        const int REFINE_COUNT = 6;
        
        var task = new TaskElectro();
        var prob = new ProblemLine(task, SRC_DIR + "InputElectro");

        for (int i = 0; i < REFINE_COUNT; i++)
        {
            prob.buildType = GlobalMatrixImplType.Host;
            prob.Build<DiagSlaeBuilder<RZ>>();
            // ElectroSolveHost<CgmHost>(prob);
            var sens = ElectroSolveOCL<BicgStabOCL>(prob);

            for (int k = 0; k < sens.Length; k++)
            {
                Console.WriteLine($"Sensor {k}: {sens[k]}.");
            }

            prob.MeshDouble();
        }
    }
    
    static PairF64[] ElectroSolveOCL<T>(ProblemLine prob)
    where T: SparkAlgos.SlaeSolver.ISlaeSolver
    {
        var (ans, iters, rr) = prob.SolveOCL<T>();
        
        Console.WriteLine($"(iters {iters}) (discrep: {rr})");

        // sensors
        Real[] sensR = [250, 500, 750, 1000];
        
        var sensVals = new PairF64[4];
        var res = new Real[4];
        var resNew = new Real[4];
        for (int i = 0; i < sensR.Length; i++)
        {
            sensVals[i] = -prob.GradAt(ans, sensR[i], 0);
            res[i] = prob.ResultAt(ans, sensR[i], 0);
            resNew[i] = prob.ResultAtNew(ans, sensR[i], 0);
        }
        return sensVals;
    }
    
    static void ElectroSolveHost<T>(ProblemLine prob)
    where T: SlaeSolver.ISlaeSolver
    {
        var (ans, iters, rr) = prob.SolveHost<T>();
        
        Console.WriteLine($"(iters {iters}) (discrep: {rr})");

        // sensors
        Real[] sensR = [250, 500, 750, 1000];

        for (int i = 0; i < sensR.Length; i++)
        {
            // var res = ans[sensIdx[i]];
            var res = prob.ResultAt(ans, sensR[i], 0);
            Console.WriteLine($"Sensor {i+1}: {res}");
        }

        Console.WriteLine();
    }
    
    static void Iterate()
    {
        const int REFINE_COUNT = 3;

        var task = new TaskRect4x5RZ2();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");

        // prob.MeshRefine(new()
        // {
        //     XSplitCount = [32],
        //     YSplitCount = [32],
        //     XStretchRatio = [1],
        //     YStretchRatio = [1],
        // });

        for (int i = 0; i < REFINE_COUNT; i++)
        {
            prob.buildType = GlobalMatrixImplType.Host;
            prob.Build<DiagSlaeBuilder<RZ>>();
            SolveOCL<BicgStabOCL>(prob);

            prob.MeshDouble();
        }
    }
    
    static Real[] SolveOCL<T>(ProblemLine prob)
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

        return ans;
    }
}
