using NUnit.Framework.Internal;
using SparkAlgos;
using SparkCompute;
using System.Diagnostics;
using System.Globalization;

namespace Main.Tests.ComputeTests;

[TestFixture]
class ComputeTests
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

    [Test]
    public static void Ids()
    {
        var prg = ComputeProgram.FromFilename("Test.cl");
        var test1d = prg.GetKernel(
            "write_ids_1d",
            globalWork: new NDRange(128).PadTo(16),
            localWork: new NDRange(16)
        );
        int[] zeros = Enumerable.Repeat(0, 200).ToArray();
        var buf_ids_1d = new ComputeBuffer<int>(zeros, BufferFlags.OnHostAndDevice);

        test1d.SetArg(0, buf_ids_1d);
        test1d.SetArg(1, 128);

        test1d.Execute();

        buf_ids_1d.ToHost();
        var accessor = buf_ids_1d.MapHost(MapFlags.Read);
        for (int i = 0; i < 128; i++)
        {
            Assert.That(accessor[i], Is.EqualTo(i));
        }
        for (int i = 128; i < buf_ids_1d.Length; i++)
        {
            Assert.That(accessor[i], Is.EqualTo(0));
        }
    }

    [Test]
    public static void Fill()
    {
        var prg = ComputeProgram.FromFilename("Test.cl");
        var test1d = prg.GetKernel(
            "fill",
            globalWork: new NDRange(128).PadTo(16),
            localWork: new NDRange(16)
        );
        var buf = new ComputeBuffer<Real>((int)test1d.GlobalWork[0], BufferFlags.OnHostAndDevice);

        test1d.SetArg(0, buf);
        test1d.SetArg(1, (Real)3);
        test1d.SetArg(2, buf.Length);

        test1d.Execute();

        buf.ToHost();
        var accessor = buf.MapHost(MapFlags.Read);
        for (int i = 0; i < buf.Length; i++)
        {
            Assert.That(accessor[i], Is.EqualTo((Real)3));
        }
    }

    [Test]
    public static void BlasScale()
    {
        var host = new Real[10];
        for (int i = 0; i < host.Length; i++)
        {
            host[i] = i % 10;
        }
        var dev = new ComputeBuffer<Real>(host, BufferFlags.OnHostAndDevice);

        var blas = Blas.GetInstance();
        blas.Scale(0.5, dev);

        dev.ToHost();
        for (int i = 0; i < host.Length; i++)
        {
            var acc = dev.MapHost(MapFlags.Read);
            Assert.That(acc[i]*2, Is.EqualTo(i % 10).Within(1e-5));
        }
    }


    [Test]
    public static void BlasDot()
    {
        var host0 = new Real[100];
        var host1 = new Real[100];
        Real ans = 0;
        for (int i = 0; i < host0.Length; i++)
        {
            host0[i] = i % 10;
            host1[i] = (i+5) % 10;
            ans += host0[i] * host1[i];
        }          
        var dev0 = new ComputeBuffer<Real>(host0, BufferFlags.OnHostAndDevice);
        var dev1 = new ComputeBuffer<Real>(host1, BufferFlags.OnHostAndDevice);

        var blas = Blas.GetInstance();
        blas.Scratch1 = new ComputeBuffer<Real>(1, BufferFlags.OnDevice);
        blas.Scratch64 = new ComputeBuffer<Real>(64, BufferFlags.OnDevice);
        var dot = blas.Dot(dev1, dev0);

        Assert.That(dot, Is.EqualTo(ans).Within(1e-5));
    }
}
