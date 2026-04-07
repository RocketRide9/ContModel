using NUnit.Framework.Internal;
using Silk.NET.Core.Native;
using SparkAlgos;
using SparkCompute;
using SparkCU;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
        var prg = ComputeProgram.FromFilename("Test");
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
    public static void SimpleFill()
    {
        var prg = ComputeProgram.FromFilename("Test");
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
    public static void BlasTests()
    {
        var host = new Real[100];
        for (int i = 0; i < host.Length; i++)
        {
            host[i] = i % 10;
        }
        var dev = new ComputeBuffer<Real>(host, BufferFlags.OnHostAndDevice);
        dev.ToDevice();

        var blas = Blas.GetInstance();
        blas.Scale(0.5, dev);

        dev.ToHost();
        for (int i = 0; i < host.Length; i++)
        {
            var acc = dev.MapHost(MapFlags.Read);
            Assert.That(acc[i]*2, Is.EqualTo(i % 10).Within(1e-5));
        }
    }
}
