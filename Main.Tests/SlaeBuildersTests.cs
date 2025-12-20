using static MathShards.Quadrature.Gauss;
using MathShards.TelmaCore;
using MathShards.SlaeBuilder.Fem;
using MathShards.CoordSystem.Dim2;

namespace Main.Tests.SlaeBuildersTests;

[TestFixture]
class SlaeBuilders
{
    const string MAIN_SRC_DIR = "../../../../Main/";

    [SetUp]
    public void Setup() {
        Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.InvariantCulture;

    }
    
    [Test]
    public void MsrVsDiag ()
    {
        const int REPEAT_COUNT = 3;
        const int REFINE_COUNT = 1;
        
        var task = new TaskRect4x5XY1();
        var prob = new ProblemLine(task, MAIN_SRC_DIR + "InputRect4x5");
        prob.MeshRefine(new()
        {
            XSplitCount = [1024/64],
            YSplitCount = [1024/64],
            XStretchRatio = [1],
            YStretchRatio = [1],
        });
        
        prob.Build<DiagSlaeBuilder<XY>>();
        var nzd = prob.matrix.FlatNonZero();
        
        prob.Build<MsrSlaeBuilder>();
        var nzm = prob.matrix.FlatNonZero();
        
        foreach (var (first, second) in nzd.Zip(nzm))
        {
            Assert.That(first, Is.EqualTo(second).Within(1e-12));
        }
    }

}
