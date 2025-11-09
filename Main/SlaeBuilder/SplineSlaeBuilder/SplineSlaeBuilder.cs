using Real = double;

using System.Numerics;

using Types;

namespace SplineSlaeBuilder;

public interface ISplineSlaeBuilder
{
    // RectMesh Mesh { get; }

    /// <param name="inMesh">Сетка из решённой задачи</param>
    /// <param name="inValues">Результат задачи</param>
    /// <param name="mesh">Сетка для построения сплайна</param>
    static abstract ISplineSlaeBuilder Construct(
        RectMesh inMesh, Real[] inValues,
        RectMesh mesh
    );
    (IMatrix matrix, Real[] right) Build();
}
