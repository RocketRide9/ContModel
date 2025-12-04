using Real = double;

public class LineMesh
{
    RefineParams1D _refineParams;
    
    public Real[] Xw;
    public Real[] X { get; private set; }
    public int[] IXw { get; private set; }
    
    public int nodesCount;
    public int feCount;

    public LineMesh(
        Real[] axis
    ) {
        Xw = axis;

        X = (Real[])Xw.Clone();
        IXw = Enumerable.Range(0, X.Length).ToArray();

        nodesCount = X.Length;
        feCount = (X.Length - 1);
        
        _refineParams = new()
        {
            XSplitCount = Enumerable.Repeat(1, X.Length - 1).ToArray(),
            XStretchRatio = Enumerable.Repeat((Real)1, X.Length - 1).ToArray(),
        };
    }
    
    public void RefineDiv2()
    {
        for (int i = 0; i < _refineParams.XSplitCount.Length; i++)
        {
            _refineParams.XSplitCount[i] *= 2;
            _refineParams.XStretchRatio[i] = (Real)Math.Sqrt(_refineParams.XStretchRatio[i]);
        }

        Refine(_refineParams);
    }
    
    static Real FirstStepSize(Real stretch, int seg_count, Real gap)
    {
        Real sum;
        if (stretch != 1d)
        {
            sum = (Real)(1 - Math.Pow(stretch, seg_count)) / (1 - stretch);
        } else {
            sum = seg_count;
        }

        return gap / sum;
    }
    
    public void Refine(RefineParams1D refineParams)
    {
        _refineParams = refineParams;

        { // ось X
            var xLength = _refineParams.XSplitCount.Sum() + 1;
            IXw = new int[Xw.Length];
            X = new Real[xLength];

            IXw[0] = 0;
            X[0] = Xw[0];
            int xCount = 1;
            for (int i = 1; i < Xw.Length; i++)
            {
                Real gap = Xw[i] - Xw[i - 1];

                int seg_count = _refineParams.XSplitCount[i - 1];
                Real stretch = _refineParams.XStretchRatio[i - 1];

                var step = FirstStepSize(stretch, seg_count, gap);
                var step_n = step;
                var stretch_n = stretch;
                int idx = xCount - 1;
                for (int j = 0; j < seg_count - 1; j++)
                {
                    X[xCount] = X[idx] + step_n;
                    xCount++;
                    stretch_n *= stretch;
                    if (stretch != 1d)
                    {
                        step_n = step * (stretch_n - 1) / (stretch - 1);
                    } else {
                        step_n = step * (j + 2);
                    }
                }
                IXw[i] = xCount;
                X[xCount] = Xw[i];
                xCount++;
            }
        }

        nodesCount = X.Length;
        feCount = (X.Length - 1);
    }
}
