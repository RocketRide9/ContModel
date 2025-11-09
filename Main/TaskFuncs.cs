using Real = double;

public interface TaskFuncs
{
    string Description { get; }

    Real Answer(int subdom, Real x, Real y);
    Real Lambda(int subdom, Real x, Real y);
    Real Gamma(int subdom, Real x, Real y);
    Real F(int subdom, Real x, Real y);
    
    // I
    Real Ug(int bcNum, Real x, Real y);
    // II
    Real Theta(int bcNum, Real x, Real y);
    // III
    Real Beta(int bcNum);
    Real uBeta(int bcNum, Real x, Real y);
}
