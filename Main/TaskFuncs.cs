/*
MathShards
Copyright (C) 2025 Afonin Anton

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

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
