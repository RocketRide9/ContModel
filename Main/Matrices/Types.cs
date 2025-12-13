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

namespace Types;
public interface IMatrix
{
    int Size { get; }
    // TODO: нужно для предобуславливания.
    // Надо придумать что-то более разумное
    Span<Real> Di { get; }

    SparkAlgos.Types.Matrix GetComputeMatrix();
    void Mul(ReadOnlySpan<Real> vec, Span<Real> res);
    // не нулевый, потому что так проще
    IEnumerable<Real> FlatNonZero();
}

// public interface IPatchable<T>
// {
    
// }
