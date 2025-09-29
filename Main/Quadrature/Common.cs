namespace Quadrature;

public record Node<T>
(T Point, double Weight)
where T : notnull;

public record Quadrature<T> (Node<T>[] Nodes)
where T : notnull;
