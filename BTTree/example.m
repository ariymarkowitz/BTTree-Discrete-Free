Attach("BTTree.m");

Qp := pAdicField(2, 20);
M := MatrixAlgebra(Qp, 2);

// Create the Bruhat-Tits tree over the 2-adic field.
tree := BruhatTitsTree(Qp);
print Origin(tree);

// Create a vertex of the tree corresponding to the matrix [1 2; 3 22].
v := BTTVertex(tree, M ! [1, 2, 3, 22]);
print v;
// Get the distance between two vertices.
print Distance(v, Origin(tree));
// Get the standard form matrix corresponding to `v`.
print Matrix(v);
// Get the vertices neighbouring `v`.
n := Neighbors(v);
print n;
// Get the quotient v/w as a representative vector over the residue field.
print TypeOfNeighbor(v, n[1]);

A := M ! [2, 0, 1, 3];
// Apply an isometry to a vertex.
print v*A;
print TranslationLength(A);
print Distance(v, v*A^2);
print Path(v, v*A^2);
print IsHyperbolic(A);
print TranslationAxisBoundary(A);
w := ProjectionOntoMinTranslationSet(tree, A, Vector(Qp, [1, 9]));
print w;
print IsInMinTranslationSet(A, w);

B := M ! [5, 3, 1, 3];
print TranslationLength(B);
print IsInMinTranslationSetBoundary(B, Vector(Qp, [1, 1]));
print FixSet(tree, B);

print FixSet(tree, M ! [1, 0, 0, 1]);
print FixSet(tree, M ! [0, 1, 1, 0]);
print FixSet(tree, M ! [0, 1, -1, -2]);
print FixSet(tree, M ! [2, 1, 1/8, 5]);
print FixSet(tree, M ! [2, 1, 1/16, 5]);
print FixSet(tree, M ! [-3227, -4789, 4661, 7227]);

print Midpoint(tree, Vector(Qp, [1, 1]), Vector(Qp, [1, 0]), Vector(Qp, [1, 12]));

A := IsometryBetweenAxes(tree, Vector(Qp, [1, -1]), Vector(Qp, [1, -1/2]));
B := IsometryBetweenAxes(tree, Vector(Qp, [1, 2]), Vector(Qp, [1, 0]));
b1, length := Intersects(A, B);
b2, path := CrossPath(tree, A, B);
print b2, path, length;
print b1 eq b2 and #path eq length + 1;