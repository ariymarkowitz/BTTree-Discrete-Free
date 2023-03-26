AttachSpec("spec");

/**
 * Take a sequence of axes and determine whether they satisfy the conditions of the
 * Ping-Pong Lemma.
 */
SatisfiesPPL := function(tree, gens)
  for i -> g in gens do
    l := TranslationLength(g);
    if l eq 0 then return false; end if;
    
    dist := 0;
    init := false;

    for j -> h in gens do
      if i eq j then continue; end if;

      // Find the endpoints of the projection of the axis of h onto the axis of g.
      intersects, pathEnds := CrossPath(tree, g, h);
      if intersects then
        points := pathEnds;
      else
        points := <pathEnds[1]>;
      end if;

      // Update `min` or `max` depending on where `point` lies on the axis.
      for point in points do
        if init eq false then
          min := point;
          max := point;
          init := true;
        else
          ldist := Distance(point, min);
          rdist := Distance(point, max);
          if rdist ge dist and ldist ge rdist then
            // `point` is on the far side of `min`.
            min := point;
            dist := rdist;
          elif ldist ge dist and rdist ge ldist then
            // `point` is on the far side of `max`.
            max := point;
            dist := ldist;
          end if;

          if dist ge l then
            print i, j;
            return false;
          end if;
        end if;
      end for;
    end for;
  end for;

  return true;
end function;

RandomGens := function(n, K, B)
  p := UniformizingElement(K);
  Z := Integers(); gens := [];
  M := MatrixAlgebra(K, 2);
  i := 1;
  while i le n do
    repeat a := Random(Z,B) * p^Random(Z,B); until a ne 0;
    b := Random(Z,B) * p^Random(Z,B);
    c := Random(Z,B) * p^Random(Z,B);
    d := (1 + b * c) / a;
    A := M![a, b, c, d];
    if not IsElliptic(A) then
      gens[i] := A;
    end if;
    i +:= 1;
  end while;
  return gens;
end function;

p := 5;
Qp := pAdicField(p,100);
T := BruhatTitsTree(Qp);

// Verify that the algorithm is correct.
count := 0;
while true do
  gens := RandomGens(5,Qp,5);
  count +:= 1;
  result, reducedGens := IsFreeDiscrete(T, gens);
  print count, result;
  if result and not SatisfiesPPL(T, reducedGens) then
    break;
  end if;
end while;

// Time the program
times := 100;
for j in [1..10] do
  gensList := [RandomGens(2*j-1, Qp, 3) : i in [1..times]];
  t := Cputime();
  for gens in gensList do
    b, result := IsDiscreteFree(T, gens);
  end for;
  print Cputime(t)/times;
end for;

// Bug!
p := 5;
Qp := pAdicField(p,109);
T := BruhatTitsTree(Qp);

g := [ Matrix(2, [Qp!(-2/3125), Qp!(-62500), Qp!(62500), Qp!(-3944304526105059027058642826413931148366032175545115023845291137696875)]),
Matrix(2, [Qp!( -9765625), Qp!(0), Qp!(3125000), Qp!(-1/9765625)]),
Matrix(2, [Qp!(-97656250), Qp!(8/3125), Qp!(-7/78125), Qp!3944304526105059027058642826413931148366032175545115023851394531250028/11920928955078125]),
Matrix(2, [Qp!(2/625), Qp!(781250), Qp!(-6250), Qp!(3944304526105059027058642826413931148366032175545115023849868774414375)]),
Matrix(2, [Qp!( -97656250), Qp!(225), Qp!(3/390625), Qp!(-7826/762939453125)]) ];

G := [];
X := g;
b := 1;
while b eq 1 do
  Append(~G, X);
  b, X := ReduceGenerators(T, G[#G]);
end while;
