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
    repeat
      repeat a := Random(Z,B) * p^Random(Z,B); until a ne 0;
      b := Random(Z,B) * p^Random(Z,B);
      c := Random(Z,B) * p^Random(Z,B);
      d := (1 + b * c) / a;
      A := M![a, b, c, d];
    until IsHyperbolic(A);
    gens[i] := A;
    i +:= 1;
  end while;
  return gens;
end function;

RandomAut := function(X, nAut)
  if #X eq 1 then
    return X;
  end if;
  for i in [1..nAut] do
    j := Random(1, #X);
    k := random{k : k in [1..#X] | k ne j};
    x := X[j];
    y := X[k];
    X[j] := eval(Random(["x*y", "x*y^-1", "y*x", "y^-1*x"]));
  end for;
  return X;
end function;

// Same as IsDiscreteFree, but returns how many reduction steps were performed.
IsDiscreteFreeCount := function(T, X)
  steps := 1;
  Y := X;
  while true do
    a, b := ReduceGenerators(T, Y);
    if a eq 0 then
      // b is the index of a scalar matrix.
      // Check that Y approximates a scalar matrix.
      A := Y[b];
      if IsApproximatelyIdentity(A) then
        Remove(~Y, b);
      else
        return false, Y[b], steps;
      end if;
    elif a eq 1 then
      Y := b;
    else
      return true, Y, steps;
    end if;
    steps +:= 1;
  end while;
end function;

p := 5;
Qp := pAdicField(p,1000);
T := BruhatTitsTree(Qp);

// Verify that the algorithm is correct.
// Breaks if a discrete free isometry does not satisfy the Ping Pong Lemma.
count := 0;
while true do
  gens := RandomGens(5,Qp,10);
  count +:= 1;
  t := Cputime();
  result, reducedGens, steps := IsDiscreteFreeCount(T, gens);
  t := Cputime(t);
  print count, result, steps, t/steps;
  if result and not SatisfiesPPL(T, reducedGens) then
    print result;
    break;
  end if;
end while;

// Verify that the algorithm correctly reduces a discrete free group.
count := 0;
while true do
  count +:= 1;
  // Find a random discrete free group X.
  repeat
    gens := RandomGens(5,Qp,10);
    result, reducedGens := IsDiscreteFree(T, gens);
  until result;
  X := RandomAut(reducedGens, 10);
  print count, IsSameGroup(T, reducedGens, X);
end while;

// Time the program against number of generators.
times := 1000;
for j in [2, 3, 5, 10, 20, 50, 100] do
  gensList := [RandomGens(j, Qp, 10) : i in [1..times]];
  nsteps := 0;
  t := Cputime();
  for gens in gensList do
    b, result, steps := IsDiscreteFreeCount(T, gens);
    nsteps +:= steps;
  end for;
  t := Cputime(t);
  print j, t/times, t/nsteps;
end for;

// Time the program against max valuation of entries.
times := 1000;
for j in [3, 5, 7, 10] do
  gensList := [RandomGens(5, Qp, j) : i in [1..times]];
  nsteps := 0;
  t := Cputime();
  for gens in gensList do
    b, result, steps := IsDiscreteFreeCount(T, gens);
    nsteps +:= steps;
  end for;
  t := Cputime(t);
  print j, t/times, t/nsteps;
end for;

// Time the program against precision.
times := 100;
for j in [100..1000 by 100] do
  Qp2 := pAdicField(997,j);
  T2 := BruhatTitsTree(Qp2);
  gensList := [RandomGens(5, Qp2, 3) : i in [1..times]];
  nsteps := 0;
  t := Cputime();
  for gens in gensList do
    b, result, steps := IsDiscreteFreeCount(T2, gens);
    nsteps +:= steps;
  end for;
  t := Cputime(t);
  print j, t/times, t/nsteps;
end for;

// Time the program against number of automorphisms of a discrete free group.
times := 10;
for nAuts in [0..10] do
  gensList := [];
  stepsum := 0;
  for i in [1..times] do
    repeat
      gens := RandomGens(5,Qp,5);
      result, reducedGens := IsDiscreteFree(T, gens);
    until result;
    gensList[i] := RandomAut(reducedGens, nAuts);
  end for;
  t := Cputime();
  for gens in gensList do
    b, result, steps := IsDiscreteFreeCount(T, gens);
    stepsum +:= steps;
  end for;
  print Cputime(t)/times, Sprintf("%.1o", Real(stepsum) / times);
end for;

// Time the average number of steps.
times := 100;
stepsum := 0;
gensList := [RandomGens(5, Qp, 10) : i in [1..times]];
t := Cputime();
for gens in gensList do
  b, result, steps := IsDiscreteFreeCount(T, gens);
  stepsum +:= steps;
end for;
print Cputime(t)/stepsum;
end for;